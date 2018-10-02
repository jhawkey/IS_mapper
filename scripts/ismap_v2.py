#!/usr/bin/env python

import logging
import sys, os
import argparse
import pathlib
import datetime
from Bio import SeqIO
import time
import read_grouping
from run_commands import run_command, CommandError, BedtoolsError, make_directories, check_command, remove_files
from mapping_to_query import map_to_is_query
from mapping_to_ref import map_to_ref_seq, create_bed_files
from create_output import create_typing_output


def parse_args():
    """
    Parse the input arguments, use -h for help.
    """
    # go back and fix these when running on command line
    parser_parent = argparse.ArgumentParser()
    parser_script = parser_parent.add_argument_group('Basic ISMapper options')
    parser_hits = parser_parent.add_argument_group('Parameters for defining hits')
    parser_bwa = parser_parent.add_argument_group('BWA parameters')
    parser_table = parser_parent.add_argument_group('Parameters for output table')
    parser_output = parser_parent.add_argument_group('Reporting parameters')

    #parser_parent.add_argument("--version", action='version', version='%(prog)s ' + ismap_version)
    # Inputs
    parser_script.add_argument('--reads', nargs='+', type=pathlib.Path, required=True,
                               help='Paired end reads for analysing (can be gzipped)')
    parser_script.add_argument('--queries', type=str, nargs='+', required=True,
                               help='Multifasta file for query gene(s) (eg: insertion sequence) that will be mapped to.')
    parser_script.add_argument('--reference', type=str, nargs='+', required=True,
                               help='Reference genome for typing against in genbank format')
    parser_script.add_argument('--output', type=str, required=False, default=os.getcwd(),
                               help='Location for all output files (default is current directory).')
    parser_script.add_argument('--log', action='store_true', required=False,
                               help='Switch on logging to file (otherwise log to stdout)')
    parser_script.add_argument('--log_name', type=str, required=False,
                                 help='Prefix for log file. If not supplied, prefix will be current date and time.',
                                 default=datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))

    parser_script.add_argument('--help_all', required=False, action=None, help='Display extended help')

    # List of arguments to show in quick help
    quick_help_args = ('--reads', '--queries', '--reference', '--output', '--log', '--log_name')

    # Parameters for defining hits
    parser_hits.add_argument('--min_clip', type=int, required=False, default=10,
                                 help='Minimum size for softclipped region to be extracted from initial mapping (default 10).')
    parser_hits.add_argument('--max_clip', type=int, required=False, default=30,
                                 help='Maximum size for softclipped regions to be included (default 30).')
    parser_hits.add_argument('--cutoff', type=int, required=False, default=6,
                                 help='Minimum depth for mapped region to be kept in bed file (default 6)')
    parser_hits.add_argument('--min_range', type=str, required=False, default=0.9,
                                 help='Minimum percent size of the gap to be called a known hit (default 0.9, or 90 percent)')
    parser_hits.add_argument('--max_range', type=str, required=False, default=1.1,
                                 help='Maximum percent size of the gap to be called a known hit (default 1.1, or 110 percent)')
    parser_hits.add_argument('--merging', type=str, required=False, default='100',
                                 help='Value for merging left and right hits in bed files together to simply calculation of closest and intersecting regions (default 100).')

    # BWA/read grouping parameters
    parser_bwa.add_argument('--a', action='store_true', required=False,
                                 help='Switch on all alignment reporting for bwa.')
    parser_bwa.add_argument('--T', type=str, required=False, default='30',
                                 help='Mapping quality score for bwa (default 30).')
    parser_bwa.add_argument('--t', type=str, required=False, default='1',
                                 help='Number of threads for bwa (default 1).')
    parser_bwa.add_argument('--forward', type=str, required=False, default='_1',
                                 help='Identifier for forward reads if ISMapper is unable to pair (default is Miseq format _1)')
    parser_bwa.add_argument('--reverse', type=str, required=False, default='_2',
                                 help='Identifier for forward reads if ISMapper is unable to pair (default is Miseq format _2)')

    # Output table parameters
    parser_table.add_argument('--cds', nargs='+', type=str, required=False, default=['locus_tag', 'gene', 'product'],
                                 help='qualifiers to look for in reference genbank for CDS features (default locus_tag gene product)')
    parser_table.add_argument('--trna', nargs='+', type=str, required=False, default=['locus_tag', 'product'],
                                 help='qualifiers to look for in reference genbank for tRNA features (default locus_tag product)')
    parser_table.add_argument('--rrna', nargs='+', type=str, required=False, default=['locus_tag', 'product'],
                                 help='qualifiers to look for in reference genbank for rRNA features (default locus_tag product)')
    parser_table.add_argument('--igv', action='store_true',
                                 help='format of output bedfile - if True, adds IGV trackline and formats 4th column for hovertext display')
    parser_table.add_argument('--chr_name', type=str, required=False, default='not_specified',
                                 help='chromosome name for bedfile - must match genome name to load in IGV (default = genbank accession)')

    # Reporting parameters
    parser_output.add_argument('--temp', action='store_true', required=False,
                               help='Switch on keeping the temp folder instead of deleting it at the end of the run')
    parser_output.add_argument('--bam', action='store_true', required=False,
                               help='Switch on keeping the final bam files instead of deleting them at the end of the run')

    # Suppress all arguments other than quick help unless --help_all is parsed
    if '--help_all' in sys.argv[1:]:
        parser_parent.print_help()
        sys.exit(0)
    else:
        for arg in parser_parent._actions:
            if not any(qarg in arg.option_strings for qarg in quick_help_args):
                arg.help = argparse.SUPPRESS

    args = parser_parent.parse_args()

    return args


class NoSeqError(Exception):
    def __init__(self, message):
        self.message = message

def get_sequences(seq_files, seq_format):
    """
    Takes a list of sequence files and their respective format (all must be the same format).
    Returns a list of SeqRecord objects for downstream use.
    """

    seq_records = []
    for seq in seq_files:
        seq_parse = SeqIO.parse(seq, seq_format)
        for record in seq_parse:
            seq_records.append(record)

    # if the list is empty, raise an error
    if len(seq_records) == 0 and seq_format == 'fasta':
        logging.error('One of your IS queries contained no sequence, please check')
        logging.info('ISMapper exiting')
        raise NoSeqError('No sequence was found in one of your IS query file(s).')

    # if the list is empty, raise an error
    if len(seq_records) == 0 and seq_format == 'genbank':
        logging.error('No entry was found in your reference genbank(s), please check')
        logging.info('ISMapper exiting')
        raise NoSeqError('No entry was found in your reference genbank file(s).')

    return seq_records

def main():

    # intialise the start time
    start_time = time.time()

    # get arguments
    args = parse_args()

    if args.output != '':
        working_dir = os.path.expanduser(args.output)
    else:
        working_dir = os.getcwd()

    # TODO: don't set up logfile if log isn't set to true - default should be logfile, turning off logfile prints to stdout
    # set up logfile
    logging.basicConfig(
        filename=args.log_name + '.log',
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%d/%m/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))
    logging.info(working_dir)

    # Checks that the correct programs are installed
    check_command('bwa', 'BWA')
    check_command('samtools', 'SAMtools')
    check_command('makeblastdb', 'BLAST')
    check_command('bedtools', 'BedTools')

    # group reads
    read_groups = read_grouping.group_reads(args.reads)

    # print out unpaired reads
    if read_groups.unpaired:
        # report unpaired reads
        logging.info('ISMapper can only use paired reads. The following reads were identified as unpaired:')
        for unpaired in read_groups.unpaired:
            logging.info(unpaired.prefix)
    # print out paired reads
    if read_groups.paired:
        if len(read_groups.paired) == 1:
            logging.info('Found %s set of paired reads', len(read_groups.paired))
        else:
            logging.info('Found %s sets of paired reads', len(read_groups.paired))
        for paired in read_groups.paired:
            logging.info('Found paired reads for sample %s: %s %s', paired.prefix, str(paired.forward.filepath), str(paired.reverse.filepath))
            # check none of the input read files are empty
            # if they are, raise an error
            if os.stat(str(paired.forward.filepath))[6] <= 100:
                logging.error('Forward read %s is empty! ISMapper exiting', str(paired.forward.filepath))
                raise NoSeqError('Forward read file is empty, ISMapper exiting')
            if os.stat(str(paired.reverse.filepath))[6] <= 100:
                logging.error('Reverse read %s is empty! ISMapper exiting', str(paired.forward.filepath))
                raise NoSeqError('Reverse read file is empty, ISMapper exiting')
    # if there were no paired read sets found, raise an exception and quit
    else:
        logging.error('No paired read sets found, ISMapper exiting')
        raise NoSeqError('No paired read sets found, ISMapper exiting')

    # parse queries
    query_records = get_sequences(args.queries, 'fasta')
    for record in query_records:
        logging.info('Found query %s', record.id)

    # read in the reference genomes to type against
    reference_seqs = get_sequences(args.reference, 'genbank')
    for record in reference_seqs:
        logging.info('Found reference sequence %s', record.id)

    for sample in read_groups.paired:

        logging.info('Processing sample %s', sample.prefix)
        # set up output folder for each strain:
        output_sample = os.path.join(working_dir, sample.prefix)

        # we need to do the following, for each query:
        for is_query in query_records:
            logging.info('Processing IS query %s against sample %s', is_query.id, sample.prefix)
            left_flanking_reads, right_flanking_reads, is_output_folder, tmp_output_folder = map_to_is_query(sample, is_query, output_sample, args.min_clip, args.max_clip, args.t)


            # we need to loop through each reference:
            for ref_seq in reference_seqs:
                # map our flanking reads to this
                filenames = map_to_ref_seq(ref_seq, sample.prefix, left_flanking_reads, right_flanking_reads, tmp_output_folder, is_output_folder, args.t, args.a)

                # make the bed files, find intersects and closest points of regions
                filenames_bedfiles = create_bed_files(filenames, args.cutoff, args.merging)

                # Create table and annotated genbank with hits
                create_typing_output(filenames_bedfiles, ref_seq, is_query, args.min_range, args.max_range, tmp_output_folder, sample.prefix)
                logging.info('ISMapper has completed successfully for sample %s', sample.prefix)

                # TODO: remove temp directory unless specifically asked not to
                if not args.temp:
                    remove_files([tmp_output_folder])
                if not args.bam:
                    remove_files([filenames['left_sorted'], filenames['right_sorted']])

    total_time = time.time() - start_time
    time_mins = float(total_time) / 60
    logging.info('ISMapper finished in ' + str(time_mins) + ' mins.')


if __name__ == '__main__':
    main()