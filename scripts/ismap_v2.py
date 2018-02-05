#!/usr/bin/env python

import logging
import sys, re, os
from argparse import ArgumentParser
from subprocess import call, check_output, CalledProcessError, STDOUT, Popen, PIPE
import pathlib
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import resource
import time
import read_grouping
from run_commands import run_command, CommandError, BedtoolsError, make_directories
from mapping_to_query import map_to_is_query
from mapping_to_ref import map_to_ref_seq, create_bed_files


def parse_args():
    """
    Parse the input arguments, use -h for help.
    """
    ## TODO:
    # go back and fix these when running on command line
    parser = ArgumentParser(description='IS mapper')

    #parser.add_argument("--version", action='version', version='%(prog)s ' + ismap_version)
    # Inputs
    parser.add_argument('--runtype', type=str, required=False, help='"typing" or "improvement"', default='typing')
    parser.add_argument('--reads', nargs='+', type=pathlib.Path, required=False, help='Paired end reads for analysing (can be gzipped)')
    parser.add_argument('--forward', type=str, required=False, default='_1', help='Identifier for forward reads if not in MiSeq format (default _1)')
    parser.add_argument('--reverse', type=str, required=False, default='_2', help='Identifier for reverse reads if not in MiSeq format (default _2)')
    parser.add_argument('--queries', type=str, required=False, help='Multifasta file for query gene(s) (eg: insertion sequence) that will be mapped to.')
    parser.add_argument('--assemblies', nargs='+', type=str, required=False, help='Contig assemblies, one for each read set')
    parser.add_argument('--assemblyid', type=str, required=False, help='Identifier for assemblies eg: sampleName_contigs (specify _contigs) or sampleName_assembly (specify _assembly). Do not specify extension.')
    parser.add_argument('--extension', type=str, required=False, help='Extension for assemblies (eg: .fasta, .fa, .gbk, default is .fasta)', default='.fasta')
    parser.add_argument('--typingRef', type=str, required=False, help='Reference genome for typing against in genbank format')
    parser.add_argument('--type', type=str, required=False, default='fasta', help='Indicator for contig assembly type, genbank or fasta (default fasta)')
    parser.add_argument('--path', type=str, required=False, default='', help='Path to folder where scripts are (only required for development, default is VLSCI path).')
    # Parameters
    parser.add_argument('--cutoff', type=int, required=False, default=6, help='Minimum depth for mapped region to be kept in bed file (default 6)')
    parser.add_argument('--min_range', type=str, required=False, default='0.2', help='Minimum percent size of the gap to be called a known hit (default 0.2, or 20 percent)')
    parser.add_argument('--max_range', type=str, required=False, default='1.1', help='Maximum percent size of the gap to be called a known hit (default 1.1, or 110 percent)')
    parser.add_argument('--merging', type=str, required=False, default='100', help='Value for merging left and right hits in bed files together to simply calculation of closest and intersecting regions (default 100).')
    parser.add_argument('--a', action='store_true', required=False, help='Switch on all alignment reporting for bwa.')
    parser.add_argument('--T', type=str, required=False, default='30', help='Mapping quality score for bwa (default 30).')
    parser.add_argument('--t', type=str, required=False, default='1', help='Number of threads for bwa (default 1).')
    parser.add_argument('--min_clip', type=int, required=False, default='10', help='Minimum size for softclipped region to be extracted from initial mapping (default 10).')
    parser.add_argument('--max_clip', type=int, required=False, default=30, help='Maximum size for softclipped regions to be included (default 30).')
    # Options for table output (typing)
    parser.add_argument('--cds', nargs='+', type=str, required=False, default=['locus_tag', 'gene', 'product'], help='qualifiers to look for in reference genbank for CDS features (default locus_tag gene product)')
    parser.add_argument('--trna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for tRNA features (default locus_tag product)')
    parser.add_argument('--rrna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for rRNA features (default locus_tag product)')
    parser.add_argument('--igv', action='store_true', help='format of output bedfile - if True, adds IGV trackline and formats 4th column for hovertext display')
    parser.add_argument('--chr_name', type=str, required=False, default='not_specified', help='chromosome name for bedfile - must match genome name to load in IGV (default = genbank accession)')
    # Reporting options
    parser.add_argument('--log', action='store_true', required=False, help='Switch on logging to file (otherwise log to stdout')
    parser.add_argument('--output', type=str, required=False, help='Prefix for output files. If not supplied, prefix will be current date and time.', default='~/Desktop/ismap_v2/')
    parser.add_argument('--temp', action='store_true', required=False, help='Switch on keeping the temp folder instead of deleting it at the end of the program')
    parser.add_argument('--bam', action='store_true', required=False, help='Switch on keeping the final bam files instead of deleting them at the end of the program')
    parser.add_argument('--directory', type=str, required=False, default='', help='Output directory for all output files.')

    return parser.parse_args()


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
    if len(seq_records) == 0:
        #TODO: make this error more informative as to which files contained no sequence
        raise NoSeqError('No sequences were found in one of your input files, please check both your IS queries file(s) and your reference genome file(s).')

    return seq_records

def main():
    # get arguments
    args = parse_args()

    if args.output != '':
        working_dir = os.path.expanduser(args.output)
   # else:
        #working_dir = os.getcwd()

    print(working_dir)

    # set up logfile
    logging.basicConfig(
        filename=working_dir + 'tmp_log.log',
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))

    #TODO: make a function that checks the input files are correct

    #TODO: make a function that checks all dependencies are correct

    # group reads
    #read_groups = read_grouping.group_reads(args.reads)
    test_reads = [pathlib.Path("/Users/jane/Desktop/ismap_v2/reads/9262_1#29_1.fastq.gz"), pathlib.Path("/Users/jane/Desktop/ismap_v2/reads/9262_1#29_2.fastq.gz")]
    read_groups = read_grouping.group_reads(test_reads)

    # print out unpaired reads
    if read_groups.unpaired:
        # report unpaired reads
        logging.info('The following reads were identified as unpaired: ISMapper can only use paired reads.')
        for unpaired in read_groups.unpaired:
            logging.info(unpaired.prefix)
    # print out paired reads
    for paired in read_groups.paired:
        print(paired.prefix)
        forward_read = paired.forward
        reverse_read = paired.reverse

    # parse queries
    query_single = ["/Users/jane/Desktop/ismap_v2/queries/ISAba1.fasta"]
    query_records = get_sequences(query_single, 'fasta')
    #query_list = ["/Users/jane/Desktop/ismap_v2/queries/IS26.fasta", "/Users/jane/Desktop/ismap_v2/queries/ISAba1.fasta"]
    #query_records = get_queries(query_list)
    #query_multi = ["/Users/jane/Desktop/ismap_v2/queries/two_ISqueries.fasta"]
    #query_records = get_queries(query_multi)

    single_ref = ["/Users/jane/Desktop/ismap_v2/refs/CP010781.gbk"]

    tmp_folder = "/Users/jane/Desktop/ismap_v2/test_results/tmp"

    # read in the reference genomes to type against
    reference_seqs = get_sequences(single_ref, 'genbank')

    for sample in read_groups.paired:

        # set up output folder for each strain:
        output_sample = os.path.join(working_dir, sample.prefix)
        #make_directories(output_sample)

        # regardless of which mode we're in, we need to do the following, for each query:
        for is_query in query_records:

            left_flanking_reads, right_flanking_reads = map_to_is_query(sample, is_query, output_sample)

            # typing mode first
            if args.runtype == 'typing':
                # we need to loop through each reference:
                for ref_seq in reference_seqs:
                    # map our flanking reads to this
                    map_to_ref_seq(ref_seq, left_flanking_reads, right_flanking_reads)

                    # make the bed files
                    create_bed_files()


                # Set up file names for output files


                # if it's a genbank, convert to fasta
                # can just be a fasta also - we just won't get gene annotations

                # Map reads to reference, sort

                # Create BED files with coverage information

                # Filter coveraged BED files on coverage cutoff (so only take
                # high coverage regions for further analysis)

                # Find intersects and closest points of regions

                # Create all possible closest bed files for checking unpaired hits
                # If any of these fail, just make empty unapired files to pass to create_typing_out

                # Create table and annotate genbank with hits
                pass



if __name__ == '__main__':
    main()