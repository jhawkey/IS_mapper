#!/usr/bin/env python

# ISMapper
# Python Version 2.7.5
#
# Authors - Jane Hawkey (jhawkey@student.unimelb.edu.au), Kathryn Holt (kholt@unimelb.edu.au)
#
# see LICENSE.txt for the license
#
# Dependencies:
#   BWA v0.7.5a - http://bio-bwa.sourceforge.net/
#   Samtools v0.1.19 - http://samtools.sourceforge.net/
#   Bedtools v2.20.1 - http://bedtools.readthedocs.org/en/latest/content/installation.html
#   BioPython v1.63 - http://biopython.org/wiki/Main_Page
#   BLAST+ v2.2.28 - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/
#
# Git repository: https://github.com/jhawkey/IS_mapper
# README: https://github.com/jhawkey/IS_mapper/blob/master/README.txt
# Questions or feature requests: https://github.com/jhawkey/IS_mapper/issues

import logging
import sys, re, os
from argparse import ArgumentParser
from subprocess import call, check_output, CalledProcessError, STDOUT
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import resource
import time
try:
    from version import ismap_version
except:
    ismap_version = 'version unknown'

def parse_args():
    '''
    Parse the input arguments, use -h for help.
    '''

    parser = ArgumentParser(description='IS mapper')

    parser.add_argument("--version", action='version', version='v1.2')
    # Inputs
    parser.add_argument('--runtype', type=str, required=True, help='"typing" or "improvement"')
    parser.add_argument('--reads', nargs='+', type=str, required=False, help='Paired end reads for analysing (can be gzipped)')
    parser.add_argument('--forward', type=str, required=False, default='_1', help='Identifier for forward reads if not in MiSeq format (default _1)')
    parser.add_argument('--reverse', type=str, required=False, default='_2', help='Identifier for reverse reads if not in MiSeq format (default _2)')
    parser.add_argument('--queries', nargs='+', type=str, required=True, help='Fasta files for query genes (eg: insertion sequence) that will be mapped to')
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
    parser.add_argument('--output', type=str, required=True, help='prefix for output files')
    parser.add_argument('--temp', action='store_true', required=False, help='Switch on keeping the temp folder instead of deleting it at the end of the program')
    parser.add_argument('--bam', action='store_true', required=False, help='Switch on keeping the final bam files instead of deleting them at the end of the program')
    parser.add_argument('--directory', type=str, required=False, default='', help='Output directory for all output files.')

    return parser.parse_args()

# Exception to raise if the command we try to run fails for some reason
class CommandError(Exception):
    pass

def run_command(command, **kwargs):
    '''
    Execute a shell command and check the exit status and any O/S exceptions.
    '''

    command_str = ' '.join(command)
    logging.info('Running: {}'.format(command_str))
    try:
        exit_status = call(command_str, **kwargs)
    except OSError as e:
        message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
        raise CommandError({"message": message})
    if exit_status != 0:
        message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
        raise CommandError({"message": message})

def bwa_index(fasta):
    '''
    Check to see if bwa index for given input fasta exists.
    If it doesn't, build an index from the given input fasta.
    '''

    #check_command_version('bwa')
    built_index = fasta + '.bwt'
    print built_index
    if os.path.exists(built_index):
        logging.info('Index for {} is already built...'.format(fasta))
    else:
        logging.info('Building bwa index for {}...'.format(fasta))
        run_command(['bwa', 'index', fasta], shell=True)

def check_command_version(command_list, version_identifier, command_name, required_version):
    '''
    Check that an acceptable version of a command is installed.
    Exits the program if it can't be found.
        - command_list is the command to run to determine the version.
        - version_identifier is the unique string we look for in the stdout of the program.
        - command_name is the name of the command to show in error messages.
        - required_version is the version number to show in error messages.
    '''

    try:
        command_stdout = check_output(command_list, stderr=STDOUT)
    except OSError as e:
        logging.error("Failed command: {}".format(' '.join(command_list)))
        logging.error(str(e))
        logging.error("Could not determine the version of {}.".format(command_name))
        logging.error("Do you have {} installed in your PATH?".format(command_name))
        exit(-1)
    except CalledProcessError as e:
        # some programs such as samtools return a non-zero exit status
        # when you ask for the version (sigh). We ignore it here.
        command_stdout = e.output

    if version_identifier not in command_stdout:
        logging.error("Incorrect version of {} installed.".format(command_name))
        logging.error("{} version {} is required by ISmapper.".format(command_name, required_version))
        exit(-1)

def check_command(command_list, command_name):
    '''
    Check that the dependency is installed.
    Exits the program if it can't be found.
        - command_list is the command to run to determine the version.
        - command_name is the name of the command to show in the error message.
    '''
    try:
        command_stdout = check_output(command_list, stderr=STDOUT)
    except OSError as e:
        logging.error("Failed command: {}".format(' '.join(command_list)))
        logging.error(str(e))
        logging.error("Do you have {} installed in your PATH?".format(command_name))
        exit(-1)
    except CalledProcessError as e:
        # some programs such as samtools return a non zero exit status
        # when you ask for the version. We ignore it here.
        command_stdout = e.output

def get_readFile_components(full_file_path):
    '''
    Takes the path to the read file and splits it into
    its different parts.
    Returns the file path, the file name and the file extension.
    '''

    (file_path, file_name) = os.path.split(full_file_path)
    m1 = re.match('(.*).gz', file_name)
    ext = ''
    if m1 != None:
        # gzipped
        ext = '.gz'
        file_name = m1.groups()[0]
    (file_name_before_ext, ext2) = os.path.splitext(file_name)
    full_ext = ext2+ext

    return(file_path, file_name_before_ext, full_ext)

def read_file_sets(args):
    '''
    Takes the read files and pairs them together.
    If the improvement pathway is selected, also finds their
    respecitive assemblies puts each set together.
    Returns a dictionary where the key is the id fo the sample, and the value
    is a list of files that belong to that sample.
    '''

    fileSets = {} # key = id, value = list of files for that sample
    num_paired_readsets = 0

    # paired end
    forward_reads = {} # key = sample, value = full path to file
    reverse_reads = {} # key = sample, value = full path to file
    assemblies = {}
    num_paired_readsets = 0
    num_single_readsets = 0
    num_assemblies = 0

    for fastq in args.reads:
        (file_path, file_name_before_ext, full_ext) = get_readFile_components(fastq)
        # try to match to MiSeq format:
        m = re.match('(.*)(_S.*)(_L.*)(_R.*)(_.*)', file_name_before_ext)
        if m == None:
            # not default Illumina file naming format, expect simple/ENA format
            m = re.match('(.*)(' + args.forward + ')$', file_name_before_ext)
            if m != None:
                # store as forward read
                (baseName, read) = m.groups()
                forward_reads[baseName] = fastq
            else:
                m = re.match('(.*)(' + args.reverse + ')$', file_name_before_ext)
                if m != None:
                # store as reverse read
                    (baseName, read) = m.groups()
                    reverse_reads[baseName] = fastq
                else:
                    print 'Could not determine forward/reverse read status for input file ' + fastq
        else:
            # matches default Illumina file naming format, e.g. m.groups() = ('samplename', '_S1', '_L001', '_R1', '_001')
            baseName, read = m.groups()[0], m.groups()[3]
            if read == '_R1':
                forward_reads[baseName] = fastq
            elif read == '_R2':
                reverse_reads[baseName] = fastq
            else:
                print 'Could not determine forward/reverse read status for input file ' + fastq
                print '  this file appears to match the MiSeq file naming convention (samplename_S1_L001_[R1]_001), but we were expecting [R1] or [R2] to designate read as forward or reverse?'
                fileSets[file_name_before_ext] = fastq
                num_single_readsets += 1

    if args.runtype == 'improvement':
        for assembly in args.assemblies:
            (file_path, file_name_before_ext, full_ext) = get_readFile_components(assembly)
            identifier = file_name_before_ext.split(args.assemblyid)[0]
            assemblies[identifier] = assembly

    # store in pairs with assemblies
    for sample in forward_reads:
        if sample in reverse_reads and sample in assemblies:
            fileSets[sample] = [forward_reads[sample],reverse_reads[sample], assemblies[sample]] # store pair and assembly
            num_paired_readsets += 1
            num_assemblies += 1
        elif sample in reverse_reads:
            fileSets[sample] = [forward_reads[sample], reverse_reads[sample]] #store pair
            num_paired_readsets += 1
        else:
            fileSets[sample] = [forward_reads[sample]] # no reverse found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + forward_reads[sample])
    for sample in reverse_reads:
        if sample not in fileSets:
            fileSets[sample] = reverse_reads[sample] # no forward found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + reverse_reads[sample])
    for sample in assemblies:
        if sample not in fileSets:
            fileSets[sample] = assemblies[sample]
            num_assemblies += 1
            logging.info('Warning, could not find reads for assembly:' + assemblies[sample])

    if num_paired_readsets > 0:
        logging.info('Total paired readsets found:' + str(num_paired_readsets))
    if num_single_readsets > 0:
        logging.info('Total single reads found:' + str(num_single_readsets))
    if num_assemblies > 0:
        logging.info('Total number of assemblies found:' + str(num_assemblies))

    return fileSets

def make_directories(dir_list):
    '''
    Makes the directories specified in the list.
    Checks to make sure they exist first.
    '''
    for directory in dir_list:
        run_command(['mkdir', '-p', directory], shell=True)

def filter_on_depth(cov_file, out_bed, cov_cutoff):
    '''
    Takes a bed coverage file and removes lines that
    do not meet the coverage cutoff.
    Saves output to a new file.
    '''

    output = file(out_bed, 'w')
    with open(cov_file) as depth_info:
        for line in depth_info:
            if int(line.strip().split('\t')[3]) >= cov_cutoff:
                output.write(line)
    output.close()

def check_blast_database(fasta):
    '''
    Checks to make sure the BLAST database exists, and creates it
    if it does not.
    '''

    database_path = fasta + ".nin"

    if os.path.exists(database_path):
        logging.info('Index for {} is already built...'.format(fasta))
    else:
        logging.info('Building blast index for {}...'.format(fasta))
        run_command(['makeblastdb -in', fasta, '-dbtype nucl'], shell=True)

def gbk_to_fasta(genbank, fasta):
    '''
    Converts a genbank to a fasta using BioPython
    '''

    sequences = SeqIO.parse(genbank, "genbank")
    SeqIO.write(sequences, fasta, "fasta")

def multi_to_single(genbank, name, output):
    '''
    Converts a multi entry genbank (where each entry is a contig)
    into a single entry genbank, preserving all annotations.
    '''

    # total bases
    total = 0

    handle = open(genbank, "rU")
    records = list(SeqIO.parse(handle, "genbank"))
    feature_count = 0
    colour_count = 0

    # make header genbank format friendly
    if len(name) >= 10:
        name = name[:9]
    for r in records:
        length = len(r)
        id = r.name
        seq = r.seq
        seq.alphabet=generic_dna
        if total > 0:
            newrecord.seq = newrecord.seq + seq
        else:
            # first sequence, initialise seqrecord
            newrecord = SeqRecord(seq=r.seq,name=name,id=name)
            newrecord.seq.alphabet=generic_dna
        # create feature for contig
        if colour_count % 2 == 0:
            newrecord.features.append(SeqFeature(FeatureLocation(total, total + length), type="fasta_record", qualifiers = {'note' : [r.name], 'colour':'11'}))
            colour_count = colour_count + 1
        else:
            newrecord.features.append(SeqFeature(FeatureLocation(total, total + length), type="fasta_record", qualifiers = {'note' : [r.name], 'colour':'10'}))
            colour_count = colour_count + 1
        # copy CDS features
        for f in r.features:
            feature_count += 1
            f.qualifiers["locus_tag"] = str(feature_count)
            newrecord.features.append(SeqFeature(FeatureLocation(f.location.nofuzzy_start + total, f.location.nofuzzy_end + total), strand = f.strand, type=f.type, qualifiers = f.qualifiers))
        total += length
    handle.close()
    #write out new single entry genbank
    SeqIO.write(newrecord, output, "genbank")

def extract_clipped_reads(sam_file, min_size, max_size, out_five_file, out_three_file):
    with open(sam_file, 'r') as in_file, open(out_five_file, 'w') as out_five, open(out_three_file, 'w') as out_three:
        print "extracting clip reads into" + out_five_file + out_three_file
        for line in in_file:
            # split fields
            entries = line.split('\t')
            # check if it's a header line, and if so, skip it
            if re.search('^@[A-Z][A-Z]$', entries[0]):
                continue
            # check the SAM flag - if it's unmapped, skip this read.
            sam_flag = int(entries [1])
            if sam_flag & 4:
                continue
            # check if the read has been reverse complemented
            reverse_complement = sam_flag & 16
            #grab the read name and cigar string
            read_name, cigar = entries[0], entries[5]
            #parse the cigar info, exclude hard-clipped regions (these can only be on the very edges, outside of soft-clips if soft-clipping is present)
            map_regions = re.findall('[0-9]+[MIDNSP=X]', cigar)

            # Get the first and last items from the cigar string and see if it's soft-clipped (last letter = S). Soft clips will only ever be at the ends (inside would be I/D/N)
            # If so, find out how many bases are soft-clipped
            # If it's the right size, add the read and it's quality score to the appropriate fastq file, reverse-complementing if needed
            if map_regions[0][-1] == 'S':
                num_soft_clipped = int(map_regions[0][:-1])
                if min_size <= num_soft_clipped <= max_size:
                    soft_clipped_seq = Seq(entries[9][:num_soft_clipped], generic_dna)
                    qual_scores = entries[10][:num_soft_clipped]
                    if reverse_complement:
                        out_five.write('@' + read_name + '\n' + str(soft_clipped_seq.reverse_complement()) + '\n+\n' + qual_scores[::-1] + '\n')
                    else:
                        out_five.write('@' + read_name + '\n' + str(soft_clipped_seq) + '\n+\n' + qual_scores + '\n')
            if map_regions[-1][-1] == 'S':
                num_soft_clipped = int(map_regions[-1][:-1])
                if min_size <= num_soft_clipped <= max_size:
                    soft_clipped_seq = Seq(entries[9][-num_soft_clipped:], generic_dna)
                    qual_scores = entries[10][-num_soft_clipped:]
                    if reverse_complement:
                        out_three.write('@' + read_name + '\n' + str(soft_clipped_seq.reverse_complement()) + '\n+\n' + qual_scores[::-1] + '\n')
                    else:
                        out_three.write('@' + read_name + '\n' + str(soft_clipped_seq) + '\n+\n' + qual_scores + '\n')

def main():

    start_time = time.time()

    args = parse_args()

    # Set up logfile
    if args.log is True:
        logfile = os.path.join(args.directory, args.output + ".log")
    else:
        logfile = None
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))

    # Checks that the correct programs are installed
    check_command(['bwa'], 'bwa')
    check_command(['samtools'], 'samtools')
    check_command(['makeblastdb'], 'blast')
    check_command(['bedtools'], 'bedtools')

    # Checks to make sure the runtype is valid and provides an error if not
    if args.runtype != "improvement" and args.runtype != "typing":
        logging.info('Invalid runtype selected: {}'.format(args.runtype))
        logging.info('Runtype should be improvement or typing (see instructions for further details)')
        exit(-1)

    # Get feature types in correct format
    args.cds = ' '.join(args.cds)
    args.trna = ' '.join(args.trna)
    args.rrna = ' '.join(args.rrna)

    # Gather together the reads in pairs with their corresponding
    # assemblies (if required)
    fileSets = read_file_sets(args)
    # Start analysing each read set specified
    for sample in fileSets:
        forward_read = fileSets[sample][0]
        reverse_read = fileSets[sample][1]
        try:
            assembly = fileSets[sample][2]
        except IndexError:
            pass

        # Cycle through each query on its own before moving onto the next one
        for query in args.queries:

            # Index the IS query for BWA
            bwa_index(query)
            # Get query name
            query_name = os.path.split(query)[1].split('.f')[0]

            # Create the output file and folder names,
            # make the folders where necessary
            if args.directory == '':
                current_dir = os.getcwd() + '/'
            else:
                current_dir = args.directory
            if current_dir[-1] != '/':
                current_dir = current_dir + '/'

            temp_folder = current_dir + sample + '_' + query_name + '_temp/'
            output_sam = temp_folder + sample + '_' + query_name + '.sam'
            five_bam = temp_folder + sample + '_' + query_name + '_5.bam'
            three_bam = temp_folder + sample + '_' + query_name + '_3.bam'
            five_reads = temp_folder + sample + '_' + query_name + '_5.fastq'
            three_reads = temp_folder + sample + '_' + query_name + '_3.fastq'
            left_clipped_reads = temp_folder + sample + '_' + query_name + '_left_clipped.fastq'
            right_clipped_reads = temp_folder + sample + '_' + query_name + '_right_clipped.fastq'
            final_left_reads = temp_folder + sample + '_' + query_name + '_LeftFinal.fastq'
            final_right_reads = temp_folder + sample + '_' + query_name + '_RightFinal.fastq'
            no_hits_table = current_dir + sample + '_' + query_name + '_table.txt'
            make_directories([temp_folder])

            # Map to IS query
            run_command(['bwa', 'mem', '-t', args.t, query, forward_read, reverse_read, '>', output_sam], shell=True)
            # Pull unmapped reads flanking IS
            run_command(['samtools view', '-Sb', '-f 36', output_sam, '>', five_bam], shell=True)
            run_command(['samtools view', '-Sb', '-f 4', '-F 40', output_sam, '>', three_bam], shell=True)
            # Turn bams to reads for mapping
            run_command(['bedtools', 'bamtofastq', '-i', five_bam, '-fq', five_reads], shell=True)
            run_command(['bedtools', 'bamtofastq', '-i', three_bam, '-fq', three_reads], shell=True)
            # Add corresponding clipped reads to their respective left and right ends
            print 'Usage before extracting soft-clipped reads'
            print ('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            logging.info('Extracting soft clipped reads, selecting reads that are <= ' + str(args.max_clip) + 'bp and >= ' + str(args.min_clip) + 'bp')
            extract_clipped_reads(output_sam, args.min_clip, args.max_clip, left_clipped_reads, right_clipped_reads)
            print 'Usage after reads written out, before concatentation'
            print ('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            run_command(['cat', left_clipped_reads, five_reads, '>', final_left_reads], shell=True)
            run_command(['cat', right_clipped_reads, three_reads, '>', final_right_reads], shell=True)
            print 'Usage after reads concatenated onto previous reads'
            print ('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

            # Create BLAST database for IS query
            check_blast_database(query)
            if os.stat(final_left_reads)[6] == 0 or os.stat(final_right_reads)[6] == 0:
                logging.info('One or both read files are empty. This is probably due to no copies of the IS of interest being present in this sample. Program quitting.')
                with open(no_hits_table, 'w') as f:
                    if args.runtype == 'typing':
                        header = ["region", "orientation", "x", "y", "gap", "call", "%ID", "%Cov", "left_gene", "left_strand", "left_distance", "right_gene", "right_strand", "right_distance", "functional_prediction"]
                        f.write('\t'.join(header) + '\nNo hits found')
                    else:
                        header = ['contig', 'end', 'x', 'y']
                        f.write('\t'.join(header) + '\nNo hits found')
                continue

            # Improvement mode
            if args.runtype == "improvement":

                # Get prefix for output filenames
                five_header = sample + '_5'
                three_header = sample + '_3'
                five_to_ref_sam = temp_folder + five_header + '_' + query_name + '.sam'
                three_to_ref_sam = temp_folder + three_header + '_' + query_name + '.sam'
                five_to_ref_bam = temp_folder + five_header + '_' + query_name + '.bam'
                three_to_ref_bam = temp_folder + three_header + '_' + query_name + '.bam'
                five_bam_sorted = five_header + '_' + query_name + '.sorted'
                three_bam_sorted = three_header + '_' + query_name + '.sorted'
                five_cov_bed = temp_folder + five_header + '_' + query_name + '_cov.bed'
                three_cov_bed = temp_folder + three_header + '_' + query_name + '_cov.bed'
                five_final_cov =  current_dir + five_header + '_' + query_name + '_finalcov.bed'
                three_final_cov = current_dir + three_header + '_' + query_name + '_finalcov.bed'
                five_merged_bed = current_dir + five_header + '_' + query_name + '_merged.sorted.bed'
                three_merged_bed = current_dir + three_header + '_' + query_name + '_merged.sorted.bed'
                final_genbankSingle = current_dir + sample + '_' + query_name + '_annotatedSingle.gbk'

                # create fasta file from genbank if required
                if args.extension == '.gbk':
                    assembly_gbk = assembly
                    (file_path, file_name_before_ext, full_ext) = get_readFile_components(assembly_gbk)
                    assembly_fasta = os.path.join(temp_folder, file_name_before_ext) + '.fasta'
                    gbk_to_fasta(assembly, assembly_fasta)
                    assembly = assembly_fasta
                # Map ends back to contigs
                bwa_index(assembly)
                if args.a == True:
                    run_command(['bwa', 'mem', 'a', '-T', args.T, '-t', args.t, assembly, final_left_reads, '>', five_to_ref_sam], shell=True)
                    run_command(['bwa', 'mem', 'a', '-T', args.T, '-t', args.t, assembly, final_right_reads, '>', three_to_ref_sam], shell=True)
                else:
                    run_command(['bwa', 'mem', '-t', args.t, assembly, final_left_reads, '>', five_to_ref_sam], shell=True)
                    run_command(['bwa', 'mem', '-t', args.t, assembly, final_right_reads, '>', three_to_ref_sam], shell=True)

                run_command(['samtools', 'view', '-Sb', five_to_ref_sam, '>', five_to_ref_bam], shell=True)
                run_command(['samtools', 'view', '-Sb', three_to_ref_sam, '>', three_to_ref_bam], shell=True)
                run_command(['samtools', 'sort', five_to_ref_bam, five_bam_sorted], shell=True)
                run_command(['samtools', 'sort', three_to_ref_bam, three_bam_sorted], shell=True)
                run_command(['samtools', 'index', five_bam_sorted + '.bam'], shell=True)
                run_command(['samtools', 'index', three_bam_sorted + '.bam'], shell=True)
                # Create BED file with coverage information
                run_command(['bedtools', 'genomecov', '-ibam', five_bam_sorted + '.bam', '-bg', '>', five_cov_bed], shell=True)
                run_command(['bedtools', 'genomecov', '-ibam', three_bam_sorted + '.bam', '-bg', '>', three_cov_bed], shell=True)
                filter_on_depth(five_cov_bed, five_final_cov, args.cutoff)
                filter_on_depth(three_cov_bed, three_final_cov, args.cutoff)
                run_command(['bedtools', 'merge', '-i', five_final_cov, '-d', args.merging, '>', five_merged_bed], shell=True)
                run_command(['bedtools', 'merge', '-i', three_final_cov, '-d', args.merging, '>', three_merged_bed], shell=True)
                # Create table and genbank
                if args.extension == '.fasta':
                    run_command([args.path + 'create_genbank_table.py', '--five_bed', five_merged_bed, '--three_bed', three_merged_bed, '--assembly', assembly, '--type fasta', '--output', current_dir + sample + '_' + query_name], shell=True)
                elif args.extension == '.gbk':
                    run_command([args.path + 'create_genbank_table.py', '--five_bed', five_merged_bed, '--three_bed', three_merged_bed, '--assembly', assembly_gbk, '--type genbank', '--output', current_dir + sample + '_' + query_name], shell=True)
                #create single entry genbank
                multi_to_single(sample + '_' + query_name + '_annotated.gbk', sample, final_genbankSingle)

            # Typing mode
            if args.runtype == "typing":

                # Get prefix of typing reference for output filenames
                (file_path, file_name) = os.path.split(args.typingRef)
                typingName = file_name.split('.g')[0]
                typingRefFasta = temp_folder + typingName + '.fasta'
                # Create reference fasta from genbank
                gbk_to_fasta(args.typingRef, typingRefFasta)
                # Create bwa index file for typing reference
                bwa_index(typingRefFasta)
                # Set up file names for output files
                five_header = sample + '_5_' + typingName
                three_header = sample + '_3_' + typingName
                five_to_ref_sam = temp_folder + five_header + '_' + query_name + '.sam'
                three_to_ref_sam = temp_folder + three_header + '_' + query_name + '.sam'
                five_to_ref_bam = temp_folder + five_header + '_' + query_name + '.bam'
                three_to_ref_bam = temp_folder + three_header + '_' + query_name + '.bam'
                five_bam_sorted = five_header + '_' + query_name + '.sorted'
                three_bam_sorted = three_header + '_' + query_name + '.sorted'
                five_cov_bed = temp_folder + five_header + '_' + query_name + '_cov.bed'
                three_cov_bed = temp_folder + three_header + '_' + query_name + '_cov.bed'
                five_cov_merged = temp_folder + five_header + '_' + query_name + '_cov_merged.sorted.bed'
                three_cov_merged = temp_folder + three_header + '_' + query_name + '_cov_merged.sorted.bed'
                five_final_cov = current_dir + five_header + '_' + query_name + '_finalcov.bed'
                three_final_cov = current_dir + three_header + '_' + query_name + '_finalcov.bed'
                five_merged_bed = current_dir + five_header + '_' + query_name + '_merged.sorted.bed'
                three_merged_bed = current_dir + three_header + '_' + query_name + '_merged.sorted.bed'
                bed_intersect = current_dir + sample + '_' + typingName + '_' + query_name + '_intersect.bed'
                bed_closest = current_dir + sample + '_' + typingName + '_' + query_name + '_closest.bed'
                bed_unpaired_five = current_dir + sample + '_' + typingName + '_' + query_name + '_left_unpaired.bed'
                bed_unpaired_three = current_dir + sample + '_' + typingName + '_' + query_name + '_right_unpaired.bed'

                # Map reads to reference, sort
                if args.a == True:
                    run_command(['bwa', 'mem', '-a', '-T', args.T, '-t', args.t, typingRefFasta, final_left_reads, '>', five_to_ref_sam], shell=True)
                    run_command(['bwa', 'mem', '-a', '-T', args.T, '-t', args.t,typingRefFasta, final_right_reads, '>', three_to_ref_sam], shell=True)
                else:
                    run_command(['bwa', 'mem', '-t', args.t, typingRefFasta, final_left_reads, '>', five_to_ref_sam], shell=True)
                    run_command(['bwa', 'mem', '-t', args.t, typingRefFasta, final_right_reads, '>', three_to_ref_sam], shell=True)
                run_command(['samtools', 'view', '-Sb', five_to_ref_sam, '>', five_to_ref_bam], shell=True)
                run_command(['samtools', 'view', '-Sb', three_to_ref_sam, '>', three_to_ref_bam], shell=True)
                run_command(['samtools', 'sort', five_to_ref_bam, five_bam_sorted], shell=True)
                run_command(['samtools', 'sort', three_to_ref_bam, three_bam_sorted], shell=True)
                run_command(['samtools', 'index', five_bam_sorted + '.bam'], shell=True)
                run_command(['samtools', 'index', three_bam_sorted + '.bam'], shell=True)
                # Create BED files with coverage information
                run_command(['bedtools', 'genomecov', '-ibam', five_bam_sorted + '.bam', '-bg', '>', five_cov_bed], shell=True)
                run_command(['bedtools', 'genomecov', '-ibam', three_bam_sorted + '.bam', '-bg', '>', three_cov_bed], shell=True)
                run_command(['bedtools', 'merge', '-d', args.merging, '-i', five_cov_bed, '>', five_cov_merged], shell=True)
                run_command(['bedtools', 'merge', '-d', args.merging, '-i', three_cov_bed, '>', three_cov_merged], shell=True)
                # Filter coveraged BED files on coverage cutoff (so only take
                # high coverage regions for further analysis)
                filter_on_depth(five_cov_bed, five_final_cov, args.cutoff)
                filter_on_depth(three_cov_bed, three_final_cov, args.cutoff)
                run_command(['bedtools', 'merge', '-d', args.merging, '-i', five_final_cov, '>', five_merged_bed], shell=True)
                run_command(['bedtools', 'merge', '-d', args.merging, '-i', three_final_cov, '>', three_merged_bed], shell=True)
                # Find intersects and closest points of regions
                run_command(['bedtools', 'intersect', '-a', five_merged_bed, '-b', three_merged_bed, '-wo', '>', bed_intersect], shell=True)
                run_command(['closestBed', '-a', five_merged_bed, '-b', three_merged_bed, '-d', '>', bed_closest], shell=True)
                # Create all possible closest bed files for checking unpaired hits
                run_command(['closestBed', '-a', five_merged_bed, '-b', three_cov_merged, '-d', '>', bed_unpaired_five], shell=True)
                run_command(['closestBed', '-a', five_cov_merged, '-b', three_merged_bed, '-d', '>', bed_unpaired_three], shell=True)
                # Create table and annotate genbank with hits
                if args.igv:
                    igv_flag = '1'
                else:
                    igv_flag = '0'
                run_command([args.path + 'create_typing_out.py', '--intersect', bed_intersect, '--closest', bed_closest,
                    '--left_bed', five_merged_bed, '--right_bed', three_merged_bed,
                    '--left_unpaired', bed_unpaired_five, '--right_unpaired', bed_unpaired_three,
                    '--seq', query, '--ref', args.typingRef, '--temp', temp_folder,
                    '--cds', args.cds, '--trna', args.trna, '--rrna', args.rrna, '--min_range', args.min_range,
                    '--max_range', args.max_range, '--output', current_dir + sample + '_' + query_name, '--igv', igv_flag, '--chr_name', args.chr_name], shell=True)

            # remove temp folder if required
            if args.temp == False:
                run_command(['rm', '-rf', temp_folder], shell=True)
            if args.bam == False:
                run_command(['rm', five_bam_sorted + '.bam', three_bam_sorted + '.bam', five_bam_sorted + '.bam.bai', three_bam_sorted + '.bam.bai'], shell=True)
    total_time = time.time() - start_time
    time_mins = float(total_time) / 60
    logging.info('ISMapper finished in ' + str(time_mins) + ' mins.')

if __name__ == '__main__':
    main()
