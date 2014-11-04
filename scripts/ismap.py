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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
try:
    from version import ismap_version
except:
    ismap_version = 'version unknown'

def parse_args():
    '''
    Parse the input arguments, use -h for help.
    '''

    parser = ArgumentParser(description='IS mapper')

    parser.add_argument("--version", action='version', version='v1.0')
    # Inputs
    parser.add_argument('--runtype', type=str, required=True, help='"typing" or "improvement"')
    parser.add_argument('--reads', nargs='+', type=str, required=False, help='Paired end reads for analysing (can be gzipped)')
    parser.add_argument('--forward', type=str, required=False, default='_1', help='Identifier for forward reads if not in MiSeq format (default _1)')
    parser.add_argument('--reverse', type=str, required=False, default='_2', help='Identifier for reverse reads if not in MiSeq format (default _2)')
    parser.add_argument('--query', type=str, required=True, help='Fasta file for query gene (eg: insertion sequence) that will be mapped to')
    parser.add_argument('--assemblies', nargs='+', type=str, required=False, help='Contig assemblies, one for each read set')
    parser.add_argument('--assemblyid', type=str, required=False, help='Identifier for assemblies eg: sampleName_contigs (specify _contigs) or sampleName_assembly (specify _assembly). Do not specify extension.')
    parser.add_argument('--extension', type=str, required=False, help='Extension for assemblies (eg: .fasta, .fa, .gbk, default is .fasta)', default='.fasta')
    parser.add_argument('--typingRef', type=str, required=False, help='Reference genome for typing against in genbank format')
    parser.add_argument('--type', type=str, required=False, default='fasta', help='Indicator for contig assembly type, genbank or fasta (default fasta)')
    parser.add_argument('--path', type=str, required=False, default='', help='Path to folder where scripts are (only required for development).')
    # Cutoffs for annotation
    parser.add_argument('--cutoff', type=int, required=False, default=6, help='Minimum depth for mapped region to be kept in bed file (default 6)')
    parser.add_argument('--percentid', type=float, required=False, default=90.0, help='Minimum percent ID for hit to be annotated (default 90.0')
    parser.add_argument('--merging', type=str, required=False, default='100', help='Value for merging left and right hits in bed files together to simply calculation of closest and intersecting regions (default 100).')
    parser.add_argument('--a', action='store_true', required=False, help='Switch on all alignment reporting for bwa')
    parser.add_argument('--T', type=str, required=False, default=30, help='Mapping quality score for bwa')
    # Options for table output (typing)
    parser.add_argument('--cds', nargs='+', type=str, required=False, default='locus_tag gene product', help='qualifiers to look for in reference genbank for CDS features')
    parser.add_argument('--trna', nargs='+', type=str, required=False, default='locus_tag product', help='qualifiers to look for in reference genbank for tRNA features')
    parser.add_argument('--rrna', nargs='+', type=str, required=False, default='locus_tag product', help='qualifiers to look for in reference genbank for rRNA features')
    # Reporting options
    parser.add_argument('--log', action='store_true', required=False, help='Switch on logging to file (otherwise log to stdout')
    parser.add_argument('--output', type=str, required=True, help='prefix for output files')
    parser.add_argument('--temp', action='store_true', required=False, help='Switch on keeping the temp folder instead of deleting it at the end of the program')

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
    If the directory is a VO directory, REMOVE IT before making it again,
    as this will cause Velvet to give an error.
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

    sequences = SeqIO.parse(genbank, "genbank")
    SeqIO.write(sequences, fasta, "fasta")

def multi_to_single(genbank, name, output):
    total = 0 # total bases

    handle = open(genbank, "rU")
    records = list(SeqIO.parse(handle, "genbank"))
    feature_count = 0
    colour_count = 0

    #make header genbank format friendly
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
            newrecord = SeqRecord(seq=r.seq,name=options.name,id=options.name)
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

def main():

    args = parse_args()

    if args.path != '' and args.path[-1] != '/':
        args.path = args.path + "/"
    # Set up logfile
    if args.log is True:
        logfile = args.output + ".log"
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

    # Checks to make sure the runtype is valid and provides an error
    # if it's not
    if args.runtype != "improvement" and args.runtype != "typing":
        logging.info('Invalid runtype selected: {}'.format(args.runtype))
        logging.info('Runtype should be improvement or typing (see instructions for further details)')
        exit(-1)

    # Gather together the reads in pairs with their corresponding
    # assemblies (if required)
    fileSets = read_file_sets(args)
    # Index the IS query for BWA
    bwa_index(args.query)
    # Start analysing each read set specified
    for sample in fileSets:
        forward_read = fileSets[sample][0]
        reverse_read = fileSets[sample][1]
        try:
            assembly = fileSets[sample][2]
        except IndexError:
            pass

        # Create the output file and folder names,
        # make the folders where necessary
        current_dir = os.getcwd() + '/'
        temp_folder = current_dir + sample + '_temp/'
        output_sam = temp_folder + sample + '.sam'
        five_bam = temp_folder + sample + '_5.bam'
        three_bam = temp_folder + sample + '_3.bam'
        five_reads = temp_folder + sample + '_5.fastq'
        three_reads = temp_folder + sample + '_3.fastq'
        no_hits_table = sample + '_table.txt'
        make_directories([temp_folder])

        # Map to IS query
        run_command(['bwa', 'mem', args.query, forward_read, reverse_read, '>', output_sam], shell=True)
        # Pull unmapped reads flanking IS
        run_command(['samtools view', '-Sb', '-f 36', output_sam, '>', five_bam], shell=True)
        run_command(['samtools view', '-Sb', '-f 4', '-F 40', output_sam, '>', three_bam], shell=True)
        # Turn bams to reads for mapping
        run_command(['bedtools', 'bamtofastq', '-i', five_bam, '-fq', five_reads], shell=True)
        run_command(['bedtools', 'bamtofastq', '-i', three_bam, '-fq', three_reads], shell=True)
        # Create BLAST database for IS query
        check_blast_database(args.query)
        if os.stat(five_reads)[6] == 0 or os.stat(three_reads)[6] == 0:
            logging.info('One or both read files are empty. This is probably due to no copies of the IS of interest being present in this sample. Program quitting.')
            with open(no_hits_table, 'w') as f:
                if args.runtype == 'typing':
                    header = ["region", "orientation", "x", "y", "gap", "call", "%ID", "%Cov", "left_gene", "left_strand", "left_distance", "right_gene", "right_strand", "right_distance", "functional_prediction"]
                    f.write('\t'.join(header) + '\nNo hits found')
                else:
                    header = ['contig', 'end', 'x', 'y']
                    f.write('\t'.join(header) + '\nNo hits found')
            sys.exit()

        if args.runtype == "improvement":

            # Get prefix for output filenames
            five_header = sample + '_5_assembly'
            three_header = sample + '_3_assembly'
            five_to_ref_sam = temp_folder + five_header + '.sam'
            three_to_ref_sam = temp_folder + three_header + '.sam'
            five_to_ref_bam = temp_folder + five_header + '.bam'
            three_to_ref_bam = temp_folder + three_header + '.bam'
            five_bam_sorted = five_header + '.sorted'
            three_bam_sorted = three_header + '.sorted'
            five_cov_bed = temp_folder + five_header + '_cov.bed'
            three_cov_bed = temp_folder + three_header + '_cov.bed'
            five_final_cov = five_header + '_finalcov.bed'
            three_final_cov = three_header + '_finalcov.bed'
            five_merged_bed = five_header + '_merged.sorted.bed'
            three_merged_bed = three_header + '_merged.sorted.bed'
            final_genbankSingle = sample + '_annotatedSingle.gbk'

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
                run_command(['bwa', 'mem', 'a', '-T', args.T, assembly, five_reads, '>', five_to_ref_sam], shell=True)
                run_command(['bwa', 'mem', 'a', '-T', args.T, assembly, three_reads, '>', three_to_ref_sam], shell=True)
            else:
                run_command(['bwa', 'mem', assembly, five_reads, '>', five_to_ref_sam], shell=True)
                run_command(['bwa', 'mem', assembly, three_reads, '>', three_to_ref_sam], shell=True)
            
            run_command(['samtools', 'view', '-Sb', five_to_ref_sam, '>', five_to_ref_bam], shell=True)
            run_command(['samtools', 'view', '-Sb', three_to_ref_sam, '>', three_to_ref_bam], shell=True)
            run_command(['samtools', 'sort', five_to_ref_bam, five_bam_sorted], shell=True)
            run_command(['samtools', 'sort', three_to_ref_bam, three_bam_sorted], shell=True)
            # Create BED file with coverage information
            run_command(['bedtools', 'genomecov', '-ibam', five_bam_sorted + '.bam', '-bg', '>', five_cov_bed], shell=True)
            run_command(['bedtools', 'genomecov', '-ibam', three_bam_sorted + '.bam', '-bg', '>', three_cov_bed], shell=True)
            filter_on_depth(five_cov_bed, five_final_cov, args.cutoff)
            filter_on_depth(three_cov_bed, three_final_cov, args.cutoff)
            run_command(['bedtools', 'merge', '-i', five_final_cov, '-d', args.merging, '>', five_merged_bed], shell=True)
            run_command(['bedtools', 'merge', '-i', three_final_cov, '-d', args.merging, '>', three_merged_bed], shell=True)       
            # Create table and genbank
            if args.extension == '.fasta':
                run_command([args.path + 'create_genbank_table.py', '--five_bed', five_merged_bed, '--three_bed', three_merged_bed, '--assembly', assembly, '--type fasta', '--output', sample], shell=True)
            elif args.extension == '.gbk':
                run_command([args.path + 'create_genbank_table.py', '--five_bed', five_merged_bed, '--three_bed', three_merged_bed, '--assembly', assembly_gbk, '--type genbank', '--output', sample], shell=True)
            #create single entry genbank
            multi_to_single(sample + '_annotated.gbk', sample, final_genbankSingle)

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
            five_to_ref_sam = temp_folder + five_header + '.sam'
            three_to_ref_sam = temp_folder + three_header + '.sam'
            five_to_ref_bam = temp_folder + five_header + '.bam'
            three_to_ref_bam = temp_folder + three_header + '.bam'
            five_bam_sorted = five_header + '.sorted'
            three_bam_sorted = three_header + '.sorted'
            five_cov_bed = temp_folder + five_header + '_cov.bed'
            three_cov_bed = temp_folder + three_header + '_cov.bed'
            five_final_cov = five_header + '_finalcov.bed'
            three_final_cov = three_header + '_finalcov.bed'
            five_merged_bed = five_header + '_merged.sorted.bed'
            three_merged_bed = three_header + '_merged.sorted.bed'
            bed_intersect = sample + '_' + typingName + '_intersect.bed'
            bed_closest = sample + '_' + typingName + '_closest.bed'

            # Map reads to reference, sort
            if args.a == True:
                run_command(['bwa', 'mem', '-a', '-T', args.T, typingRefFasta, five_reads, '>', five_to_ref_sam], shell=True)
                run_command(['bwa', 'mem', '-a', '-T', args.T, typingRefFasta, three_reads, '>', three_to_ref_sam], shell=True)
            else:
                run_command(['bwa', 'mem', typingRefFasta, five_reads, '>', five_to_ref_sam], shell=True)
                run_command(['bwa', 'mem', typingRefFasta, three_reads, '>', three_to_ref_sam], shell=True)

            run_command(['samtools', 'view', '-Sb', five_to_ref_sam, '>', five_to_ref_bam], shell=True)
            run_command(['samtools', 'view', '-Sb', three_to_ref_sam, '>', three_to_ref_bam], shell=True)
            run_command(['samtools', 'sort', five_to_ref_bam, five_bam_sorted], shell=True)
            run_command(['samtools', 'sort', three_to_ref_bam, three_bam_sorted], shell=True)
            run_command(['samtools', 'index', five_bam_sorted], shell=True)
            run_command(['samtools', 'index', three_bam_sorted], shell=True)
            # Create BED files with coverage information
            run_command(['bedtools', 'genomecov', '-ibam', five_bam_sorted + '.bam', '-bg', '>', five_cov_bed], shell=True)
            run_command(['bedtools', 'genomecov', '-ibam', three_bam_sorted + '.bam', '-bg', '>', three_cov_bed], shell=True)
            # Filter coveraged BED files on coverage cutoff (so only take 
            # high coverage regions for further analysis)
            filter_on_depth(five_cov_bed, five_final_cov, args.cutoff)
            filter_on_depth(three_cov_bed, three_final_cov, args.cutoff)
            run_command(['bedtools', 'merge', '-d', args.merging, '-i', five_final_cov, '>', five_merged_bed], shell=True)
            run_command(['bedtools', 'merge', '-d', args.merging, '-i', three_final_cov, '>', three_merged_bed], shell=True)
            # Find intersects and closest points of regions
            run_command(['bedtools', 'intersect', '-a', five_merged_bed, '-b', three_merged_bed, '-wo', '>', bed_intersect], shell=True)
            run_command(['closestBed', '-a', five_merged_bed, '-b', three_merged_bed, '-d', '>', bed_closest], shell=True)
            # Create table and annotate genbank with hits
            run_command([args.path + 'create_typing_out.py', '--intersect_bed', bed_intersect, '--closest_bed', bed_closest, '--insertion_seq', args.query, '--reference_genbank', args.typingRef, '--temp_folder', temp_folder, '--cds', args.cds, '--trna', args.trna, '--rrna', args.rrna, '--output', sample], shell=True)

        # remove temp folder if required
        if args.temp == False:
            run_command(['rm', '-rf', temp_folder], shell=True)

if __name__ == '__main__':
    main()
