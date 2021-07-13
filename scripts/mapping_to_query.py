import os
import re
import shlex
from Bio import SeqIO
from Bio.Seq import Seq
import logging
from subprocess import Popen, PIPE
from run_commands import run_command, CommandError, make_directories

class RunSamtools:
    def __init__(self):
        # check iF SAMTOOLS environment variable exists, then use that command
        try:
            self.samtools_cmd = os.environ['SAMTOOLS']
        except:
            self.samtools_cmd = 'samtools'
        # check version:

        p = Popen([self.samtools_cmd], stderr = PIPE)
        out,err = p.communicate()
        version_string = [l for l in err.decode('UTF-8').split('\n') if re.search('^Version', l)][0]
        if len(version_string) == 0:
            logging.error("Could not find Samtools")
            raise IOError
        if len(re.findall('1\.[0-9]\.?[0-9]{0,2}', version_string)):
            version_id=re.findall('1\.[0-9]\.?[0-9]{0,2}', version_string)[0]
            logging.info("Found samtools version {}".format(version_id))
            self.version=1
        else:
            version_id=re.findall('0\.[0-9]\.?[0-9]{1,2}', version_string)[0]
            logging.info("Found samtools version {}".format(version_id))
            self.version=0
    def view(self, output_bam, input_sam, bigF=None, smallF=None):
        cmd = self.samtools_cmd + ' view -Sb'
        if bigF !=None:
            cmd = cmd + ' -F {}'.format(bigF)
        if smallF != None:
            cmd = cmd + ' -f {}'.format(smallF)
        cmd = cmd + ' -o {} {}'.format(output_bam, input_sam)
        return(shlex.split(cmd))
    def sort(self, output_bam, input_bam):
        cmd = self.samtools_cmd + ' sort'
        if self.version == 1:
            cmd = cmd + ' -T tmp -o {} {}'.format(output_bam, input_bam)
        else:
            cmd = cmd + ' {} {}'.format(input_bam, output_bam)
        return(shlex.split(cmd))
    def index(self, input_bam):
        cmd = self.samtools_cmd + ' index {}'.format(input_bam)
        return(shlex.split(cmd))
    def view_bam_to_sam(self, input_bam, output_sam):
        cmd = self.samtools_cmd + ' view'
        if self.version == 1:
            cmd = cmd + ' -O SAM {} > {}'.format(input_bam, output_sam)
            return(shlex.split(cmd))
        else:
            cmd = cmd + ' {} > {}'.format(input_bam, output_sam)
            return(shlex.split(cmd))

def bwa_index(fasta):

    """
    Check to see if bwa index for given input fasta exists.
    If it doesn't, build an index from the given input fasta.
    """

    built_index = fasta + '.bwt'
    if os.path.exists(built_index):
        logging.info('Index for {} is already built...'.format(fasta))
    else:
        logging.info('Building bwa index for {}...'.format(fasta))
        run_command(['bwa', 'index', fasta], shell=True)

def create_tmp_file(seq_object, file_path, file_type):

    """
    Take a folder name, a file name, and a file type.
    Concatenate the folder name and file name together, and write out
    the sequence with the correct file type (either genbank or fasta).
    Return the name of the file for use downstream.
    """
    SeqIO.write(seq_object, file_path, file_type)
    return file_path

def set_output_filenames(tmp_folder, prefix, query, out_dir):
    """
    Sets output filenames used in the mapping command.
    tmp_folder is the location of the temp folder for outputs which will be deleted
    prefix is the sample prefix
    query is the name of the IS query
    out_dir is the location of the folder where final outputs are kept
    Returns a dictionary, where the key is a shorthand for the output file, and the value is the file path
    """

    output_filenames = {}

    # set up tmp fasta
    output_filenames['query_tmp'] = os.path.join(tmp_folder, query + '.fasta')

    # set up sam
    output_filenames['sam'] = os.path.join(tmp_folder, prefix + '_' + query + '.sam')

    extensions = {'bam':'.bam', 'reads':'.fastq', 'clipped':'_clipped.fastq'}

    for type, ext in extensions.items():
        left = os.path.join(tmp_folder, prefix + '_' + query +'_left' + ext)
        right = os.path.join(tmp_folder, prefix + '_' + query + '_right' + ext)
        output_filenames['left_' + type] = left
        output_filenames['right_' + type] = right

    # set up final reads output
    output_filenames['left_final'] = os.path.join(out_dir, prefix + '_' + query + '_left_final.fastq')
    output_filenames['right_final'] = os.path.join(out_dir, prefix + '_' + query + '_right_final.fastq')

    return(output_filenames)

def extract_clipped_reads(sam_file, min_size, max_size, out_left_file, out_right_file):
    with open(sam_file, 'r') as in_file, open(out_left_file, 'w') as out_left, open(out_right_file, 'w') as out_right:
        logging.info("Extracting clip reads into" + out_left_file + out_right_file)
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
                    soft_clipped_seq = Seq(entries[9][:num_soft_clipped])
                    qual_scores = entries[10][:num_soft_clipped]
                    if reverse_complement:
                        out_left.write('@' + read_name + '\n' + str(soft_clipped_seq.reverse_complement()) + '\n+\n' + qual_scores[::-1] + '\n')
                    else:
                        out_left.write('@' + read_name + '\n' + str(soft_clipped_seq) + '\n+\n' + qual_scores + '\n')
            if map_regions[-1][-1] == 'S':
                num_soft_clipped = int(map_regions[-1][:-1])
                if min_size <= num_soft_clipped <= max_size:
                    soft_clipped_seq = Seq(entries[9][-num_soft_clipped:])
                    qual_scores = entries[10][-num_soft_clipped:]
                    if reverse_complement:
                        out_right.write('@' + read_name + '\n' + str(soft_clipped_seq.reverse_complement()) + '\n+\n' + qual_scores[::-1] + '\n')
                    else:
                        out_right.write('@' + read_name + '\n' + str(soft_clipped_seq) + '\n+\n' + qual_scores + '\n')

def map_to_is_query(sample, is_query, output_sample, min_clip, max_clip, threads):

    """
    Take the sample object (containing paths to reads and read prefix), the IS query (fasta file) and the
    output folder.
    - Create output folders for this IS within the sample folder
    - Set up output files (both temporary and final)
    - Create a temp file for the query
    - Index IS query and map reads to it
    - Extract unmapped reads flanking the IS query
    - Create fastq files from these resulting bam files
    - Extract reads which are clipped (partially mapped to the IS query)
    - Add these clipped reads to the fastq files

    Return the file names of the clipped reads for subsequent analysis.
    """
    samtools_runner = RunSamtools()

    # set up output folders
    is_query_out = os.path.join(output_sample, is_query.id)
    is_query_tmp_folder = os.path.join(output_sample, is_query.id, 'tmp')
    make_directories([is_query_out, is_query_tmp_folder])
    logging.info('Created output folder %s', is_query_out)

    # set up output file names
    filenames = set_output_filenames(is_query_tmp_folder, sample.prefix, is_query.id, is_query_out)

    # create temp file of IS query
    is_query_tmp = create_tmp_file(is_query, filenames['query_tmp'], 'fasta')

    # index the query
    bwa_index(is_query_tmp)
    # map to the query
    run_command(['bwa', 'mem', '-t', threads, is_query_tmp, str(sample.forward), str(sample.reverse), '>', filenames['sam']], shell=True)

    # pull out unmapped reads flanking IS
    run_command(samtools_runner.view(filenames['left_bam'], filenames['sam'], smallF=36), shell=True)
    run_command(samtools_runner.view(filenames['right_bam'], filenames['sam'], smallF=4, bigF=40), shell=True)

    # Turn bams to reads for mapping
    run_command(['bedtools', 'bamtofastq', '-i', filenames['left_bam'], '-fq', filenames['left_reads']], shell=True)
    run_command(['bedtools', 'bamtofastq', '-i', filenames['right_bam'], '-fq', filenames['right_reads']], shell=True)

    # Extract clipped reads
    logging.info('Extracting soft clipped reads that are <= %s bp and >= %s bp', str(max_clip), str(min_clip))
    extract_clipped_reads(filenames['sam'], min_clip, max_clip, filenames['left_clipped'], filenames['right_clipped'])

    # Add clipped reads to the final fastq files
    run_command(['cat', filenames['left_clipped'], filenames['left_reads'], '>', filenames['left_final']], shell=True)
    run_command(['cat', filenames['right_clipped'], filenames['right_reads'], '>', filenames['right_final']], shell=True)

    logging.info('Successfully extracted reads flanking left and right end of IS query')

    # return the paths to these reads
    return filenames['left_final'], filenames['right_final'], is_query_out, is_query_tmp_folder