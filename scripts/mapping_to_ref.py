import os
import re
import shlex
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import logging
from subprocess import Popen, PIPE
from run_commands import run_command, BedtoolsError, CommandError, make_directories
from mapping_to_query import create_tmp_file, bwa_index, RunSamtools

def set_ref_output_filenames(prefix, ref_name, tmp_folder, out_dir):

    output_filenames = {}

    # temp file for ref sequence
    output_filenames['ref_tmp'] = os.path.join(tmp_folder, ref_name + '.fasta')

    extensions_tmp = {'sam': '.sam', 'bam': '.bam', 'cov': '_cov.bed', 'merged': '_cov_merged.sorted.bed'}
    extensions_final = {'sorted': '.sorted.bam', 'final_cov': '_finalcov.bed', 'merged_bed': '_merged.sorted.bed',
                        'unpaired': '_unpaired.bed'}

    for type, ext in extensions_tmp.items():
        left = os.path.join(tmp_folder, prefix + '_left_' + ref_name + ext)
        right = os.path.join(tmp_folder, prefix + '_right_' +  ref_name + ext)
        output_filenames['left_' + type] = left
        output_filenames['right_' + type] = right

    for type, ext in extensions_final.items():
        left = os.path.join(out_dir, prefix + '_left_' + ref_name + ext)
        right = os.path.join(out_dir, prefix + '_right_' +  ref_name + ext)
        output_filenames['left_' + type] = left
        output_filenames['right_' + type] = right


    output_filenames['intersect'] = os.path.join(out_dir, prefix + '_' + ref_name + '_intersect.bed')
    output_filenames['closest'] = os.path.join(out_dir, prefix + '_' + ref_name + '_intersect.bed')

    return(output_filenames)

    '''
    # Set up file names for output files
    left_header = sample + '_left_' + typingName
    right_header = sample + '_right_' + typingName
    left_to_ref_sam = temp_folder + left_header + '_' + query_name + '.sam'
    right_to_ref_sam = temp_folder + right_header + '_' + query_name + '.sam'
    left_to_ref_bam = temp_folder + left_header + '_' + query_name + '.bam'
    right_to_ref_bam = temp_folder + right_header + '_' + query_name + '.bam'
    left_bam_sorted = current_dir + left_header + '_' + query_name + '.sorted'
    right_bam_sorted = current_dir + right_header + '_' + query_name + '.sorted'
    left_cov_bed = temp_folder + left_header + '_' + query_name + '_cov.bed'
    right_cov_bed = temp_folder + right_header + '_' + query_name + '_cov.bed'
    left_cov_merged = temp_folder + left_header + '_' + query_name + '_cov_merged.sorted.bed'
    right_cov_merged = temp_folder + right_header + '_' + query_name + '_cov_merged.sorted.bed'
    left_final_cov = current_dir + left_header + '_' + query_name + '_finalcov.bed'
    right_final_cov = current_dir + right_header + '_' + query_name + '_finalcov.bed'
    left_merged_bed = current_dir + left_header + '_' + query_name + '_merged.sorted.bed'
    right_merged_bed = current_dir + right_header + '_' + query_name + '_merged.sorted.bed'
    bed_unpaired_left = current_dir + sample + '_' + typingName + '_' + query_name + '_left_unpaired.bed'
    bed_unpaired_right = current_dir + sample + '_' + typingName + '_' + query_name + '_right_unpaired.bed'
    bed_intersect = current_dir + sample + '_' + typingName + '_' + query_name + '_intersect.bed'
    bed_closest = current_dir + sample + '_' + typingName + '_' + query_name + '_closest.bed'
    '''


def map_to_ref_seq(ref_seq, sample_name, left_flanking, right_flanking, tmp, out, bwa_threads):

    filenames = set_ref_output_filenames(sample_name, ref_seq.id, tmp, out)
    print(filenames)

    # make temp file
    ref_seq_file = create_tmp_file(ref_seq, filenames['ref_tmp'], 'fasta')

    # index the ref seq
    bwa_index(ref_seq_file)

    # set up samtools
    samtools_runner = RunSamtools()

    # Map reads to reference, sort
    #TODO: add bwa -a option

    # map reads to the reference sequence
    run_command(['bwa', 'mem', '-t', bwa_threads, ref_seq_file, left_flanking, '>', filenames['left_sam']], shell=True)
    run_command(['bwa', 'mem', '-t', bwa_threads, ref_seq_file, right_flanking, '>', filenames['right_sam']], shell=True)

    # convert sams to bams
    run_command(samtools_runner.view(filenames['left_bam'], filenames['left_sam']), shell=True)
    run_command(samtools_runner.view(filenames['right_bam'], filenames['right_sam']), shell=True)
    # sort bams
    run_command(samtools_runner.sort(filenames['left_sorted'], filenames['left_bam']), shell=True)
    run_command(samtools_runner.sort(filenames['right_sorted'], filenames['right_bam']), shell=True)
    # index sorted bams
    run_command(samtools_runner.index(filenames['left_sorted']), shell=True)
    run_command(samtools_runner.index(filenames['right_sorted']), shell=True)

    return(filenames['left_sorted'], filenames['right_sorted'])

def create_bed_files(left_sorted, right_sorted, filenames, cutoff, merging):

    # Create BED files with coverage information
    run_command(['bedtools', 'genomecov', '-ibam', left_sorted + '.bam', '-bg', '>', filenames['left_cov']], shell=True)
    run_command(['bedtools', 'genomecov', '-ibam', right_sorted + '.bam', '-bg', '>', filenames['right_cov']], shell=True)
    run_command(['bedtools', 'merge', '-d', args.merging, '-i', filenames['left_cov'], '>', filenames['left_merged']], shell=True)
    run_command(['bedtools', 'merge', '-d', args.merging, '-i', filenames['right_cov'], '>', filenames['right_merged']], shell=True)
    # Filter coveraged BED files on coverage cutoff (so only take
    # high coverage regions for further analysis)
    filter_on_depth(filenames['left_cov'], filenames['left_final_cov'], cutoff)
    filter_on_depth(filenames['right_cov'], filenames['right_final_cov'], cutoff)

    run_command(['bedtools', 'merge', '-d', merging, '-i', filenames['left_final_cov'], '>', filenames['left_merged_bed']], shell=True)
    run_command(['bedtools', 'merge', '-d', merging, '-i', filenames['right_final_cov'], '>', filenames['right_merged_bed']], shell=True)

    # Find intersects and closest points of regions
    run_command(['bedtools', 'intersect', '-a', filenames['left_merged_bed'], '-b', filenames['right_merged_bed'], '-wo', '>',
                 filenames['intersect']], shell=True)

    return(filenames['intersect'])
