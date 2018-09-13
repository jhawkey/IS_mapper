import os
from run_commands import run_command, BedtoolsError, CommandError, make_directories
from mapping_to_query import create_tmp_file, bwa_index, RunSamtools
from create_output import create_typing_output

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
    output_filenames['closest'] = os.path.join(out_dir, prefix + '_' + ref_name + '_closest.bed')

    # set up table.txt
    output_filenames['table'] = os.path.join(out_dir, prefix + '_' + ref_name + '_table.txt')

    return(output_filenames)

def filter_on_depth(cov_file, out_bed, cov_cutoff):
    '''
    Takes a bed coverage file and removes lines that
    do not meet the coverage cutoff.
    Saves output to a new file.
    '''

    output = open(out_bed, 'w')
    with open(cov_file) as depth_info:
        for line in depth_info:
            if int(line.strip().split('\t')[3]) >= cov_cutoff:
                output.write(line)
    output.close()

def map_to_ref_seq(ref_seq, sample_name, left_flanking, right_flanking, tmp, out, bwa_threads, bwa_all):

    filenames = set_ref_output_filenames(sample_name, ref_seq.id, tmp, out)

    # make temp file
    ref_seq_file = create_tmp_file(ref_seq, filenames['ref_tmp'], 'fasta')

    # index the ref seq
    bwa_index(ref_seq_file)

    # set up samtools
    samtools_runner = RunSamtools()

    # Map reads to reference, reporting all alignments
    if bwa_all:
        run_command(['bwa', 'mem', '-t', bwa_threads, '-a', ref_seq_file, left_flanking, '>', filenames['left_sam']],
                    shell=True)
        run_command(['bwa', 'mem', '-t', bwa_threads, '-a', ref_seq_file, right_flanking, '>', filenames['right_sam']],
                    shell=True)

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

    return(filenames)

def create_bed_files(filenames, cutoff, merging):

    left_sorted = filenames['left_sorted']
    right_sorted = filenames['right_sorted']
    # Create BED files with coverage information
    run_command(['bedtools', 'genomecov', '-ibam', left_sorted, '-bg', '>', filenames['left_cov']], shell=True)
    run_command(['bedtools', 'genomecov', '-ibam', right_sorted, '-bg', '>', filenames['right_cov']], shell=True)
    run_command(['bedtools', 'merge', '-d', merging, '-i', filenames['left_cov'], '>', filenames['left_merged']], shell=True)
    run_command(['bedtools', 'merge', '-d', merging, '-i', filenames['right_cov'], '>', filenames['right_merged']], shell=True)
    # Filter coveraged BED files on coverage cutoff (so only take
    # high coverage regions for further analysis)
    filter_on_depth(filenames['left_cov'], filenames['left_final_cov'], cutoff)
    filter_on_depth(filenames['right_cov'], filenames['right_final_cov'], cutoff)

    run_command(['bedtools', 'merge', '-d', merging, '-i', filenames['left_final_cov'], '>', filenames['left_merged_bed']], shell=True)
    run_command(['bedtools', 'merge', '-d', merging, '-i', filenames['right_final_cov'], '>', filenames['right_merged_bed']], shell=True)

    # Find intersects of regions
    run_command(['bedtools', 'intersect', '-a', filenames['left_merged_bed'], '-b', filenames['right_merged_bed'], '-wo', '>',
                 filenames['intersect']], shell=True)
    
    # Find regions that are close but not overlapping
    try:
        run_command(['closestBed', '-a', filenames['left_merged_bed'], '-b', filenames['right_merged_bed'], '-d', '>', filenames['closest']], shell=True)
    # One or more of these files are empty so we need to quit and report no hits
    except BedtoolsError:
        #TODO: TEST THIS WORKS WHEN THIS ERROR IS THROWN
        #(filenames, ref_gbk_obj, is_query_obj, min_range, max_range, tmp_output_folder)
        create_typing_output(filenames, None, None, None, None, None)
        return(filenames)
        #with open(filenames['table'], 'w') as f:
            #header = ["region", "orientation", "x", "y", "gap", "call", "%ID", "%Cov", "left_gene", "left_strand",
            #          "left_distance", "right_gene", "right_strand", "right_distance", "functional_prediction"]
            #f.write('\t'.join(header) + '\nNo hits found')
            #continue
    
    # Check all unpaired hits to see if there are any that should be paired up
    # If any of these fail because there are no hits, just make empty unapired files to pass to create_typing_out
    try:
        run_command(['closestBed', '-a', filenames['left_merged_bed'], '-b', filenames['right_merged'], '-d', '>', filenames['left_unpaired']],
                    shell=True)
    except BedtoolsError:
        if not os.path.isfile(filenames['left_unpaired']) or os.stat(filenames['left_unpaired'])[6] == 0:
            open(filenames['left_unpaired'], 'w').close()
    try:
        run_command(['closestBed', '-a', filenames['left_merged'], '-b', filenames['right_merged_bed'], '-d', '>', filenames['right_unpaired']],
                    shell=True)
    except BedtoolsError:
        if not os.path.isfile(filenames['right_unpaired']) or os.stat(filenames['right_unpaired'])[6] == 0:
            open(filenames['right_unpaired'], 'w').close()
            
    # return the filepaths for all the output file names
    return(filenames)
