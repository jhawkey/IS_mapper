#!/usr/bin/env python

from argparse import (ArgumentParser, FileType)
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from operator import itemgetter
import os, sys, re, collections, operator
import numpy as np
from collections import OrderedDict
from compiled_table import get_flanking_genes, get_other_gene, get_qualifiers

def parse_args():

    parser = ArgumentParser(description="Create a table of features for ISMapper (typing pathway)")
    # Input files
    parser.add_argument('--intersect', type=str, required=True, help='intersection bed file')
    parser.add_argument('--closest', type=str, required=True, help='closest bed file')
    parser.add_argument('--left_bed', type=str, required=True, help='merged bed file for left end (5)')
    parser.add_argument('--right_bed', type=str, required=True, help='merged bed file for right end (3)')
    parser.add_argument('--left_unpaired', type=str, required=True, help='closest bed file where left end is full coverage')
    parser.add_argument('--right_unpaired', type=str, required=True, help='closest bed file where right end is full coverage')
    parser.add_argument('--ref', type=str, required=True, help='reference genbank file to find flanking genes of regions')
    parser.add_argument('--seq', type=str, required=True, help='insertion sequence reference in fasta format')
    # Flanking gene parameters
    parser.add_argument('--cds', nargs='+', type=str, required=False, default=['locus_tag', 'gene', 'product'], help='qualifiers to look for in reference genbank for CDS features (default locus_tag gene product)')
    parser.add_argument('--trna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for tRNA features (default locus_tag product)')
    parser.add_argument('--rrna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for rRNA features (default locus_tag product)')
    # Parameters for determining known genes
    parser.add_argument('--min_range', type=float, required=False, default=0.2, help='Minimum percent size of the gap to be called a known hit (default 0.5, or 50 percent)')
    parser.add_argument('--max_range', type=float, required=False, default=1.1, help='Maximum percent size of the gap to be called a known hit (default 1.1, or 110 percent)')
    # Output parameters
    parser.add_argument('--temp', type=str, required=True, help='location of temp folder to place intermediate blast files in')
    parser.add_argument('--output', type=str, required=True, help='name for output file')
    return parser.parse_args()

def insertion_length(insertion):
    '''
    Find the size of the IS query.
    '''

    sequence = SeqIO.read(insertion, "fasta")
    length = len(sequence.seq)

    return length

def doBlast(blast_input, blast_output, database):
    '''
    Perform a BLAST using the NCBI command line tools 
    in BioPython.
    '''
    blastn_cline = NcbiblastnCommandline(query=blast_input, db=database, outfmt="'6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs'", out=blast_output)
    stdout, stderr = blastn_cline()

def check_seq_between(gb, insertion, start, end, name, temp):
    '''
    Check the sequence between two ends to see
    if it matches the IS query or not, and what
    the coverage and %ID to the query.
    '''

    genbank = SeqIO.read(gb, 'genbank')
    # Get sequence between left and right ends
    seq_between = genbank.seq[start:end]
    # Turn the sequence into a fasta file
    seq_between = SeqRecord(Seq(str(seq_between), generic_dna), id=name)
    SeqIO.write(seq_between, temp + name + '.fasta', 'fasta')
    # Perform the BLAST
    doBlast(temp + name + '.fasta', temp + name + '_out.txt', insertion)
    # Only want the top hit, so set count variable to 0
    first_result = 0
    # Open the BLAST output file 
    with open(temp + name + '_out.txt') as summary:
        for line in summary:
            # Get coverage and % ID for top hit
            if first_result == 0:
                info = line.strip().split('\t')
                coverage = float(info[4]) / float(info[5]) * 100
                hit = [info[3], coverage]
                first_result += 1
            return hit
    # If there is not hit, just return an empty list
    hit = []
    return []

def createFeature(hits, orient):
    '''
    Create a feature for the hit to
    be added to the genbank.
    '''

    # Get coordinates
    x_L = hits[0]
    y_L = hits[1]
    x_R = hits[2]
    y_R = hits[3]
    quals = {}

    # Set locations for left and right hits
    left_location = SeqFeature.FeatureLocation(x_L, y_L)
    right_location = SeqFeature.FeatureLocation(x_R, y_R)
    if orient == 'F':
        #then in forward orientation, set colour to be red
        quals['colour'] = '2'
        quals['orientation'] = 'forward'
    elif orient == 'R':
        #then in reverse orientation, set colour to be yellow
        quals['colour'] = '7'
        quals['orientation'] = 'reverse'
    # Create features
    left_feature = SeqFeature.SeqFeature(left_location, type='left_end', qualifiers=quals)
    right_feature = SeqFeature.SeqFeature(right_location, type='right_end', qualifiers=quals)

    # Return features
    return left_feature, right_feature

def novel_hit(x_L, y_L, x_R, y_R, x, y, genbank, ref, cds, trna, rrna, gap, orient, feature_count, region, results, unpaired=False, star=False):
    '''
    Get flanking gene information for novel hits.
    '''
    
    # Create features for genbank
    left_feature, right_feature = createFeature([x_L, y_L, x_R, y_R], orient)
    # Add features to genbank
    genbank.features.append(left_feature)
    genbank.features.append(right_feature)
    # Increment the feature count
    feature_count += 2
    
    # Get the genes flanking the left and right ends
    gene_left, gene_right = get_flanking_genes(ref, x, y, cds, trna, rrna)
    #print gene_left
    #print gene_right
    # If the genes are the same, then hit is inside the gene
    if gene_left[-1] == gene_right[-1]:
        func_pred = 'Gene interrupted'
    # Otherwise need to do some more looking to determine functional predition
    else:
        func_pred = functional_prediction(gene_left, gene_right)
    
    # This is a confident hit
    if unpaired == False:
        call = 'Novel'
    # Hit is paired with a low coverage end, so an unconfident hit
    elif unpaired == True:
        call = 'Novel?'
    # This hit is imprecise, as gap size is larger than expected
    if star == True:
        call = 'Novel*'
    
    # Store all information for final table output
    results['region_' + str(region)] = [orient, str(x), str(y), gap, call, '', '', gene_left[-1][:-1], gene_left[-1][-1], gene_left[1], gene_right[-1][:-1], gene_right[-1][-1], gene_right[1], func_pred]

def functional_prediction(gene_left, gene_right):
    '''
    Determine how far upstream/downstream the IS insertion
    is from flanking genes. 
    '''

    # Get distance to left gene (start codon)
    bases = gene_left[1][1:]
    #print bases
    # If left gene is on + strand, we're upstream, otherwise we're downstream
    if '+' in gene_left[1]:
        prediction = 'Upstream of ' + gene_left[-1][0] + ' by ' + bases + 'bp, '
    elif '-' in gene_left[1]:
        prediction = 'Downstream of ' + gene_left[-1][0] + ' by ' + bases + 'bp, '

    # Get distance to right gene (start codon)
    bases = gene_right[1][1:]
    #print bases
    # If right gene is on + strand, we're upstream, otherwise we're downstream
    if '+' in gene_right[1]:
        prediction += 'upstream of ' + gene_right[-1][0] + ' by ' + bases + 'bp'
    elif '-' in gene_right[1]:
        prediction += 'downstream of ' + gene_right[-1][0] + ' by ' + bases + 'bp'

    if '+' not in gene_left[1] and '-' not in gene_left[1] and '+' not in gene_right[1] and '-' not in gene_right[1]:
        prediction = 'Gene interrupted'
    #print prediction

    return prediction

def add_known(x_L, x_R, y_L, y_R, gap, genbank, ref, seq, temp, cds, trna, rrna, region, feature_count, results, removed_results, line, file_loc):
    '''
    Adds a value to the table that is a known hit
    '''
    # Get orientation
    if y_L < x_R:
        start = y_L
        end = x_R
        orient = 'F'
    else:
        start = y_R
        end = x_L
        orient = 'R'
    # Get features and append to genbank
    left_feature, right_feature = createFeature([x_L, y_L, x_R, y_R], orient)
    genbank.features.append(left_feature)
    genbank.features.append(right_feature)
    # Increment number of features found
    feature_count += 2
    # Check to see if the sequence between actually belongs to the IS query
    seq_results = check_seq_between(ref, seq, start, end, 'region_' + str(region), temp)
    # This is a known site of coverage and %ID above 80
    if len(seq_results) != 0 and seq_results[0] >= 80 and seq_results[1] >= 80:
        # Taking all four coordinates and finding min and max to avoid coordinates 
        # that overlap the actual IS (don't want to return those in gene calls)
        # Mark as a known call to improve accuracy of gene calling
        #print 'setting known to true'
        gene_left = get_other_gene(ref, min(y_L, y_R, x_R, x_L), "left", cds, trna, rrna, known=True)
        gene_right = get_other_gene(ref, max(y_L, y_R, x_R, x_L), "right", cds, trna, rrna, known=True)
        #print gene_left
        #print gene_right
        # If the genes are the same, then this gene must be interrupted by the known site
        if gene_left[0] == gene_right[0]:
            func_pred == 'Gene interrupted'
            # Remove + and - from distance as the gene is interrupted
            gene_right[1] = gene_right[1][:-1]
            gene_left[1] = gene_left[1][:-1]
        # Otherwise we need to determine who is upstream/downstream of what
        else:
            func_pred = functional_prediction(gene_left, gene_right)
        # Add to the final results
        if 'unpaired' in file_loc:
            call = 'Known?'
        else:
            call = 'Known'
        results['region_' + str(region)] = [orient, str(start), str(end), gap, call, str(seq_results[0]), str('%.2f' % seq_results[1]), gene_left[-1][:-1], gene_left[-1][-1], gene_left[1], gene_right[-1][:-1], gene_right[-1][-1], gene_right[1], func_pred]
    else:   
        # Then I'm not sure what this is
        # Get flanking genes anyway
        gene_left, gene_right = get_flanking_genes(ref, start, end, cds, trna, rrna)
        func_pred = functional_prediction(gene_left, gene_right)
        if 'unpaired' in file_loc:
            call = 'Possible related IS?'
        else:
            call = 'Possible releated IS'
        if len(seq_results) !=0:
            results['region_' + str(region)] = [orient, str(start), str(end), gap, call, str(seq_results[0]), str('%.2f' % seq_results[1]), gene_left[-1][:-1], gene_left[-1][-1], gene_left[1], gene_right[-1][:-1], gene_right[-1][-1], gene_right[1], func_pred]
        else:
            removed_results['region_' + str(region)] = line.strip() + '\t' + file_loc +'\n'                

def main():

    args = parse_args()

    # Setup variables: results - for final table, removed_results - table showing
    # results which didn't pass cutoff tests, 
    # region - , lines - , header - header for final table
    results = {}
    removed_results = {}
    region = 1
    lines = 0
    header = ["region", "orientation", "x", "y", "gap", "call", "%ID", "%Cov", "left_gene", "left_strand", "left_distance", "right_gene", "right_strand", "right_distance", "functional_prediction"]
    
    # If both intersect and bed files are empty, there are no hits
    if os.stat(args.intersect)[6] == 0 and os.stat(args.closest)[6] == 0:
        output = open(args.output + '_table.txt', 'w')
        output.write('\t'.join(header) + '\n')
        output.write('No hits found')
        output.close()
        # Exit ISMapper
        sys.exit()

    # Read in genbank and intialise feature count
    genbank = SeqIO.read(args.ref, 'genbank')
    feature_count = 0

    intersect_left = []
    intersect_right = []
    closest_left = []
    closest_right = []
    # Start with the intersect file (novel hits)
    if os.stat(args.intersect)[6] != 0:
        with open(args.intersect) as bed_merged:
            for line in bed_merged:
                info = line.strip().split('\t')
                intersect_left.append(info[0:3])
                intersect_right.append(info[3:6])  
                # Set up coordinates for checking: L is the left end of the IS (5') 
                # and R is the right end of the IS (3')
                # Eg: x_L and y_L are the x and y coordinates of the bed block that 
                # matches to the region which is flanking the left end or 5' of the IS
                x_L = int(info[1])
                y_L = int(info[2])
                x_R = int(info[4])
                y_R = int(info[5])
                # Check to see if the gap is reasonable
                if int(info[6]) <= 15:
                    R_range = range(min(x_R, y_R), max(x_R, y_R))
                    L_range = range(min(x_L, y_L), max(x_L, y_L))
                    # if one hit is inside the other hit, remove it - don't know what to do with these
                    if (x_L in R_range and y_L in R_range) or (x_R in L_range and y_R in L_range):
                        removed_results['region_' + str(lines)] = line.strip() + '\tOne hit inside the other, intersect.bed\n'
                    else:
                        # Get orientation
                        if x_L < x_R or y_L < y_R:
                            orient = 'F'
                            x = x_R
                            y = y_L
                        elif x_L > x_R or y_L > y_R:
                            orient = 'R'
                            x = x_L
                            y = y_R
                        else:
                            #print 'neither if statement were correct'
                            pass
                        # Create result
                        novel_hit(x_L, y_L, x_R, y_R, x, y, genbank, args.ref, args.cds, args.trna, args.rrna, info[6], orient, feature_count, region, results, unpaired=False)
                        region += 1
                    # Otherwise we're removing this region, but keeping the information
                    # so the user can check later
                    else:
                        removed_results['region_' + str(lines)] = line.strip() + '\tintersect.bed\n'
                lines += 1
    
    # Get size of IS query
    is_length = insertion_length(args.seq)
    # Look inside the closest file (known or imprecise hits)
    with open(args.closest) as bed_closest:
        for line in bed_closest:
            info = line.strip().split('\t')
            closest_left.append(info[0:3])
            closest_right.append(info[3:6])
            # If the fourth column contains -1, there are no closest hits
            if info[3] == '-1':
                output = open(args.output, 'w')
                output.write('\t'.join(header) + '\n')
                output.write('No hits found')
                output.close()
                # Exit ISMapper
                sys.exit()
            # Get coordinate info
            x_L = int(info[1])
            y_L = int(info[2])
            x_R = int(info[4])
            y_R = int(info[5])
            R_range = range(min(x_R, y_R), max(x_R, y_R))
            L_range = range(min(x_L, y_L), max(x_L, y_L))
            # if one hit is inside the other hit, remove it - don't know what to do with these
            if (x_L in R_range and y_L in R_range) or (x_R in L_range and y_R in L_range):
                removed_results['region_' + str(region)] = line.strip() + '\tOne hit inside the other, closest.bed\n'
            else:
                # Set orientation
                if x_L < x_R and y_L < y_R:
                    orient = 'F'
                    x = y_L
                    y = x_R
                elif x_L > x_R and y_L > y_R:
                    orient = 'R'
                    x = x_L
                    y = y_R
                # If the gap column = 0, then they are intersected
                # This will be in the intersect file, so ignore
                if int(info[6]) == 0:
                    pass
                # This is probably a novel hit where there was no overlap detected
                elif int(info[6]) <= 10:
                    novel_hit(x_L, y_L, x_R, y_R, x, y, genbank, args.ref, args.cds, args.trna, args.rrna, info[6], orient, feature_count, region, results, unpaired=False)
                    region += 1
                # This is probably a known hit, but need to check with BLAST
                # Only a known hit if we're in the a range between (default 0.5 and 1.5) the size
                # of the IS query
                elif float(info[6]) / is_length >= args.min_range and float(info[6]) / is_length <= args.max_range:
                    add_known(x_L, x_R, y_L, y_R, info[6], genbank, args.ref, args.seq, args.temp, args.cds, args.trna, args.rrna, region, feature_count, results, removed_results, line, 'closest.bed')
                    region += 1
                # Could possibly be a novel hit but the gap size is too large
                elif float(info[6]) / is_length <= args.min_range and float(info[6]) / is_length < args.max_range:
                    novel_hit(x_L, y_L, x_R, y_R, x, y, genbank, args.ref, args.cds, args.trna, args.rrna, info[6], orient, feature_count, region, results, unpaired=False,star=True)
                    region +=1
                # This is something else altogether - either the gap 
                # is really large or something, place it in removed_results
                else:
                    removed_results['region_' + str(region)] = line.strip() + '\tclosest.bed\n'
                    region += 1

    # Looking for unpaired hits which are not in the merged/closest bed files
    # Possibly unpaired because the pair is low coverage and didn't pass
    # depth cutoff
    line_check = []
    with open(args.left_bed) as left_bed:
        for line in left_bed:
            if line.strip().split('\t') not in intersect_left and line.strip().split('\t') not in closest_left:
                line_check.append(line.strip().split('\t'))
    if len(line_check) != 0:
        with open(args.left_unpaired) as left_unpaired:
            for line in left_unpaired:
                info = line.strip().split('\t')
                # This is an unpaired hit
                if line.strip().split('\t')[0:3] in line_check:
                    # Get coordinates
                    x_L = int(info[1])
                    y_L = int(info[2])
                    x_R = int(info[4])
                    y_R = int(info[5])
                    R_range = range(min(x_R, y_R), max(x_R, y_R))
                    L_range = range(min(x_L, y_L), max(x_L, y_L))
                    # if one hit is inside the other hit, remove it - don't know what to do with these
                    if (x_L in R_range and y_L in R_range) or (x_R in L_range and y_R in L_range):
                        removed_results['region_' + str(lines)] = line.strip() + '\tOne hit inside the other, intersect.bed\n'
                    else:
                        # Get orientation
                        if x_L < x_R and y_L < y_R:
                            orient = 'F'
                            x = x_R
                            y = y_L
                        elif x_L > x_R and y_L > y_R:
                            orient = 'R'
                            x = x_L
                            y = y_R
                        # This ia novel hit
                        if float(info[6]) <= 10:
                            novel_hit(x_L, y_L, x_R, y_R, x, y, genbank, args.ref, args.cds, args.trna, args.rrna, info[6], orient, feature_count, region, results, unpaired=True)
                            region += 1
                        # This is a known hit
                        elif float(info[6]) / is_length >= args.min_range and float(info[6]) / is_length <= args.max_range:
                            add_known(x_L, x_R, y_L, y_R, info[6], genbank, args.ref, args.seq, args.temp, args.cds, args.trna, args.rrna, region, feature_count, results, removed_results, line, 'left_unpaired.bed')
                            region += 1
                        # Could possibly be a novel hit but the gap size is too large
                        elif float(info[6]) / is_length <= args.min_range and float(info[6]) / is_length < args.max_range:
                            # Add to results file
                            novel_hit(x_L, y_L, x_R, y_R, x, y, genbank, args.ref, args.cds, args.trna, args.rrna, info[6], orient, feature_count, region, results, unpaired=True)
                            region +=1
                        # This is something else altogether - either the gap is
                        # really large or something, place it in removed_results
                        else:
                            removed_results['region_' + str(region)] = line.strip() + '\tleft_unpaired.bed\n'
                            region += 1
    line_check = []
    with open(args.right_bed) as right_bed:
        for line in right_bed:
            if line.strip().split('\t') not in intersect_right and line.strip().split('\t') not in closest_right:
                line_check.append(line.strip().split('\t'))
    if len(line_check) != 0:
        with open(args.right_unpaired) as right_unpaired:
            for line in right_unpaired:
                info = line.strip().split('\t')
                #this is an unpaired hit
                if line.strip().split('\t')[3:6] in line_check:
                    #get coordinate info
                    x_L = int(info[1])
                    y_L = int(info[2])
                    x_R = int(info[4])
                    y_R = int(info[5])
                    R_range = range(min(x_R, y_R), max(x_R, y_R))
                    L_range = range(min(x_L, y_L), max(x_L, y_L))
                    # if one hit is inside the other hit, remove it - don't know what to do with these
                    if (x_L in R_range and y_L in R_range) or (x_R in L_range and y_R in L_range):
                        removed_results['region_' + str(lines)] = line.strip() + '\tOne hit inside the other, intersect.bed\n'
                    else:
                        #get orientation
                        if x_L < x_R and y_L < y_R:
                            orient = 'F'
                            x = x_R
                            y = y_L
                        elif x_L > x_R and y_L > y_R:
                            orient = 'R'
                            x = x_L
                            y = y_R
                        #a novel hit
                        if float(info[6]) <= 10:
                            novel_hit(x_L, y_L, x_R, y_R, x, y, genbank, args.ref, args.cds, args.trna, args.rrna, info[6], orient, feature_count, region, results, unpaired=True)
                            region += 1
                        #a known hit
                        elif float(info[6]) / is_length >= args.min_range and float(info[6]) / is_length <= args.max_range:
                            add_known(x_L, x_R, y_L, y_R, info[6], genbank, args.ref, args.seq, args.temp, args.cds, args.trna, args.rrna, region, feature_count, results, removed_results, line, 'right_unpaired.bed')               
                            region += 1
                        #could possibly be a novel hit but the gap size is too large
                        elif float(info[6]) / is_length <= args.min_range and float(info[6]) / is_length < args.max_range:
                            novel_hit(x_L, y_L, x_R, y_R, x, y, genbank, args.ref, args.cds, args.trna, args.rrna, info[6], orient, feature_count, region, results, unpaired=True)
                            region +=1
                        #this is something else altogether - either the gap is really large or something, place it in removed_results
                        else:
                            removed_results['region_' + str(region)] = line.strip() + '\tright_unpaired.bed\n'
                            region += 1

    # Sort regions into the correct order
    table_keys = []
    for key in results:
        table_keys.append(key)
    region_indexes = []
    for region in table_keys:
        region_indexes.append(region.split('region_')[1])
    arr = np.vstack((table_keys, region_indexes)).transpose()
    if arr != 0:
        sorted_keys = arr[arr[:,1].astype('int').argsort()]

    # Write out the found hits to file
    output = open(args.output + '_table.txt', 'w')
    output.write('\t'.join(header) + '\n')
    if arr != 0:
        for key in sorted_keys[:,0]:
            output.write(key + '\t' + '\t'.join(str(i) for i in results[key]) + '\n')
    if arr == 0:
        output.write('No hits found.')
    output.close()

    # Write out hits that were removed for whatever reason to file
    if len(removed_results) != 0:
        output_removed = open(args.output + '_removedHits.txt', 'w')
        for region in removed_results:
            output_removed.write(removed_results[region])
        output_removed.close()

    SeqIO.write(genbank, args.output + '_annotated.gbk', 'genbank')
    print('Added ' + str(feature_count) + ' features to ' + args.output + '_annotated.gbk')

    #return(lines, len(removed_results))

if __name__ == "__main__":
    
    main()