#!/usr/bin/env python

import string, re
import os, sys
from argparse import (ArgumentParser, FileType)
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from operator import itemgetter
import os, sys, re, collections, operator
from collections import OrderedDict
from ismap import gbk_to_fasta
import time

def parse_args():

    parser = ArgumentParser(description="Create a table of IS hits in all isolates for ISMapper")
    # Inputs
    parser.add_argument('--tables', nargs='+', type=str, required=True, help='tables to compile')
    parser.add_argument('--reference_gbk', type=str, required=True, help='gbk file of reference to report closest genes')
    parser.add_argument('--seq', type=str, required=True, help='fasta file for insertion sequence looking for in reference')
    # Parameters for hits
    parser.add_argument('--gap', type=int, required=False, default=0, help='distance between regions to call overlapping')
    parser.add_argument('--cds', nargs='+', type=str, required=False, default=['locus_tag', 'gene', 'product'], help='qualifiers to look for in reference genbank for CDS features')
    parser.add_argument('--trna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for tRNA features')
    parser.add_argument('--rrna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for rRNA features')
    # Output parameters
    parser.add_argument('--output', type=str, required=True, help='name of output file')

    return parser.parse_args()

def check_ranges(ranges, range_to_check, gap, orientation):
    '''
    Takes a list of tuples with currently known ranges, and a new range
    to check against these to see if it overlaps.
    Also takes gap variable, indicating that the ranges may be gap distance
    apart and still merged.
    Also takes orientation, as only want to merge ranges that are the same 
    orientation.

    Returns the old range (to be replaced) and the new, merged range (to replace
        the old range with) and the orientation.
    If the range_to_check can't be merged with any of the known ranges, then
    just return False, False, False.
    '''

    # From ranges, create a list of tuples (min, max, orientation)
    list_of_range_tuples = []
    for key in ranges:
        list_of_range_tuples.append((min(key[0], key[1]), max(key[0], key[1]), ranges[key]))

    # get the largest value
    largest_value = max(list_of_range_tuples, key=operator.itemgetter(1))[1] + gap + 10

    # calculate the slice size
    slice_size = largest_value / len(list_of_range_tuples)

    #create our list of boxes
    range_boxes = []
    print 'largest, slice, box length'
    print largest_value
    print slice_size
    for i in range(0, len(list_of_range_tuples) + 10):
        range_boxes.append([])
    print len(range_boxes)
    #populate boxes
    for tup in list_of_range_tuples:
        index_1 = tup[0] / slice_size
        index_2 = tup[1] / slice_size
        while index_1 <= index_2:
            try:
                range_boxes[index_1].append(tup)
            except IndexError:
                print 'index1, index2, tuple'
                print index_1
                print index_2
                print tup
                exit(-1)
            index_1 += 1

    # find box for new range to check
    start = min(range_to_check[0], range_to_check[1])
    stop = max(range_to_check[1], range_to_check[0])

    #print 'slice, start, stop'
    #print slice_size
    #print start
    #print stop
    index_start = start / slice_size
    index_stop = stop / slice_size
    #print 'index start, index stop'
    #print index_start
    #print index_stop

    #print range_boxes
    #print len(range_boxes)

    # check each potential box

    while index_start <= index_stop and stop <= ((largest_value/slice_size) + 1):
        if range_boxes[index_start] != []:
            for tup in range_boxes[index_start]:
                if orientation == tup[2]:
                    x = tup[0]
                    y = tup[1]
                    # The x value must lie between start and stop in test range 
                    # taking into account gap
                    if x in range(start - gap, stop + 1) or x in range(start, stop + gap + 1):
                        # If so, then these ranges overlap
                        new_start = min(x, start)
                        new_end = max(y, stop)
                        return (x, y), (new_start, new_end), orientation
                    # Otherwise the y value must lie between the start and stop in the test range
                    # taking into account the gap
                    elif y in range(start - gap, stop + 1) or y in range(start, stop + gap + 1):
                        # If so, then these ranges overlap
                        new_start = min(x, start)
                        new_end = min(y, stop)
                        return (x, y), (new_start, new_end), orientation

        index_start += 1

    return False, False, False

def get_ref_positions(reference, is_query, positions_dict, orientation_dict):
    '''
    Get the coordinates of known IS sites in the reference.

    Takes the reference genbank, the IS query and the dictionary to add_argument
    IS query positions into, as well as a dictionary to add orientations
    of each of these positions.
    Returns these positions and orientations, as well as the reference name
    for file naming.
    '''
    # Get the name of the IS query to create temp file
    is_name = os.path.split(is_query)[1]
    ref_name = os.path.split(reference)[1]
    blast_output = os.getcwd() + '/' + is_name + '_' + ref_name + '.tmp'

    # Create a BLAST database of the reference if there isn't one already
    # Should probably use this function - blast_db(reference)
    if not os.path.exists(reference):
        os.system('makeblastdb -in ' + reference + ' -dbtype nucl')
    # Do the BLAST
    blastn_cline = NcbiblastnCommandline(query=is_query, db=reference, outfmt="'6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs'", out=blast_output)
    stdout, stderr = blastn_cline()
    # Open the BLAST output and get IS query sites
    with open(blast_output) as out:
        for line in out:
            info = line.strip('\n').split('\t')
            # To be a known site, hast to match query at least 90% with coverage of 95
            if float(info[3]) >= 90 and float(info[4])/float(info[1]) * 100 >= 95:
                x = int(info[6])
                y = int(info[7])
                positions_dict[(min(x, y), max(x, y))][ref_name] = '+'
                # Get orientation of known site for merging purposes
                if x > y:
                    orientation_dict[(min(x, y), max(x, y))] = 'R'
                else:
                    orientation_dict[(min(x, y), max(x, y))] = 'F'
    return positions_dict, orientation_dict, ref_name

def get_qualifiers(cds_qualifiers, trna_qualifiers, rrna_qualifiers, feature):
    '''
    Takes a list of possible qualifier IDs and attempts
    to find them in the feature given.
    If the qualifier is present, appends to a list, otherwise
    skips and keeps going.
    Returns a list of qualfiers found in that feature.
    '''
    
    return_quals = []
    if feature.type == 'CDS':
        qualifier_list = cds_qualifiers
    elif feature.type == 'tRNA':
        qualifier_list = trna_qualifiers
    elif feature.type == 'rRNA':
        qualifier_list = rrna_qualifiers
    for qual in qualifier_list:
        try:
            return_quals.append(feature.qualifiers[qual][0])
        except KeyError:
            pass
    return return_quals

def get_main_gene_id(qualifier_list, feature):
    '''
    Takes a list of qualifiers and a genbank feature.
    Returns the name of the qualifier that contains
    the gene id.
    '''

    for qual in qualifier_list:
        try:
            id_name = feature.qualifiers[qual][0]
            return id_name
        except KeyError:
            pass

def binary_search(features, isPosition, direction):
    min = 0
    max = len(features) - 1

    while True:

        # If the min has exceeded the max, then the IS position is not
        # inside a feature, and m will now be pointing to a
        # feature next to the IS position.
        if max < min:
            if direction == 'R':
                return findFeatureAfterPosition(features, isPosition, m)
            else:
                return findFeatureBeforePosition(features, isPosition, m)
        
        # Find the midpoint and save the feature attributes
        m = (min + max) // 2
        featureStart = features[m][0]
        featureEnd = features[m][1]
        featureIndex = features[m][2]       

        # If the IS position is after the feature, move the minimum to
        # be after the feature.
        if featureEnd < isPosition:
            min = m + 1

        # If the IS position is before the feature, move the maximum to
        # be before the feature.
        elif featureStart > isPosition:
            max = m - 1
        
        # If the IS position is inside the feature, return only that feature
        elif isPosition >= featureStart and isPosition <= featureEnd:
            return featureIndex

        else:
            return "1 - THIS SHOULDN'T HAPPEN!"

def findFeatureBeforePosition(features, isPosition, m):
    # If we are looking for the feature to the left of the
    # IS position, then either m-1 or m is our answer

    # If the start of the m feature is after the IS position,
    # then m is after the IS and m-1 is the correct feature
    if features[m][0] > isPosition:
        return features[m-1][2]

    # If both m and m+1 features are before the IS position,
    # then m will be closer to the IS and is the correct feature
    elif features[m-1][1] < isPosition and features[m][1] < isPosition:
        return features[m][2]

    else:
        return "2 - THIS SHOULDN'T HAPPEN!"

def findFeatureAfterPosition(features, isPosition, m):
    # If we are looking for the feature to the right of the
    # IS position, then either m or m+1 is our answer

    # If the end of the m feature is before the IS position,
    # then m is before the IS and m+1 is the correct feature
    if features[m][1] < isPosition:
        index = m + 1
        if index >= len(features):
            return features[0][2]
        return features[m+1][2]

    # If both m and m+1 features are after the IS position,
    # then m will be closer to the IS and is the correct feature
    elif features[m][0] > isPosition and features[m+1][0] > isPosition:
        return features[m][2]

    else:
        return "3 - THIS SHOULDN'T HAPPEN!"

def get_flanking_genes(features, feature_list, left, right, cds_features, trna_features, rrna_features):
    
    # Find the correct indexes
    left_feature_index = binary_search(feature_list, left, 'L')
    right_feature_index = binary_search(feature_list, right, 'R')
    # Extract the SeqFeature object that corresponds to that index
    left_feature = features[left_feature_index]
    right_feature = features[right_feature_index]

    # The info we require is:
    # [geneid, distance, [locus_tag, (gene), product, strand]]
    left_values = get_qualifiers(cds_features, trna_features, rrna_features, left_feature)
    right_values = get_qualifiers(cds_features, trna_features, rrna_features, right_feature)
    # Add the strand information
    left_values.append(str(left_feature.strand))
    right_values.append(str(right_feature.strand))
    # The distance to the left gene is the endmost position of the feature - the left IS coord
    left_dist = abs(max(left_feature.location.start, left_feature.location.end) - left)
    # The distance to the right gene is the startmost position of the feature - the right IS coord
    right_dist = abs(min(right_feature.location.start, right_feature.location.end) - right)
    # The first string in this values list is the main gene id (eg locus_tag)
    left_gene = [left_values[0], str(left_dist), left_values[1:]]
    right_gene = [right_values[0], str(right_dist), right_values[1:]]

    return left_gene, right_gene

def blast_db(fasta):
    '''
    Takes a fasta file and creates a BLAST database 
    if one doesn't exist already.
    '''
    
    if not os.path.exists(fasta + '.nin'):
        os.system('makeblastdb -in ' + fasta + ' -dbtype nucl')

def main():
    
    start_time = time.time()

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))
    list_of_isolates = []

    # key1 = (start, end), key2 = isolate, value = +/*/?
    list_of_positions = collections.defaultdict(dict)
    # key1 = (start, end), key2 = ref, value = +
    list_of_ref_positions = collections.defaultdict(dict)
    # key = (start, end), value = orientation (F/R)
    position_orientation = {}

    reference_fasta = args.reference_gbk.split('.g')[0]
    # Create a fasta file of the reference for BLAST
    print 'Creating fasta file and database of reference ...'
    gbk_to_fasta(args.reference_gbk, reference_fasta)
    # Make a BLAST database
    blast_db(reference_fasta)
    # Get the reference positions and orientations for this IS query
    print '\nGetting query positions in reference ...'
    list_of_positions, position_orientation, ref_name = get_ref_positions(reference_fasta, args.seq, list_of_positions, position_orientation)

    elapsed_time = time.time() - start_time
    print 'Time taken: ' + str(elapsed_time)
    # Loop through each table give to --tables
    print 'Collating results files ...'
    for result_file in unique_results_files:
        # Get isolate name
        isolate = result_file.split('_table.txt')[0]
        list_of_isolates.append(isolate)
        # Skip the header
        header = 0
        with open(result_file) as file_open:
            for line in file_open:
                # Skip header
                if header == 0:
                    header += 1
                # Check to make sure there were actually hits
                elif 'No hits found' not in line and line != '':
                    info = line.strip('\n').split('\t')
                    # Get orientation for hit and start/end coordinates
                    orientation = info[1]
                    is_start = min(int(info[2]), int(info[3]))
                    is_end = max(int(info[3]), int(info[2]))
                    # Note whether call is Known, Novel or Possible related IS
                    call = info[5]
                    # If the position Know and therefore in the reference, 
                    # compare the hit to the known reference positions
                    if call == 'Known' or call == 'Known?':
                        if (is_start, is_end) not in list_of_positions:
                            # Looking for hits that are with 100 bp of the Known reference hits
                            old_range, new_range, new_orientation = check_ranges(position_orientation, (is_start, is_end), 100, orientation)
                            # If we can merge ranges
                            if old_range != False:
                                store_values = list_of_positions[old_range]
                                # Remove the old range, and add the new range
                                del list_of_positions[old_range]
                                list_of_positions[new_range] = store_values
                                # Note whether the hit is uncertain (?)
                                # or confident (+)
                                if '?' in call:
                                    list_of_positions[new_range][isolate] = '?'
                                else:
                                    list_of_positions[new_range][isolate] = '+'
                                # Remove the old range from the reference positions
                                # and add the new merged range
                                del position_orientation[old_range]
                                position_orientation[new_range] = new_orientation
                            # If we can't merge with a known position
                            else:
                                # Mark as uncertain if ? in call
                                if '?' in call:
                                    list_of_positions[(is_start, is_end)][isolate] = '?'
                                # Otherwise just append it as a new reference position
                                else:
                                    list_of_positions[(is_start, is_end)][isolate] = '+'
                                # Note the orientation of the hit
                                position_orientation[(is_start, is_end)] = orientation
                    # Otherwise try and merge with positions that are novel
                    elif (is_start, is_end) not in list_of_positions:
                        # If the list of positions isn't empty, then there are ranges to check against
                        if list_of_positions.keys() != []:
                            old_range, new_range, new_orientation = check_ranges(position_orientation, (is_start, is_end), args.gap, orientation)
                            # So the current range overlaps with a range we already have
                            if old_range != False:
                                # Remove the old range and add the new one
                                store_values = list_of_positions[old_range]
                                del list_of_positions[old_range]
                                list_of_positions[new_range] = store_values
                                # Mark as ? if uncertain, * if imprecise
                                # or + if confident
                                if '?' in call:
                                    list_of_positions[new_range][isolate] = '?'
                                elif '*' in call:
                                    list_of_positions[new_range][isolate] = '*'
                                else:
                                    list_of_positions[new_range][isolate] = '+'
                                # Remove the old range from the orientations 
                                # and add the new one
                                del position_orientation[old_range]
                                position_orientation[new_range] = new_orientation
                            # Otherwise this range hasn't been seen before, so all values are False
                            else:
                                if '?' in call:
                                    list_of_positions[(is_start, is_end)][isolate] = '?'
                                elif '*' in call:
                                    list_of_positions[(is_start, is_end)][isolate] = '*'
                                else:
                                    list_of_positions[(is_start, is_end)][isolate] = '+'
                                position_orientation[(is_start, is_end)] = orientation
                        # Otherwise the position list is empty, so there are no ranges to check against
                        else:
                            if '?' in call:
                                    list_of_positions[(is_start, is_end)][isolate] = '?'
                            elif '*' in call:
                                list_of_positions[(is_start, is_end)][isolate] = '*'
                            else:
                                list_of_positions[(is_start, is_end)][isolate] = '+'
                            position_orientation[(is_start, is_end)] = orientation
                    # This position is already in the list, so just append the new isolate
                    elif (is_start, is_end) in list_of_positions:
                        if '?' in call:
                            list_of_positions[(is_start, is_end)][isolate] = '?'
                        elif '*' in call:
                            list_of_positions[(is_start, is_end)][isolate] = '*'
                        else:
                            list_of_positions[(is_start, is_end)][isolate] = '+'

    elapsed_time = time.time() - start_time
    print 'Time taken: ' + str(elapsed_time)

    # key = (start, end), valye = [left_gene, right_gene]
    position_genes = {}
    # Get the flanking genes for each know position
    print 'Getting flanking genes for each position (this step is the longest and could take some time) ...'

    # Get feature list
    gb = SeqIO.read(args.reference_gbk, "genbank")
    feature_list = []
    feature_count = 0
    feature_types = ["CDS", "tRNA", "rRNA"]

    for feature in gb.features:
        if feature.type in feature_types:
            feature_list.append([int(feature.location.start), int(feature.location.end), feature_count])
            feature_count += 1
        else:
            feature_count += 1

    '''if len(list_of_ref_positions.keys()) != 0:
        for position in list_of_ref_positions.keys():
            if position[0] < position[1]:
                left_pos = position[0]
                right_pos = position[1]
            else:
                left_pos = position[1]
                right_pos = position[0]
            flanking_left, flanking_right = get_flanking_genes(gb.features, feature_list, left_pos, right_pos, args.cds, args.trna, args.rrna)
            position_genes[(position[0], position[1])] = [flanking_left, flanking_right]'''
    # Get flanking genes
    for position in list_of_positions.keys():
        left_pos = min(position[0], position[1])
        right_pos = max(position[0], position[1])
        genes_before, genes_after =  get_flanking_genes(gb.features, feature_list, left_pos, right_pos, args.cds, args.trna, args.rrna)
        position_genes[(position[0], position[1])] = [genes_before, genes_after]

    # Order positions from smallest to largest for final table output
    order_position_list = list(OrderedDict.fromkeys(list_of_positions.keys()))
    order_position_list.sort()

    elapsed_time = time.time() - start_time
    print 'Time taken: ' + str(elapsed_time)

    # Create header of table
    print 'Writing output table to ' + args.output + ' ...'
    with open(args.output, 'w') as out:
        header = ['isolate']
        for position in order_position_list:
            if position_orientation[position] == 'F':
                header.append(str(position[0]) + '-' + str(position[1]))
            else:
                header.append(str(position[1]) + '-' + str(position[0]))
        out.write('\t'.join(header) + '\n')
        # Add the values for the reference positions
        row = [ref_name]
        for position in order_position_list:
            if position in list_of_positions:
                if ref_name in list_of_positions[position]:
                    row.append(list_of_positions[position][ref_name])
                else:
                    row.append('-')
        out.write('\t'.join(row) + '\n')
        
        # Loop through each isoalte
        # and create each row
        for isolate in list_of_isolates:
            row = [isolate]
            for position in order_position_list:
                if position in list_of_positions:
                    if isolate in list_of_positions[position]:
                        row.append(list_of_positions[position][isolate])
                    else:
                        row.append('-')
            out.write('\t'.join(row) + '\n')
        # Set up flanking genes
        row_orientation = ['orientation']
        row_l_locus = ['left ID']
        row_r_locus = ['right ID']
        row_l_dist = ['left distance']
        row_r_dist = ['right distance']
        row_l_strand = ['left strand']
        row_r_strand = ['right strand']
        row_l_prod = ['left info']
        row_r_prod = ['right info']

        # Print flanking genes for each position
        for position in order_position_list:
            #   genes_before, genes_after = get_flanking_genes(args.reference_gbk, position[0], position[1], args.cds, args.trna, args.rrna)
            row_orientation.append(position_orientation[position])
            if position in position_genes:
                #print position_genes[position]
                row_l_locus.append(position_genes[position][0][0])
                row_r_locus.append(position_genes[position][1][0])
                row_l_dist.append(position_genes[position][0][1])
                row_r_dist.append(position_genes[position][1][1])
                row_l_strand.append(position_genes[position][0][2][-1])
                row_r_strand.append(position_genes[position][1][2][-1])
                row_l_prod.append(position_genes[position][0][2])
                row_r_prod.append(position_genes[position][1][2])
        out.write('\t'.join(row_orientation) + '\n')
        out.write('\t'.join(row_l_locus) + '\n')
        out.write('\t'.join(row_l_dist) + '\n')
        out.write('\t'.join(row_l_strand) + '\n')
        out.write('\t'.join(str(i) for i in row_l_prod) + '\n')
        out.write('\t'.join(row_r_locus) + '\n')
        out.write('\t'.join(row_r_dist) + '\n')
        out.write('\t'.join(row_r_strand) + '\n')
        out.write('\t'.join(str(i) for i in row_r_prod) + '\n')

    elapsed_time = time.time() - start_time
    print 'Table compilation finished in ' + str(elapsed_time) 

if __name__ == "__main__":
    main()
