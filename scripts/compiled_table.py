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
import time

class Position(object):
    def __init__(self, x, y, orientation, isolate_dict, left_feature, right_feature):
        self.x = x
        self.y = y
        self.orientation = orientation
        self.isolate_dict = isolate_dict
        self.left_feature = left_feature
        self.right_feature = right_feature
    def __eq__(self, other):
        return self.__dict__ == other.__dict__

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

def check_ranges(positions, range_to_check, gap, orientation):
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

    # From positions, create a list of tuples (min, max, orientation)
    list_of_range_tuples = []
    for pos in positions:
        list_of_range_tuples.append((pos.x, pos.y, pos.orientation))

    # get the largest value
    largest_value = max(list_of_range_tuples, key=operator.itemgetter(1))[1] + gap + 10

    # calculate the slice size
    slice_size = largest_value / len(list_of_range_tuples)

    #create our list of boxes
    range_boxes = []
    for i in range(0, len(list_of_range_tuples) + 10):
        range_boxes.append([])
    #populate boxes
    for tup in list_of_range_tuples:
        index_1 = tup[0] / slice_size
        index_2 = tup[1] / slice_size
        while index_1 <= index_2:
            try:
                range_boxes[index_1].append(tup)
            except IndexError:
                return False, False
            index_1 += 1

    # find box for new range to check
    start = min(range_to_check[0], range_to_check[1])
    stop = max(range_to_check[1], range_to_check[0])

    index_start = start / slice_size
    index_stop = stop / slice_size

    # check each potential box
    while index_start <= index_stop and index_stop <= ((largest_value/slice_size) + 1):
        if range_boxes[index_start] != []:
            for tup in range_boxes[index_start]:
                if orientation == tup[2]:
                    #print tup
                    x = tup[0]
                    y = tup[1]
                    # The x value must lie between start and stop in test range 
                    # taking into account gap
                    if (x in range(start - gap, stop + 1)) or (x in range(start, stop + gap + 1)):
                        # If so, then these ranges overlap
                        new_start = min(x, start)
                        new_end = max(y, stop)
                        for pos in positions:
                            if x == pos.x and y == pos.y:
                                matched_pos = pos
                        return matched_pos, (new_start, new_end)
                    # Otherwise the y value must lie between the start and stop in the test range
                    # taking into account the gap
                    elif (y in range(start - gap, stop + 1)) or (y in range(start, stop + gap + 1)):
                        # If so, then these ranges overlap
                        new_start = min(x, start)
                        new_end = max(y, stop)
                        for pos in positions:
                            if x == pos.x and y == pos.y:
                                matched_pos = pos
                        return matched_pos, (new_start, new_end)
                    # Also need to check if start and top lie within x and y
                    elif start in range(x - gap, y + 1) or start in range(x, y + gap + 1):
                        new_start = min(x, start)
                        new_end = max(y, stop)
                        for pos in positions:
                            if x == pos.x and y == pos.y:
                                matched_pos = pos
                        return matched_pos, (new_start, new_end)
                    elif stop in range(x - gap, y + 1) or stop in range(x, y + gap + 1):
                        new_start = min(x, start)
                        new_end = max(y, stop)
                        for pos in positions:
                            if x == pos.x and y == pos.y:
                                matched_pos = pos
                        return matched_pos, (new_start, new_end)
                    '''else:
                        print 'This is x and y, and then start and stop, then stop + 1'
                        print x
                        print y
                        print start
                        print stop
                        print stop + 1
                        print 'This is the gap'
                        print gap
                        print 'These are the x ranges to check'
                        print range(start - gap, stop + 1)
                        print range(start, stop + gap + 1)
                        print 'These are the y ranges to check'
                        print range(start - gap, stop + 1)
                        print range(start, stop + gap + 1)'''


        index_start += 1

    return False, False

def get_ref_positions(reference, is_query, positions_list):
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
                if x > y:
                    orientation = 'R'
                else:
                    orientation = 'F'
                pos_dict = {ref_name: '+'}
                new_pos = Position(min(x, y), max(x, y), orientation, pos_dict, None, None)
                #positions_dict[(min(x, y), max(x, y))][ref_name] = '+'
                # Get orientation of known site for merging purposes
                #if x > y:
                #    orientation_dict[(min(x, y), max(x, y))] = 'R'
                #else:
                #    orientation_dict[(min(x, y), max(x, y))] = 'F'
                positions_list.append(new_pos)
    return positions_list, ref_name

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
    print m
    print len(features)
    print features[m]
    print isPosition
    # an index error will occur if m is the final feature, so just check that the first part is true
    # and return m
    try:
        features[m+1]
    except IndexError:
        if features[m][0] > isPosition:
            return features[m][2]
        else:
            return "4 - THIS SHOULDN'T HAPPEN!"
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
    if type(left_feature_index) != int or type(right_feature_index) != int:
        print 'left index'
        print left_feature_index
        print 'right index'
        print right_feature_index
        print 'left position: ' + str(left)
        print 'right position: ' + str(right)
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

def gbk_to_fasta(genbank, fasta):
    '''
    Converts a genbank to a fasta using BioPython
    '''

    sequences = SeqIO.parse(genbank, "genbank")
    SeqIO.write(sequences, fasta, "fasta")

def final_ranges_check():

    # From ranges, create a list of tuples (min, max, orientation)
    list_of_range_tuples = []
    for key in ranges:
        list_of_range_tuples.append((min(key[0], key[1]), max(key[0], key[1]), ranges[key]))

    # sort the list into order based on x value
    list_of_range_tuples.sort()

    # start at the first index
    current_index = 0
    while current_index < len(list_of_range_tuples):
        # get x and y coordinates of range
        x = list_of_range_tuples[current_index][0]
        y = list_of_range_tuples[current_index][1]
        try:
            # check to see that the very next range is the same orientation,
            # otherwise we're not going to merge
            if list_of_range_tuples[current_index][2] == list_of_range_tuples[current_index + 1][2]:
                x2 = list_of_range_tuples[current_index + 1][0]
                y2 = list_of_range_tuples[current_index + 1][1]
                # if our previous y value is greater than or equal to the next-door
                # x value, then we should merge these columns
                if y >= x2:
                    new_position = (x, y2, list_of_range_tuples[current_index][2])
                    # need to check to see if this new position is able to be merged with any others

                # otherwise just store it as is
                else:
                    new_position = (x, y, list_of_range_tuples[current_index][2])
        except IndexError:
            # then we're at the last value so we just add that to the final list
            pass
        # increment our index
        current_index += 1



def main():
    
    start_time = time.time()

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))
    list_of_isolates = []

    # key1 = (start, end), key2 = isolate, value = +/*/?
    #list_of_positions = collections.defaultdict(dict)
    list_of_positions = []
    # key1 = (start, end), key2 = ref, value = +
    #list_of_ref_positions = collections.defaultdict(dict)
    # key = (start, end), value = orientation (F/R)
    #position_orientation = {}

    reference_fasta = args.reference_gbk.split('.g')[0]
    # Create a fasta file of the reference for BLAST
    print 'Creating fasta file and database of reference ...'
    gbk_to_fasta(args.reference_gbk, reference_fasta)
    # Make a BLAST database
    blast_db(reference_fasta)
    # Get the reference positions and orientations for this IS query
    print '\nGetting query positions in reference ...'
    list_of_positions, ref_name = get_ref_positions(reference_fasta, args.seq, list_of_positions)

    elapsed_time = time.time() - start_time
    print 'Time taken: ' + str(elapsed_time)
    #print list_of_positions
    #print ref_name
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
                    # See if this position is already in the list of positions
                    match = False
                    isolate_dict = {}
                    for pos in list_of_positions:
                        if pos.x == is_start and pos.y == is_end and pos.orientation == orientation:
                            # Then this position already exists
                            match = True
                            # And we want to retreive the position to which it is exactly the same
                            matching_pos = pos
                            # Then we want to add the info about this new position to the list
                            if '?' in call:
                                matching_pos.isolate_dict[isolate] = '?'
                            elif '*' in call:
                                matching_pos.isolate_dict[isolate] = '*'
                            else:
                                matching_pos.isolate_dict[isolate] = '+'
                    
                    # So we haven't seen this position before
                    if match == False:
                        # The position list is empty, so there's nothing to check against, so just add
                        # this new position
                        if list_of_positions == []:
                            if '?' in call:
                                isolate_dict[isolate] = '?'
                            elif '*' in call:
                                isolate_dict[isolate] = '*'
                            else:
                                isolate_dict[isolate] = '+'
                            new_pos = Position(is_start, is_end, orientation, isolate_dict, None, None)
                            list_of_positions.append(new_pos)
                        
                        # If the list of positions isn't empty, then there are ranges to check against
                        else:
                            old_position, new_range = check_ranges(list_of_positions, (is_start, is_end), args.gap, orientation)
                            # So the current range overlaps with a range we already have
                            if old_position != False:
                                isolate_dict = old_position.isolate_dict
                                # Add the new isolate to this dictionary
                                # Mark as ? if uncertain, * if imprecise
                                # or + if confident
                                if '?' in call:
                                    isolate_dict[isolate] = '?'
                                elif '*' in call:
                                    isolate_dict[isolate] = '*'
                                else:
                                    isolate_dict[isolate] = '+'
                                # Remove the old position from the list
                                list_of_positions.remove(old_position)
                                # Create the new position and add it
                                new_pos = Position(new_range[0], new_range[1], orientation, isolate_dict, None, None)
                                list_of_positions.append(new_pos)
                            # Otherwise this range hasn't been seen before, so all values are False
                            else:
                                if '?' in call:
                                    isolate_dict[isolate] = '?'
                                elif '*' in call:
                                    isolate_dict[isolate] = '*'
                                else:
                                    isolate_dict[isolate] = '+'
                                new_pos = Position(is_start, is_end, orientation, isolate_dict, None, None)
                                list_of_positions.append(new_pos)

    elapsed_time = time.time() - start_time
    print 'Time taken: ' + str(elapsed_time)
    
    # Get the flanking genes for each position now they've all been merged
    print 'Getting flanking genes for each position (this step is the longest and could take some time) ...'
    # key = (start, end), valye = [left_gene, right_gene]
    position_genes = {}

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
    # Sort the list just in case it's out of order (has caused issues in the past!!)
    feature_list = sorted(feature_list, key=itemgetter(0))
    # Get flanking genes
    for pos in list_of_positions:
        genes_before, genes_after =  get_flanking_genes(gb.features, feature_list, pos.x, pos.y, args.cds, args.trna, args.rrna)
        pos.left_feature = genes_before
        pos.right_feature = genes_after


    elapsed_time = time.time() - start_time
    print 'Time taken: ' + str(elapsed_time)

    # Order positions from smallest to largest for final table output
    list_of_positions.sort(key=lambda x: x.x)
    
    # Write out table
    print 'Writing output table to ' + args.output + ' ...'
    with open(args.output, 'w') as out:
        header = ['isolate']
        for pos in list_of_positions:
            if pos.orientation == 'F':
                header.append(str(pos.x) + '-' + str(pos.y))
            else:
                header.append(str(pos.y) + '-' + str(pos.x))
        out.write('\t'.join(header) + '\n')
        # Add the values for the reference positions
        row = [ref_name]
        for pos in list_of_positions:
            if ref_name in pos.isolate_dict.keys():
                row.append(pos.isolate_dict[ref_name])
            else:
                row.append('-')
        out.write('\t'.join(row) + '\n')
        
        # Loop through each isolate
        # and create each row
        for isolate in list_of_isolates:
            row = [isolate]
            for pos in list_of_positions:
                if isolate in pos.isolate_dict.keys():
                    row.append(pos.isolate_dict[isolate])
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

        # Print orientation and flanking genes for each position
        for pos in list_of_positions:
            row_orientation.append(pos.orientation)
            row_l_locus.append(pos.left_feature[0])
            row_r_locus.append(pos.right_feature[0])
            row_l_dist.append(pos.left_feature[1])
            row_r_dist.append(pos.right_feature[1])
            row_l_strand.append(pos.left_feature[2][-1])
            row_r_strand.append(pos.right_feature[2][-1])
            row_l_prod.append(pos.left_feature[2])
            row_r_prod.append(pos.right_feature[2])
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
