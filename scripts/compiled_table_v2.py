#!/usr/bin/env python

import os
import operator
from Bio import SeqIO
from create_output import ISHit, doBlast, get_features
from run_commands import run_command
from collections import OrderedDict
import time
import argparse

class Position(ISHit):

    def __init__(self, left_pos, right_pos):

        self.x = left_pos
        self.y = right_pos

        self.isolate_dict = {}

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
    blast_output = os.path.join(os.getcwd(), is_name + '_' + ref_name + '.tmp')

    # Do a BLAST of the IS query and the reference
    doBlast(is_query, blast_output, reference)

    # Open the BLAST output and get IS query sites
    with open(blast_output) as out:
        for line in out:
            fields = line.strip('\n').split('\t')
            # To be a known site, hast to match query at least 90% with coverage of 95
            if float(fields[3]) >= 90 and float(fields[10]) >= 95:
                left_pos = int(fields[6])
                right_pos = int(fields[7])
                if left_pos > right_pos:
                    orientation = 'R'
                else:
                    orientation = 'F'
                new_pos = Position(min(left_pos, right_pos), max(left_pos, right_pos))
                new_pos.orientation = orientation
                new_pos.isolate_dict = {ref_name: '+'}
                positions_list.append(new_pos)
    # remove the output file
    run_command(['rm', blast_output], shell=True)

    # return the list of positions and the name of the reference

    return positions_list, ref_name

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
    slice_size = int(largest_value / len(list_of_range_tuples))

    #create our list of boxes
    range_boxes = []
    for i in range(0, len(list_of_range_tuples) + 10):
        range_boxes.append([])
    #populate boxes
    for tup in list_of_range_tuples:
        index_1 = int(tup[0] / slice_size)
        index_2 = int(tup[1] / slice_size)
        while index_1 <= index_2:
            try:
                range_boxes[index_1].append(tup)
            except IndexError:
                return False, False
            index_1 += 1

    # find box for new range to check
    start = min(range_to_check[0], range_to_check[1])
    stop = max(range_to_check[1], range_to_check[0])

    index_start = int(start / slice_size)
    index_stop = int(stop / slice_size)

    # check each potential box
    while index_start <= index_stop and index_stop <= int(((largest_value/slice_size) + 1)):
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
        index_start += 1

    return False, False

def final_ranges_check(positions, gap):

    # if the number of positions to check is only 1, then we don't need
    # to do an overlap check. Exit the function
    if len(positions) <= 1:
        return(positions)
    # want to check if any of the final positions need to be merged
    final_positions = []
    # go through each position, and check if there are any other positions it should be merged with
    # obviously we want to ignore the position itself
    for pos in positions:
        pos_index = positions.index(pos)
        positions_without_check = positions[:pos_index] + positions[(pos_index + 1):]
        matching_position, new_range = check_ranges(positions_without_check, (pos.x, pos.y), gap, pos.orientation)
        # if the current range overlaps with a range we already have (so no False's are returned)
        if matching_position != False:
            # grab the isolate dictionary of the matching position
            matching_isolate_dict = matching_position.isolate_dict
            # grab the isolate dictionry of the current position
            current_isolate_dict = pos.isolate_dict
            # merge the two dictionaries together
            new_isolate_dict = matching_isolate_dict.copy()
            new_isolate_dict.update(current_isolate_dict)
            # Remove the matching position from the list
            positions.remove(matching_position)
            # Create the new position and add it to the list of final positions
            new_pos = Position(new_range[0], new_range[1])
            new_pos.orientation = pos.orientation
            new_pos.isolate_dict = new_isolate_dict
            final_positions.append(new_pos)
        else:
            # just append the position to the list of final positions
            final_positions.append(pos)

    return(final_positions)

def write_output(pos_list, isolate_list, out_prefix, ref_name, imprecise_value, unconfident_value, binary=False):

    if not binary:
        out = open(out_prefix + '_full_compiled.txt', 'w')
        print('Writing full output file to ' + out_prefix + '_full_compiled.txt')
    else:
        out = open(out_prefix + '_binary_compiled.txt', 'w')
        print('Writing binary output file to ' + out_prefix + '_binary_compiled.txt')

    header_line = ['isolate']
    for pos in pos_list:
        if pos.orientation == 'F':
            header_line.append(str(pos.x) + '-' + str(pos.y))
        else:
            header_line.append(str(pos.y) + '-' + str(pos.x))

    out.write('\t'.join(header_line) + '\n')

    # add the reference
    row = [ref_name]
    for pos in pos_list:
        if ref_name in pos.isolate_dict.keys():
            if not binary:
                row.append(pos.isolate_dict[ref_name])
            else:
                row.append('1')
        else:
            if not binary:
                row.append('-')
            else:
                row.append('0')
    out.write('\t'.join(row) + '\n')

    # now loop through each isolate and create each row
    for isolate in isolate_list:
        row = [isolate]
        for pos in pos_list:
            if isolate in pos.isolate_dict.keys():
                if not binary:
                    row.append(pos.isolate_dict[isolate])
                else:
                    if pos.isolate_dict[isolate] == '+':
                        row.append('1')
                    elif pos.isolate_dict[isolate] == '*':
                        row.append(imprecise_value)
                    elif pos.isolate_dict[isolate] == '?':
                        row.append(unconfident_value)
            else:
                if not binary:
                    row.append('-')
                else:
                    row.append('0')
        out.write('\t'.join(row) + '\n')

    # finally, write out the orientation and flanking gene info
    # set up these rows, only for the full output file
    # the binary output file doesn't need this section
    if not binary:
        orientation_row = ['orientation']
        left_locus_row = ['left locus']
        right_locus_row = ['right locus']
        left_distance_row = ['left distance']
        right_distance_row = ['right distance']
        left_strand_row = ['left strand']
        right_strand_row = ['right strand']
        left_product_row = ['left product']
        right_product_row = ['right product']
        for pos in pos_list:
            orientation_row.append(pos.orientation)
            left_locus_row.append(pos.gene_left)
            right_locus_row.append(pos.gene_right)
            left_distance_row.append(str(pos.left_distance))
            right_distance_row.append(str(pos.right_distance))
            left_strand_row.append(str(pos.left_strand))
            right_strand_row.append(str(pos.right_strand))
            left_product_row.append(pos.left_feature.qualifiers['product'][0])
            right_product_row.append(pos.right_feature.qualifiers['product'][0])

        out.write('\t'.join(orientation_row) + '\n')
        out.write('\t'.join(left_locus_row) + '\n')
        out.write('\t'.join(right_locus_row) + '\n')
        out.write('\t'.join(left_distance_row) + '\n')
        out.write('\t'.join(left_strand_row) + '\n')
        out.write('\t'.join(left_product_row) + '\n')
        out.write('\t'.join(right_distance_row) + '\n')
        out.write('\t'.join(right_strand_row) + '\n')
        out.write('\t'.join(right_product_row) + '\n')

    # close the output file
    out.close()

def parse_args():
    parser = argparse.ArgumentParser(description="Create a table of IS hits in all isolates for ISMapper")
    # Inputs
    parser.add_argument('--tables', nargs='+', type=str, required=True, help='tables to compile')
    parser.add_argument('--reference', type=str, required=True,
                        help='gbk file of reference')
    parser.add_argument('--query', type=str, required=True,
                        help='fasta file for insertion sequence query for compilation')
    # Parameters for hits
    parser.add_argument('--gap', type=int, required=False, default=0,
                        help='distance between regions to call overlapping, default is 0')
    # TODO: Check these parameters work in the cases where they are supposed to
    parser.add_argument('--cds', nargs='+', type=str, required=False, default=['locus_tag', 'gene', 'product'],
                        help='qualifiers to look for in reference genbank for CDS features')
    parser.add_argument('--trna', nargs='+', type=str, required=False, default=['locus_tag', 'product'],
                        help='qualifiers to look for in reference genbank for tRNA features')
    parser.add_argument('--rrna', nargs='+', type=str, required=False, default=['locus_tag', 'product'],
                        help='qualifiers to look for in reference genbank for rRNA features')
    parser.add_argument('--imprecise', type=str, required=False, default='1',
                        help='Binary value for imprecise (*) hit (can be 1, 0 or 0.5), default is 1')
    parser.add_argument('--unconfident', type=str, required=False, default='0',
                        help='Binary value for questionable (?) hit (can be 1, 0 or 0.5), default is 0')
    # Output parameters
    parser.add_argument('--out_prefix', type=str, required=True, help='Prefix for output file')

    return parser.parse_args()

def main():

    start_time = time.time()

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))
    list_of_isolates = []

    list_of_positions = []

    reference_fasta = args.reference.split('.g')[0]
    # Create a fasta file of the reference for BLAST
    print('Creating fasta file and database of reference ...')
    gbk_to_fasta(args.reference, reference_fasta)
    # Make a BLAST database
    blast_db(reference_fasta)
    # Get the reference positions and orientations for this IS query
    print('\nGetting query positions in reference ...')
    list_of_positions, ref_name = get_ref_positions(reference_fasta, args.query, list_of_positions)

    elapsed_time = time.time() - start_time
    print('Time taken: ' + str(elapsed_time))

    # Loop through each table given to the tables argument
    print('Collating results files ...')
    for result_file in unique_results_files:
        # Get isolate name
        isolate = os.path.split(result_file)[1].split('__')[0]
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
                            new_pos = Position(is_start, is_end)
                            new_pos.orientation = orientation
                            new_pos.isolate_dict = isolate_dict
                            list_of_positions.append(new_pos)

                        # If the list of positions isn't empty, then there are ranges to check against
                        else:
                            old_position, new_range = check_ranges(list_of_positions, (is_start, is_end), args.gap,
                                                                   orientation)
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
                                new_pos = Position(new_range[0], new_range[1])
                                new_pos.orientation = orientation
                                new_pos.isolate_dict = isolate_dict
                                list_of_positions.append(new_pos)
                            # Otherwise this range hasn't been seen before, so all values are False
                            else:
                                if '?' in call:
                                    isolate_dict[isolate] = '?'
                                elif '*' in call:
                                    isolate_dict[isolate] = '*'
                                else:
                                    isolate_dict[isolate] = '+'
                                new_pos = Position(is_start, is_end)
                                new_pos.orientation = orientation
                                new_pos.isolate_dict = isolate_dict
                                list_of_positions.append(new_pos)

    # do one last check for positions that should be merged
    list_of_positions = final_ranges_check(list_of_positions, args.gap)

    elapsed_time = time.time() - start_time
    print('Time taken: ' + str(elapsed_time))

    # Get the flanking genes for each position now they've all been merged
    print('Getting flanking genes for each position (this step is the longest and could take some time) ...')

    # Get feature list
    gb = SeqIO.read(args.reference, "genbank")
    feature_list = get_features(gb)

    # Get flanking genes
    for pos in list_of_positions:
        pos.get_flanking_genes(gb, feature_list)

    elapsed_time = time.time() - start_time
    print('Time taken: ' + str(elapsed_time))

    # Order positions from smallest to largest for final table output
    list_of_positions.sort(key=lambda x: x.x)
    write_output(list_of_positions, list_of_isolates, args.out_prefix, ref_name, args.imprecise, args.unconfident)
    write_output(list_of_positions, list_of_isolates, args.out_prefix, ref_name, args.imprecise, args.unconfident, binary=True)

    elapsed_time = time.time() - start_time
    print('Table compilation finished in ' + str(elapsed_time))

if __name__ == '__main__':
    main()