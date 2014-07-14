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

def parse_args():

    parser = ArgumentParser(description="create a table of features for the is mapping pipeline")
    parser.add_argument('--tables', nargs='+', type=str, required=True, help='tables to compile')
    parser.add_argument('--reference', type=str, required=True, help='reference to determine positions against')
    parser.add_argument('--gap', type=int, required=False, default=300, help='distance between regions to call overlapping')

    return parser.parse_args()

def check_ranges(ranges, range_to_check, gap, orientation, unpaired = None):

    if unpaired == None:
        start = range_to_check[0]
        stop = range_to_check[1]
        #print 'this is the start and stop'
        #print start, stop
        #print 'these are the ranges that we are checking against'
        ranges_list = ranges.keys()
        #print ranges
        #print 'this is the orientation: ' + orientation
        for i in range(0, len(ranges_list)):
            #print 'these are the x and y values'
            x = ranges_list[i][0]
            y = ranges_list[i][1]
            #print x, y
            #print 'this the orientation of the range that we are currently checking against'
            checking_orientation = ranges[(ranges_list[i][0], ranges_list[i][1])]
            #print checking_orientation
            if orientation == checking_orientation:
                #print 'so the orientation matches'
                #print 'these are the x and y values'
                #print x, y
                if orientation == "5' to 3'":
                    if start >= (x - gap) and start <= (y + 1):
                        #print 'we are the start part'
                        new_start = min(x, start)
                        new_end = max(y, stop)
                        #print 'this is the start and stop that were successful'
                        #print start, stop
                        return ranges_list[i], (new_start, new_end), orientation
                    elif stop >= x and stop <= (y + gap + 1):
                        #print 'we are in the stop part'
                        new_start = min(x, start)
                        new_end = max(y, stop)
                        #print 'this is the start and stop that were successful'
                        #print start, stop
                        return ranges_list[i], (new_start, new_end), orientation
                else:
                    if start <= (x + gap) and start > y:
                        new_start = max(x, start)
                        new_end = min(y, stop)
                        return ranges_list[i], (new_start, new_end), orientation
                    elif stop >= (y - gap) and stop < x:
                        new_start = max(x, start)
                        new_end = min(y, stop)
                        return ranges_list[i], (new_start, new_end), orientation
            #print 'the orientation did not match'

        return False, False, False
    elif unpaired == True:
        coord = range_to_check
        for i in range(0, len(ranges)):
            x = min(ranges[i][0], ranges[i][1])
            y = max(ranges[i][0], ranges[i][1])
            if coord in range(x-gap, y+1):
                return ranges[i], True, False
            elif coord in range(x, y+gap+1):
                return ranges[i], True, False
           
        return False, False, False

def get_ref_positions():


def main():

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))
    list_of_isolates = []

    list_of_positions = collections.defaultdict(dict) # key1 = pos, key2 = isolate, value = +/-
    unpaired_hits = {}
    position_orientation = {}

    ref_row = get_ref_positions()

    for result_file in unique_results_files:
        isolate = result_file.split('__')[0]
        list_of_isolates.append(isolate)
        header = 0
        with open(result_file) as file_open:
            for line in file_open:
                if header == 0:
                    header = header + 1
                elif 'No hits found' not in line and line != '':
                    #print isolate
                    info = line.strip('\n').split('\t')
                    orientation = info[1]
                    is_start = int(info[3])
                    #print is_start
                    if info[4] != '':
                        is_end = int(info[4])
                        #print is_end
                    else:
                        #print 'this is not a paired hit'
                        is_end = info[4]
                        if isolate not in unpaired_hits:
                            #print 'adding it to the unpaired hits'
                            unpaired_hits[isolate] = [is_start]
                        else:
                            #print 'adding it to the unpaired hits (else)'
                            unpaired_hits[isolate].append(is_start)
                    if (is_start, is_end) not in list_of_positions and is_end != '':
                        if list_of_positions.keys() != []:
                            old_range, new_range, new_orientation = check_ranges(position_orientation, (is_start, is_end), args.gap, orientation, unpaired=None)
                            #print old_range, new_range, new_orientation
                            if old_range != False:
                                store_values = list_of_positions[old_range]
                                del list_of_positions[old_range]
                                list_of_positions[new_range] = store_values
                                list_of_positions[new_range][isolate] = '+'
                                del position_orientation[old_range]
                                position_orientation[new_range] = new_orientation
                            else:
                                list_of_positions[(is_start, is_end)][isolate] = '+'
                                position_orientation[(is_start, is_end)] = orientation
                        else:
                            list_of_positions[(is_start, is_end)][isolate] = '+'
                            position_orientation[(is_start, is_end)] = orientation
                    elif (is_start, is_end) in list_of_positions and is_end != '':
                        list_of_positions[(is_start, is_end)][isolate] = '+'
                    #print position_orientation
        
    # dealing with unpaired hits after files have been read in
    paired_hits = list_of_positions.keys()
    for isolate in unpaired_hits:
        for hit in unpaired_hits[isolate]:
            range_hit, boolean, orientation = check_ranges(paired_hits, hit, args.gap, orientation=None, unpaired=True)
            #print range_hit, boolean
            if boolean == True:
                list_of_positions[range_hit][isolate] = '+*'
            else:
                list_of_positions[(hit, hit+1)][isolate] = '+?'

    #print list_of_positions
    #print unpaired_hits
    #print list_of_isolates

    # ordering positions from smallest to largest for final table output
    order_position_list = list(OrderedDict.fromkeys(list_of_positions.keys()))
    order_position_list.sort()
    #print order_position_list

    # create header of table
    header = ['isolate']
    for position in order_position_list:
        header.append(str(position[0]) + '-' + str(position[1]))
    print '\t'.join(header)
    
    # create each row
    for isolate in list_of_isolates:
        row = [isolate]
        for position in order_position_list:
            if isolate in list_of_positions[position]:
                row.append(list_of_positions[position][isolate])
            else:
                row.append('-')
        #row.append('\n')
        print '\t'.join(row)

if __name__ == "__main__":
    main()
