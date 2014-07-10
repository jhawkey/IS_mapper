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
    parser.add_argument('--tables', nargs='+', type=str, required=False, help='tables to compile')

    return parser.parse_args()

def check_ranges(ranges, range_to_check, gap, unpaired = None):

    if unpaired == None:
        start = range_to_check[0]
        stop = range_to_check[1]
        for i in range(0, len(ranges)):
            x = min(ranges[i][0], ranges[i][1])
            y = max(ranges[i][0], ranges[i][1])
            if start in range(x-gap, y+1):
                new_start = min(x, start)
                new_end = max(y, stop)
                print 'we are at the first return'
                return ranges[i], (new_start, new_end)
            elif stop in range(x, y+gap+1):
                new_start = min(x, start)
                new_end = max(y, stop)
                print 'we are at the second return'
                return ranges[i], (new_start, new_end)
        print 'we are at the last return'
        return False, False
    elif unpaired == True:
        coord = range_to_check
        for i in range(0, len(ranges)):
            x = min(ranges[i][0], ranges[i][1])
            y = max(ranges[i][0], ranges[i][1])
            if coord in range(x-gap, y+1):
                return ranges[i], True
            elif coord in range(x, y+gap+1):
                return ranges[i], True
            else:
                return False, False


def main():

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))

    list_of_positions = collections.defaultdict(dict) # key1 = pos, key2 = isolate, value = +/-
    unpaired_hits = {}
    for result_file in unique_results_files:
        isolate = result_file.split('__')[0]
        header = 0
        with open(result_file) as file_open:
            for line in file_open:
                if header == 0:
                    header = header + 1
                else:
                    print isolate
                    info = line.strip('\n').split('\t')
                    is_start = int(info[3])
                    print is_start
                    if info[4] != '':
                        print 'this is a paired hit'
                        is_end = int(info[4])
                        print is_end
                    else:
                        print 'this is not a paired hit'
                        is_end = info[4]
                        if isolate not in unpaired_hits:
                            print 'adding it to the unpaired hits'
                            unpaired_hits[isolate] = [is_start]
                        else:
                            print 'adding it to the unpaired hits (else)'
                            unpaired_hits[isolate].append(is_start)
                    if (is_start, is_end) not in list_of_positions and is_end != '':
                        if list_of_positions.keys() != []:
                            old_range, new_range = check_ranges(list_of_positions.keys(), (is_start, is_end), 300, unpaired=None)
                            if old_range != False:
                                store_values = list_of_positions[old_range]
                                del list_of_positions[old_range]
                                list_of_positions[new_range] = store_values)
                                list_of_positions[new_range][isolate] = '+'
                            else:
                                list_of_positions[(is_start, is_end)][isolate] = '+'
                        else:
                            list_of_positions[(is_start, is_end)][isolate] = '+'
                    elif (is_start, is_end) in list_of_positions and is_end != '':
                        list_of_positions[(is_start, is_end)][isolate] = '+'
        
        # dealing with unpaired hits after files have been read in
        paired_hits = list_of_positions.keys()
        for isolate in unpaired_hits:
            for hit in unpaired_hits[isolate]:
                range_hit, boolean = check_ranges(paired_hits, hit, 300, unpaired=True)
                if boolean == True:
                    list_of_positions[range_hit][isolate] = '+*'
                else:
                    pass

    print list_of_positions
    print unpaired_hits


if __name__ == "__main__":
    main()