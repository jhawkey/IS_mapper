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

def check_ranges(ranges, range_to_check, gap):

    start = range_to_check[0]
    if range_to_check[1] == '':
        stop = range_to_check[0]
    else:
        stop = range_to_check[1]
    for i in range(0, len(ranges)):
        x = min(ranges[i][0], ranges[i][1])
        y = max(ranges[i][0], ranges[i][1])
        if start in range(x-gap, y+1):
            new_start = min(x, start)
            if stop != '':
                new_end = max(y, stop)
            else:
                new_end = y
            return ranges[i], (new_start, new_end)
        elif stop in range(x, y+gap+1):
            new_start = min(x, start)
            new_end = max(y, stop)
            return ranges[i], (new_start, new_end)

    return False, False

def main():

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))

    list_of_positions = {} # key1 = pos, key2 = isolate, value = +/-
    for result_file in unique_results_files:
        isolate = result_file.split('__')[0]
        header = 0
        with open(result_file) as file_open:
            for line in file_open:
                if header == 0:
                    header = header + 1
                else:
                    info = line.strip('\n').split('\t')
                    is_start = int(info[3])
                    if info[4] != '':
                        is_end = int(info[4])
                    else:
                        is_end = ''
                    if (is_start, is_end) not in list_of_positions and is_end != '':
                        if list_of_positions.keys() != []:
                            old_range, new_range = check_ranges(list_of_positions.keys(), (is_start, is_end), 1000)
                            print old_range, new_range
                            if old_range != False:
                                store_values = list_of_positions[old_range]
                                del list_of_positions[old_range]
                                list_of_positions[new_range] = store_values
                                list_of_positions[new_range].append(isolate)
                        list_of_positions[(is_start, is_end)] = [isolate]
                    elif (is_start, is_end) in list_of_positions and is_end != '':
                        list_of_positions[(is_start, is_end)].append(isolate)
                    '''elif is_end == '':
                        check_ranges(list_of_positions.keys(), (is_start, is_end), 1000)
                        print old_range, new_range
                        if old_range != False:
                            store_values = list_of_positions[old_range]
                            del list_of_positions[old_range]
                            list_of_positions[new_range] = store_values
                            list_of_positions[new_range].append(isolate)
                        list_of_positions[(is_start, is_end)] = [isolate]'''


    

    print list_of_positions


if __name__ == "__main__":
    main()