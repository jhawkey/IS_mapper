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

def main():

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))

    list_of_positions = collections.defaultdict(dict) # key1 = pos, key2 = isolate, value = +/-
    for result_file in unique_results_files:
        isolate = result_file.split('__')[0]
        header = 0
        with open(result_file) as file_open:
            for line in file_open:
                if header == 0:
                    header = header + 1
                else:
                    info = line.strip('\n').split('\t')
                    if (info[3], info[4]) not in list_of_positions:
                        list_of_positions[(info[3], info[4])][isolate] = '+'
                    elif (info[3], info[4]) in list_of_positions:
                        list_of_positions[(info[3], info[4])][isolate] = '+'
    
                    #if isolate not in list_of_positions:
                    #    list_of_positions[isolate] = [(info[3], info[4])]
                    #else:
                    #    list_of_positions[isolate].append((info[3], info[4]))
    print list_of_positions


if __name__ == "__main__":
    main()