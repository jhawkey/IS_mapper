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
    parser.add_argument('--table', type=str, required=True, help='raw compiled table to compress')
    parser.add_argument('--output', type=str, required=True, help='name of output file')

    return parser.parse_args()

def collapse_row(row_dict, row_positions):

    collapsed_rows = {}

    for position in row_positions:
        x = position[0]
        y = position[1]
        xy_range = range(x, y)
        for position_2 in row_positions:
            if position_2[0] != x and position_2[1] != y:
                if position_2[0] in xy_range or position_2[1] in xy_range:
                    new_x = min(position_2[0], x)
                    new_y = max(position_2[1], y)
                    collapsed_rows[(new_x, new_y)] = [row_dict[position], row_dict[position_2]]
    return collapsed_rows

def main():

    args = parse_args()

    f_rows = {}
    r_rows = {}

    with open(args.table) as table_in:
        count = 0
        for line in table_in:
            if count == 0:
                header = line.strip().split('\t')
                count += 1
            else:
                info = line.strip().split('\t')
                if info[0] == 'F':
                    f_rows[(int(info[1]), int(info[2]))] = info[3:]
                elif info[0] == 'R':
                    r_rows[(int(info[1]), int(info[2]))] = info[3:]
    print header
    #print f_rows
    #print r_rows
    order_f_rows = f_rows.keys()
    order_f_rows.sort()
    order_r_rows = r_rows.keys()
    order_r_rows.sort()

    collapsed_f_rows = collapse_row(f_rows, order_f_rows)
    collapsed_r_rows = collapse_row(r_rows, order_r_rows)

    print collapsed_f_rows
    output_file = open(args.output, 'w')
    output_file.write('\t'.join(header) + '\n')

if __name__ == "__main__":
    main()