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
    parser.add_argument('--output', type=str, required=True, help='name of output file')

    return parser.parse_args()

def main():

    args = parse_args()

    # get a list of 
    unique_results_files = list(OrderedDict.fromkeys(args.tables))
    list_of_isolates = []

    list_of_positions = collections.defaultdict(dict)

    for result in unique_results_files:
        isolate = result.split('_table')[0]
        list_of_isolates.append(isolate)
        header = 0
        with open(result) as result_open:
            for line in result_open:
                if header == 0:
                    header += 1
                elif 'No hits found' not in line and line != '':
                    info = line.split('\t')
                    orient = info[1]
                    x = int(info[2])
                    y = int(info[3])
                    gap = info[4]
                    call = info[5]
                    gene_l = info[8]
                    gene_r = info[11]
                    if orient != "3' unpaired" and orient != "5' unpaired":
                        list_of_positions[(orient, x, y, gap, call, gene_l, gene_r)][isolate] = '+'

    order_position_list = list(OrderedDict.fromkeys(list_of_positions.keys()))
    order_position_list.sort(key=operator.itemgetter(1))

    output_header = ['orientation', 'x', 'y', 'gap', 'call', 'left_gene', 'right_gene']
    output_header.extend(list_of_isolates)

    output_file = open(args.output + '.txt', 'w')
    print '\t'.join(output_header)
    output_file.write('\t'.join(str(i) for i in output_header) + '\n')

    for position in order_position_list:
        row = []
        for i in position:
            row.append(str(i))
        
        for isolate in list_of_isolates:
            if isolate in list_of_positions[position]:
                row.append('+')
            else:
                row.append('-')
        output_file.write('\t'.join(str(i) for i in row) + '\n')
        print '\t'.join(row)



if __name__ == "__main__":
    main()