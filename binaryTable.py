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

    parser = ArgumentParser(description="create a binary table from the compiled table for IS mapper")
    parser.add_argument('--table', type=str, required=True, help='input table')
    parser.add_argument('--output', type=str, required=True, help='output binary table')

    return parser.parse_args()

def main():

    args = parse_args()

    with open(args.table) as table_in:
        for line in table_in:
            header = 0
            if header == 0:
                print line
                header = header + 1
            else:
                info = line.split('\t')
                row = []
                for element in info:
                    if element != '+' and element != '+*' and element != '-' and element != '+?':
                        row.append(element)
                    elif element == '-':
                        row.append('0')
                    elif element == '+':
                        row.append('1')
                    elif element == '+*':
                        row.append('2')
                    elif element == '+?':
                        row.append('3')
                    else:
                        DoError('unknown value in line: ' + element)
                print '\t'.join(row)

if __name__ == "__main__":
    main()