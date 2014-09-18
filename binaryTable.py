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
    parser.add_argument('--delimiter', type=str, required=False, default=',', help='delimiter for output file (default is , can also be t for tab)')

    return parser.parse_args()

def main():

    args = parse_args()

    out_file = open(args.output, 'w')

    if args.delimiter != ',' and args.delimiter != 't':
        print 'Delimiter type unknown. Must be , or t.'
        sys.exit()

    with open(args.table) as table_in:
        header = 0
        for line in table_in:
            if header == 0:
                if args.delimiter == 't':
                    out_file.write(line)
                elif args.delimiter == ',':
                    info = line.strip().split('\t')
                    out_file.write(','.join(info) + '\n')
                header += 1
            elif 'flanking genes' in line:
                pass
            else:
                info = line.strip().split('\t')
                name = info[0]
                if '.fasta' in name:
                    name = name.split('.fasta')[0]
                else:
                    name = name.split('_table.txt')[0]
                row = [name]
                for element in info[1:]:
                    if element == '+':
                        row.append('1')
                    elif element == '-':
                        row.append('0')
                    else:
                        print element
                        DoError('unknown value in line')
                if args.delimiter == 't':
                    out_file.write('\t'.join(row) + '\n')
                elif args.delimiter == ',':
                    out_file.write(','.join(row) + '\n')
    out_file.close()
if __name__ == "__main__":
    main()