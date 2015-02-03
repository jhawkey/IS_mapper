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

def parse_args():

    parser = ArgumentParser(description="Create a binary table from the compiled table for ISMapper")
    parser.add_argument('--table', type=str, required=True, help='input table')
    parser.add_argument('--output', type=str, required=True, help='output binary table')
    parser.add_argument('--imprecise', type=str, required=False, default='1', help='Binary value for imprecise (*) hit (can be 1, 0 or 0.5), default is 1')
    parser.add_argument('--question', type=str, required=False, default='0', help='Binary value for questionable (?) hit (can be 1, 0 or 0.5), default is 0')
    parser.add_argument('--delimiter', type=str, required=False, default=',', help='delimiter for output file (default is , can also be t for tab)')

    return parser.parse_args()

def main():

    args = parse_args()

    # Open the output file for writing
    out_file = open(args.output, 'w')
    # Check to make sure the delmiter is correct
    if args.delimiter != ',' and args.delimiter != 't':
        print 'Delimiter type unknown. Must be , or t.'
        # Exit the script if it's not
        sys.exit()

    header_names = {}
    ordered_header_names = []
    # Open up the compiled table and read through line by line
    with open(args.table) as table_in:
        header = 0
        for line in table_in:
            # Get the header information, append to dictionary and list
            if header == 0:
                header_info = line.strip().split('\t')
                for element in header_info:
                    header_names[element] = []
                    # Just want to record positions
                    if 'isolate' not in element:
                        ordered_header_names.append(element)
                # Write the header to the output file
                if args.delimiter == 't':
                    out_file.write(line)
                elif args.delimiter == ',':
                    info = line.strip().split('\t')
                    out_file.write(','.join(info) + '\n')
                # Increment header as we're done with header now
                header += 1
            elif header == 1:
                # Ignore the reference isolate for the time being as it's wrong
                header += 1
            # If we're at flanking genes (end of table), don't include this information
            elif 'left' in line or 'right' in line:
                pass
            # Otherwise add +, *, ? or - values to correct position header
            else:
                info = line.strip().split('\t')
                # Gather list of isolate names
                header_names['isolate'].append(info[0])
                for index in range(1, len(info)):
                    header_names[ordered_header_names[index-1]].append(info[index])
    # For each position, count the occurance of + and *
    # Compare this to occurance of ?, want to make ? a 1 (confident hit)
    # if there are many other isolates with an IS at this position
    for pos in ordered_header_names:
        conf_value = header_names[pos].count('+')
        imp_value = header_names[pos].count('*')
        quest_value = header_names[pos].count('?')
        # Repalce +, *, ? and - with correct values
        if (conf_value + imp_value) != 0 and (float(quest_value) / float((conf_value + imp_value))) < 1:
            # Then mark ? as 1's
            header_names[pos] = [h.replace('+', '1') for h in header_names[pos]]
            header_names[pos] = [h.replace('-', '0') for h in header_names[pos]]
            header_names[pos] = [h.replace('*', '1') for h in header_names[pos]]
            header_names[pos] = [h.replace('?', '1') for h in header_names[pos]]
        else:
            # Mark ? as user specificed value (default 0)
            header_names[pos] = [h.replace('+', '1') for h in header_names[pos]]
            header_names[pos] = [h.replace('-', '0') for h in header_names[pos]]
            header_names[pos] = [h.replace('*', args.imprecise) for h in header_names[pos]]
            header_names[pos] = [h.replace('?', args.question) for h in header_names[pos]]
    # Should now have each position in binary format
    # Can now loop through each row (determined by isolate) and print out
    # Header has already been printed to file
    for isolate_no in range(0, len(header_names['isolate'])):
        # Get isolate name
        row = [header_names['isolate'][isolate_no]]
        for pos in ordered_header_names:
            row.append(header_names[pos][isolate_no])
        # Join together the row with the correct delimiter
        if args.delimiter == 't':
            out_file.write('\t'.join(row) + '\n')
        elif args.delimiter == ',':
            out_file.write(','.join(row) + '\n')
    out_file.close()

if __name__ == "__main__":
    main()