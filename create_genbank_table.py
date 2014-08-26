from argparse import (ArgumentParser, FileType)
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import NcbiblastnCommandline
from operator import itemgetter
import os, sys, re, collections, operator
import numpy as np
from collections import OrderedDict

def parse_args():

    parser = ArgumentParser(description="create a table of features for the is mapping pipeline")
    parser.add_argument('--five_bed', type=str, required=True, help='five prime end bed file')
    parser.add_argument('--three_bed', type=str, required=True, help='three prime end bed file')
    parser.add_argument('--assembly', type=str, required=True, help='assembly file for annotation (can be genbank or fasta)')
    parser.add_argument('--type', type=str, required=True, help='the type of assembly file (is either a genbank or fasta)')
    parser.add_argument('--output', type=str, required=True, help='prefix for output file')
    return parser.parse_args()

def create_feature(hit, end):

    start = int(hit[1])
    stop = int(hit[2])
    quals = {}

    location = SeqFeature.FeatureLocation(start, stop)
    if end == 'five':
        quals['colour'] = '2'
        quals['end'] = 'left_end'
        feat_type = 'left end'
    elif end == 'three':
        quals['colour'] = '7'
        quals['end'] = 'right end'
        feat_type = 'right_end'

    feature = SeqFeature.SeqFeature(location, type=feat_type, qualifiers=quals)
    
    return feature

def main():

    args = parse_args()

    header = ['contig', 'end', 'x', 'y']
    results_five = []
    results_three = []
    if os.stat(args.five_bed)[6] == 0 and os.stat(args.three_bed)[6] == 0:
        output = open(args.output + '_table.txt', 'w')
        output.write('\t'.join(header + '\n'))
        output.write('No hits found')
        output.close()
        sys.exit()

    if os.stat(args.five_bed)[6] != 0:
        with open(args.five_bed) as five_bed:
            for line in five_bed:
                info = line.strip().split('\t')
                results_five.append([info[0], info[1], info[2]])
    if os.stat(args.three_bed)[6] != 0:
        with open(args.three_bed) as three_bed:
            for line in three_bed:
                info = line.strip().split('\t')
                results_three.append([info[0], info[1], info[2]])
    if args.type == 'fasta':
        # convert the assembly to a genbank file first
        print('Creating multi entry genbank for annotation...')
        count = SeqIO.convert(args.assembly, 'fasta', args.output + '.gbk', 'genbank', generic_dna)
        print('Successfully converted %i records' % count)
        # open file to write
        annotated_genbank = args.output + '_annotated.gbk'
        handle = open(annotated_genbank, 'w')
        # read in file
        record_list = SeqIO.parse(args.output + '.gbk', 'genbank')
    elif args.type == 'genbank':
        record_list = SeqIO.parse(args.assembly, 'genbank')
    # initalise a list where any edited contigs will go
    new_record_list = []
    # intialise number of features added
    feature_count = 0
    # go through each record and see if there is a five or three end hit
    for record in record_list:
        for result in results_five:
            if result[0] == record.id:
                # then we need to annotate this hit
                new_feature = create_feature(result, 'five')
                record.features.append(new_feature)
                feature_count += 1
                new_record_list.append(record)
            else:
                new_record_list.append(record)
        for result in results_three:
            if result[0] == record.id:
                # then we need to annotate this hit
                new_feature = create_feature(result, 'three')
                record.features.append(new_feature)
                feature_count += 1
                new_record_list.append(record)
            else:
                new_record_list.append(record)
    SeqIO.write(new_record_list, args.output + '_annotated.gbk', 'genbank')
    print('Added ' + str(feature_count) + ' features to ' + args.output + '_annotated.gbk')

    output = open(args.output + '_table.txt', 'w')
    output.write('\t'.join(header + '\n'))
    for result in results_five:
        output.write(result[0] + '\tfive\t' + result[1] + '\t' + result[2] + '\n')
    for result in results_three:
        output.write(result[0] + '\tthree\t' + result[1] + '\t' + result[2] + '\n')
    output.close()

if __name__ == '__main__':
    main()
