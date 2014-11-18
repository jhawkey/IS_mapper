#!/usr/bin/env python

# read in genbank file, print out coordinates & strand of features
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
from compiled_table import get_flanking_genes, get_other_gene, get_qualifiers

def parse_args():

    parser = ArgumentParser(description="create a table of features for ISMapper")
    parser.add_argument('--intersect_bed', type=str, required=True, help='intersection bed file')
    parser.add_argument('--closest_bed', type=str, required=True, help='closestBed bed file')
    parser.add_argument('--left_bed', type=str, required=True, help='merged bed file for left end (5)')
    parser.add_argument('--right_bed', type=str, required=True, help='merged bed file for right end (3)')
    parser.add_argument('--reference_genbank', type=str, required=True, help='reference genbank file to find flanking genes of regions')
    parser.add_argument('--insertion_seq', type=str, required=True, help='insertion sequence reference in fasta format')
    parser.add_argument('--cds', nargs='+', type=str, required=False, default=['locus_tag', 'gene', 'product'], help='qualifiers to look for in reference genbank for CDS features (default locus_tag gene product)')
    parser.add_argument('--trna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for tRNA features (default locus_tag product)')
    parser.add_argument('--rrna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for rRNA features (default locus_tag product)')
    parser.add_argument('--min_range', type=float, required=False, default=0.5, help='Minimum percent size of the gap to be called a known hit (default 0.5, or 50 percent)')
    parser.add_argument('--max_range', type=float, required=False, default=1.5, help='Maximum percent size of the gap to be called a known hit (default 1.5, or 150 percent)')
    parser.add_argument('--temp_folder', type=str, required=True, help='location of temp folder to place intermediate blast files in')
    parser.add_argument('--output', type=str, required=True, help='name for output file')
    return parser.parse_args()

def insertion_length(insertion):

    sequence = SeqIO.read(insertion, "fasta")
    length = len(sequence.seq)

    return length

def doBlast(blast_input, blast_output, database):
    #perform BLAST
    blastn_cline = NcbiblastnCommandline(query=blast_input, db=database, outfmt="'6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs'", out=blast_output)
    stdout, stderr = blastn_cline()

def check_seq_between(gb, insertion, start, end, name, temp):

    genbank = SeqIO.read(gb, 'genbank')
    seq_between = genbank.seq[start:end]
    seq_between = SeqRecord(Seq(str(seq_between), generic_dna), id=name)
    print name
    print len(seq_between)
    SeqIO.write(seq_between, temp + name + '.fasta', 'fasta')
    doBlast(temp + name + '.fasta', temp + name + '_out.txt', insertion)
    first_result = 0
    with open(temp + name + '_out.txt') as summary:
        for line in summary:
            if first_result == 0:
                info = line.strip().split('\t')
                coverage = float(info[4]) / float(info[5]) * 100
                #coverage = info[-1]
                hit = [info[3], coverage]
                first_result += 1
            #os.system('rm ' + name + '.fasta ' + name + '_out.txt')
            return hit
    #os.system('rm ' + name + '.fasta ' + name + '_out.txt')
    hit = []
    return []

def createFeature(hits, orient):

    x_L = hits[0]
    y_L = hits[1]
    x_R = hits[2]
    y_R = hits[3]
    quals = {}

    left_location = SeqFeature.FeatureLocation(x_L, y_L)
    right_location = SeqFeature.FeatureLocation(x_R, y_R)
    if orient == 'F':
        #then in forward orientation, set colour to be red
        quals['colour'] = '2'
        quals['orientation'] = 'forward'
    elif orient == 'R':
        #then in reverse orientation, set colour to be yellow
        quals['colour'] = '7'
        quals['orientation'] = 'reverse'

    left_feature = SeqFeature.SeqFeature(left_location, type='left_end', qualifiers=quals)
    right_feature = SeqFeature.SeqFeature(right_location, type='right_end', qualifiers=quals)

    return left_feature, right_feature

def main():

    args = parse_args()

    results = {}
    removed_results = {}
    region = 1
    lines = 0
    header = ["region", "orientation", "x", "y", "gap", "call", "%ID", "%Cov", "left_gene", "left_strand", "left_distance", "right_gene", "right_strand", "right_distance", "functional_prediction"]
    if os.stat(args.intersect_bed)[6] == 0 and os.stat(args.closest_bed)[6] == 0:
        output = open(args.output + '_table.txt', 'w')
        output.write('\t'.join(header) + '\n')
        output.write('No hits found')
        output.close()
        sys.exit()

    genbank = SeqIO.read(args.reference_genbank, 'genbank')
    feature_count = 0

    intersect_bed_lines = []
    closest_bed_lines = []

    if os.stat(args.intersect_bed)[6] != 0:
        with open(args.intersect_bed) as bed_merged:
            for line in bed_merged:
                intersect_bed_lines.append(line)
                info = line.strip().split('\t')
                
                #set up coordinates for checking: L is the left end of the IS (5') and R is the right end of the IS (3')
                #eg x_L and y_L are the x and y coordinates of the bed block that matches to the region which is flanking the left end or 5' of the IS
                x_L = int(info[1])
                y_L = int(info[2])
                x_R = int(info[4])
                y_R = int(info[5])
                
                #check to see if the gap is reasonable
                if int(info[6]) <= 15:
                    if x_L < x_R or y_L < y_R:
                        orient = 'F'
                        x = x_R
                        y = y_L
                    elif x_L > x_R or y_L > y_R:
                        orient = 'R'
                        x = x_L
                        y = y_R
                    else:
                        print 'neither if statement were correct'

                    left_feature, right_feature = createFeature([x_L, y_L, x_R, y_R], orient)
                    genbank.features.append(left_feature)
                    genbank.features.append(right_feature)
                    feature_count += 2

                    gene_left, gene_right = get_flanking_genes(args.reference_genbank, x, y, args.cds, args.trna, args.rrna)
                    if gene_left[1] == gene_right[1]:
                        funct_pred = 'Gene interrupted'
                    else:
                        funct_pred = ''
                    results['region_' + str(region)] = [orient, str(x), str(y), info[6], 'Novel', '', '', gene_left[-1][:-1], gene_left[-1][-1], gene_left[1], gene_right[-1][:-1], gene_right[-1][-1], gene_right[1], funct_pred]
                    region += 1
                else:
                    removed_results['region_' + str(lines)] = line.strip() + '\tintersect.bed\n'
                lines += 1
    
    is_length = insertion_length(args.insertion_seq)
    with open(args.closest_bed) as bed_closest:
        for line in bed_closest:
            closest_bed_lines.append(line)
            info = line.strip().split('\t')
            
            # then there are no closest regions, this is a dud file
            if info[3] == '-1':
                output = open(args.output, 'w')
                output.write('\t'.join(header) + '\n')
                output.write('No hits found')
                output.close()
                sys.exit()

            x_L = int(info[1])
            y_L = int(info[2])
            x_R = int(info[4])
            y_R = int(info[5])
            if x_L < x_R and y_L < y_R:
                orient = 'F'
                x = x_R
                y = y_L
            elif x_L > x_R and y_L > y_R:
                orient = 'R'
                x = x_L
                y = y_R
            #this is an overlap, so will be in the intersect file
            if int(info[6]) == 0:
                pass
            #this is probably a novel hit where there was no overlap detected
            elif int(info[6]) <= 10:

                left_feature, right_feature = createFeature([x_L, y_L, x_R, y_R], orient)
                genbank.features.append(left_feature)
                genbank.features.append(right_feature)
                feature_count += 2
                
                gene_left, gene_right = get_flanking_genes(args.reference_genbank, x, y, args.cds, args.trna, args.rrna)
                if gene_left[:-1] == gene_right[:-1]:
                    funct_pred = 'Gene interrupted'
                else:
                    funct_pred = ''
                results['region_' + str(region)] = [orient, str(x), str(y), info[6], 'Novel', '', '', gene_left[-1][:-1], gene_left[-1][-1], gene_left[1], gene_right[-1][:-1], gene_right[-1][-1], gene_right[1], funct_pred]
                region += 1
            #this is probably a known hit, but need to check with BLAST
            elif float(info[6]) / is_length >= args.min_range and float(info[6]) / is_length <= args.max_range:
                
                left_feature, right_feature = createFeature([x_L, y_L, x_R, y_R], orient)
                genbank.features.append(left_feature)
                genbank.features.append(right_feature)
                feature_count += 2

                seq_results = check_seq_between(args.reference_genbank, args.insertion_seq, start, end, 'region_' + str(region), args.temp_folder)
                if len(seq_results) != 0 and seq_results[0] >= 80 and seq_results[1] >= 80:
                    #then this is definitely a known site
                    gene_left = get_other_gene(args.reference_genbank, min(start, end), "left", args.cds, args.trna, args.rrna)
                    gene_right = get_other_gene(args.reference_genbank, max(start, end), "right", args.cds, args.trna, args.rrna)
                    #gene_left, gene_right = get_flanking_genes(args.reference_genbank, start, end, args.cds, args.trna, args.rrna)
                    results['region_' + str(region)] = [orient, str(start), str(end), info[6], 'Known', str(seq_results[0]), str('%.2f' % seq_results[1]), gene_left[-1][:-1], gene_left[-1][-1], gene_left[1], gene_right[-1][:-1], gene_right[-1][-1], gene_right[1]]
                else:
                   #then I'm not sure what this is
                   print 'not sure'
                   gene_left, gene_right = get_flanking_genes(args.reference_genbank, start, end, args.cds, args.trna, args.rrna)
                   if len(seq_results) !=0:
                       results['region_' + str(region)] = [orient, str(start), str(end), info[6], 'Possible related IS', str(seq_results[0]), str('%.2f' % seq_results[1]), gene_left[-1][:-1], gene_left[-1][-1], gene_left[1], gene_right[-1][:-1], gene_right[-1][-1], gene_right[1]]
                   else:
                        removed_results['region_' + str(region)] = line.strip() + '\tclosest.bed\n'                
                        region += 1
            #could possibly be a novel hit but the gap size is too large
            elif float(info[6]) / is_length <= args.min_range and float(info[6]) / is_length < args.max_range:

                left_feature, right_feature = createFeature([x_L, y_L, x_R, y_R], orient)
                genbank.features.append(left_feature)
                genbank.features.append(right_feature)
                feature_count += 2

                gene_left, gene_right = get_flanking_genes(args.reference_genbank, x, y, args.cds, args.trna, args.rrna)
                if gene_left[:-1] == gene_right[:-1]:
                    funct_pred = 'Gene interrupted'
                else:
                    funct_pred = ''
                results['region_' + str(region)] = [orient, str(x), str(y), info[6], 'Novel*', '', '', gene_left[-1][:-1], gene_left[-1][-1], gene_left[1], gene_right[-1][:-1], gene_right[-1][-1], gene_right[1], funct_pred]
                region +=1
            #this is something else altogether - either the gap is really large or something, place it in removed_results
            else:
                removed_results['region_' + str(region)] = line.strip() + '\tclosest.bed\n'
                region += 1

    #looking for unpaired hits which are not in the merged/closest bed files
    #possibly unpaired due to a repeat on one end of the IS


    #sort regions into the correct order
    table_keys = []
    for key in results:
        table_keys.append(key)
    region_indexes = []
    for region in table_keys:
        region_indexes.append(region.split('region_')[1])
    arr = np.vstack((table_keys, region_indexes)).transpose()
    if arr != 0:
        sorted_keys = arr[arr[:,1].astype('int').argsort()]

    #write out the found hits to file
    output = open(args.output + '_table.txt', 'w')
    output.write('\t'.join(header) + '\n')
    if arr != 0:
        for key in sorted_keys[:,0]:
            output.write(key + '\t' + '\t'.join(str(i) for i in results[key]) + '\n')
    if arr == 0:
        output.write('No hits found.')
    output.close()

    #write out hits that were removed for whatever reason to file
    if len(removed_results) != 0:
        output_removed = open(args.output + '_removedHits.txt', 'w')
        for region in removed_results:
            output_removed.write(removed_results[region])
        output_removed.close()

    SeqIO.write(genbank, args.output + '_annotated.gbk', 'genbank')
    print('Added ' + str(feature_count) + ' features to ' + args.output + '_annotated.gbk')

    #return(lines, len(removed_results))

if __name__ == "__main__":
    main()