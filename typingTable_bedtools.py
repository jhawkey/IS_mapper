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

def parse_args():

    parser = ArgumentParser(description="create a table of features for the is mapping pipeline")
    parser.add_argument('--intersect_bed', type=str, required=True, help='intersection bed file')
    parser.add_argument('--closest_bed', type=str, required=True, help='closestBed bed file')
    parser.add_argument('--reference_genbank', type=str, required=True, help='reference genbank file to find flanking genes of regions')
    parser.add_argument('--insertion_seq', type=str, required=True, help='insertion sequence reference in fasta format')
    parser.add_argument('--cds', type=str, required=False, default='locus_tag,gene,product', help='qualifiers to look for in reference genbank for CDS features')
    parser.add_argument('--trna', type=str, required=False, default='locus_tag,product', help='qualifiers to look for in reference genbank for tRNA features')
    parser.add_argument('--rrna', type=str, required=False, default='locus_tag,product', help='qualifiers to look for in reference genbank for rRNA features')
    parser.add_argument('--output', type=str, required=True, help='name for output file')
    return parser.parse_args()

def get_qualifiers(qualifier_list, feature):
    
    return_quals = []
    for qual in qualifier_list:
        try:
            return_quals.append(feature.qualifiers[qual][0])
        except KeyError:
            pass
    return return_quals

def get_flanking_genes(reference, pos_x, pos_y, cds_quals, trna_quals, rrna_quals):

    gb = SeqIO.read(reference, 'genbank')
    distance_l = {}
    distance_r = {}
    cds_features = cds_quals.split(',')
    trna_features = trna_quals.split(',')
    rrna_features = rrna_quals.split(',')

    # cycle through features in genbank
    for feature in gb.features:
        # only if the feature is a CDS, tRNA or rRNA do we care
        if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
            # if both positions inside gene, then just return this information (no other checking required)
            if pos_x in feature.location and pos_y in feature.location:
                # get qualifiers of interest for this feature
                if feature.type == 'CDS':
                    values = get_qualifiers(cds_features, feature)
                elif feature.type == 'tRNA':
                    values = get_qualifiers(trna_features, feature)
                elif feature.type == 'rRNA':
                    values = get_qualifiers(rrna_features, feature)
                # get the strand
                values.append(feature.strand) 
                gene_left = values
                gene_left_distance = abs(feature.location.start - pos_x)
                gene_right = values
                gene_right_distance = abs(feature.location.start - pos_y)
                # exit the function 
                return gene_left, gene_left_distance, gene_right, gene_right_distance
            
            # if x inside gene but y is not, need to report gene x is in gene further down genome from y
            elif pos_x in feature.location and pos_y not in feature.location:
                if feature.type == 'CDS':
                    values = get_qualifiers(cds_features, feature)
                elif feature.type == 'tRNA':
                    values = get_qualifiers(trna_features, feature)
                elif feature.type == 'rRNA':
                    values = get_qualifiers(rrna_features, feature)
                # get the strand
                values.append(feature.strand) 
                gene_left = values
                gene_left_distance = feature.location.start - pos_x
                # then go and find the next closest gene that y is next to
                gene_right, gene_right_distance = get_other_gene(gb.features, pos_y, cds_features, trna_features, rrna_features, 'downstream')
                # exit function
                return gene_left, gene_left_distance, gene_right, gene_right_distance

            elif pos_y in feature.location and pos_x not in feature.location:
                if feature.type == 'CDS':
                    values = get_qualifiers(cds_features, feature)
                elif feature.type == 'tRNA':
                    values = get_qualifiers(trna_features, feature)
                elif feature.type == 'rRNA':
                    values = get_qualifiers(rrna_features, feature)
                # get the strand
                values.append(feature.strand)
                gene_right = values
                gene_right_distance = feature.location.start - pos_y
                # then go and find next closest gene that x is next to
                gene_left, gene_left_distance = get_other_gene(gb.features, pos_x, cds_features, trna_features, rrna_features, 'upstream')
                # exit function
                return gene_left, gene_left_distance, gene_right, gene_right_distance

            # if x and y both aren't inside the gene, find distances from gene
            dist_l = abs(feature.location.start - pos_x)
            dist_r = abs(feature.location.start - pos_y)
            # get qualifiers of interest for feature
            if feature.type == 'CDS':
                values = get_qualifiers(cds_features, feature)
            elif feature.type == 'tRNA':
                values = get_qualifiers(trna_features, feature)
            elif feature.type == 'rRNA':
                values = get_qualifiers(rrna_features, feature)
            # get the strand
            values.append(feature.strand) 
            # append to respective dictionaries
            distance_l[dist_l] = values
            distance_r[dist_r] = values

    # for each side, get the ordered list of distances
    distance_lkeys = list(OrderedDict.fromkeys(distance_l))
    distance_rkeys = list(OrderedDict.fromkeys(distance_r))
    
    gene_left = distance_l[min(distance_lkeys)] 
    gene_right = distance_r[min(distance_rkeys)]
    gene_left_distance = min(distance_lkeys)
    gene_right_distance = min(distance_rkeys)

    # if the genes are equal, we need to fix this (as the gene is not interrupted)
    if gene_left == gene_right:
        # we found the correct left, but need to find the next closest right
        if gene_left_distance < gene_right_distance:
            gene_right, gene_right_distance = get_other_gene(gb.features, pos_y, cds_features, trna_features, rrna_features, 'downstream')
        # we found the correct right, but need to find the next closest left
        elif gene_left_distance > gene_right_distance:
            gene_left, gene_left_distance = get_other_gene(gb.features, pos_x, cds_features, trna_features, rrna_features, 'upstream')

    return gene_left, gene_left_distance, gene_right, gene_right_distance

def get_other_gene(features, pos, cds_features, trna_features, rrna_features, direction):

    distance = {}
    for feature in features:
        if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
            if pos in feature.location:
                values = []
                if feature.type == 'CDS':
                    values = get_qualifiers(cds_features, feature)
                elif feature.type == 'tRNA':
                    values = get_qualifiers(trna_features, feature)
                elif feature.type == 'rRNA':
                    values = get_qualifiers(rrna_features, feature)
                # get the strand
                values.append(feature.strand)
                gene_distance = feature.location.start - pos
                return values, gene_distance
            else:
                dist = feature.location.start - pos
                if dist > 0 and direction == 'downstream':
                    if feature.type == 'CDS':
                        values = get_qualifiers(cds_features, feature)
                    elif feature.type == 'tRNA':
                        values = get_qualifiers(trna_features, feature)
                    elif feature.type == 'rRNA':
                        values = get_qualifiers(rrna_features, feature)
                    # get the strand
                    values.append(feature.strand)
                    distance[dist] = values
                elif dist < 0 and direction == 'upstream':
                    dist = abs(dist)
                    if feature.type == 'CDS':
                        values = get_qualifiers(cds_features, feature)
                    elif feature.type == 'tRNA':
                        values = get_qualifiers(trna_features, feature)
                    elif feature.type == 'rRNA':
                        values = get_qualifiers(rrna_features, feature)
                    # get the strand
                    values.append(feature.strand)
                    distance[dist] = values

    distance_keys = list(OrderedDict.fromkeys(distance))
    gene = distance[min(distance_keys)]
    gene_distance = min(distance_keys)
    #if direction == 'upstream':
    #    gene_distance = -gene_distance
    return gene, gene_distance

def insertion_length(insertion):

    sequence = SeqIO.read(insertion, "fasta")
    length = len(sequence.seq)

    return length

def doBlast(blast_input, blast_output, database):
    #perform BLAST
    blastn_cline = NcbiblastnCommandline(query=blast_input, db=database, outfmt="'6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs'", out=blast_output)
    stdout, stderr = blastn_cline()

def check_seq_between(gb, insertion, start, end):

    genbank = SeqIO.read(gb, 'genbank')
    seq_between = genbank.seq[start:end]
    seq_between = SeqRecord(Seq(str(seq_between), generic_dna), id='temp')
    SeqIO.write(seq_between, 'temp.fasta', 'fasta')
    doBlast('temp.fasta', 'temp_out.txt', insertion)
    first_result = 0
    with open('temp_out.txt') as summary:
        for line in summary:
            if first_result == 0:
                info = line.strip().split('\t')
                coverage = (float(info[1])/float(info[4])) * 100
                hit = [info[3], coverage]
                first_result += 1
            os.system('rm temp.fasta temp_out.txt')
            return hit
    os.system('rm temp.fasta temp_out.txt')
    hit = []
    return []

def main():

    args = parse_args()

    results = {}
    removed_results = {}
    region = 1
    lines = 0
    header = ["region", "orientation", "x", "y", "gap", "call", "%ID", "%Cov", "left_gene", "left_strand", "left_distance", "right_gene", "right_strand", "right_distance", "functional_prediction"]
    if os.stat(args.intersect_bed) != 0:
        with open(args.intersect_bed) as bed_merged:
            for line in bed_merged:
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

                    gene_left, gene_left_dist, gene_right, gene_right_dist = get_flanking_genes(args.reference_genbank, x, y, args.cds, args.trna, args.rrna)
                    if gene_left[:-1] == gene_right[:-1]:
                        funct_pred = 'Gene interrupted'
                    else:
                        funct_pred = ''
                    results['region_' + str(region)] = [orient, str(x), str(y), info[6], 'Novel', '', '', gene_left[:-1], gene_left[-1], gene_left_dist, gene_right[:-1], gene_right[-1], gene_right_dist, funct_pred]
                    region += 1
                else:
                    removed_results['region_' + str(lines)] = line.strip() + '\tintersect.bed\n'
                lines += 1
    
    is_length = insertion_length(args.insertion_seq)
    with open(args.closest_bed) as bed_closest:
        for line in bed_closest:
            info = line.strip().split('\t')
            if int(info[6]) == 0:
                #this is an overlap, so will be in the intersect file
                pass
            elif int(info[6]) <= 10:
                #this is probably a novel hit where there was no overlap detected
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
                gene_left, gene_left_dist, gene_right, gene_right_dist = get_flanking_genes(args.reference_genbank, x, y, args.cds, args.trna, args.rrna)
                if gene_left[:-1] == gene_right[:-1]:
                    funct_pred = 'Gene interrupted'
                else:
                    funct_pred = ''
                results['region_' + str(region)] = [orient, str(x), str(y), info[6], 'Novel', '', '', gene_left[:-1], gene_left[-1], gene_left_dist, gene_right[:-1], gene_right[-1], gene_right_dist, funct_pred]
                region += 1
            elif float(info[6]) / is_length >= 0.8 and float(info[6]) / is_length <= 1.5:
                #this is probably a known hit, but need to check with BLAST
                y_L = int(info[2])
                x_R = int(info[4])
                if y_L < x_R:
                    start = y_L
                    end = x_R
                    orient = 'F'
                else:
                    start = x_R
                    end = y_L
                    orient = 'R'
                print orient
                seq_results = check_seq_between(args.reference_genbank, args.insertion_seq, start, end)
                if len(seq_results) != 0 and seq_results[0] >= 80 and seq_results[1] >= 80:
                    #then this is definitely a known site
                    gene_left, gene_left_dist, gene_right, gene_right_dist = get_flanking_genes(args.reference_genbank, start, end, args.cds, args.trna, args.rrna)
                    results['region_' + str(region)] = [orient, str(start), str(end), info[6], 'Known', str(seq_results[0]), str('%.2f' % seq_results[1]), gene_left[:-1], gene_left[-1], gene_left_dist, gene_right[:-1], gene_right[-1], gene_right_dist]
                else:
                   #then I'm not sure what this is
                   print 'not sure'
                   gene_left, gene_left_dist, gene_right, gene_right_dist = get_flanking_genes(args.reference_genbank, start, end, args.cds, args.trna, args.rrna)
                   if len(seq_results) !=0:
                        results['region_' + str(region)] = [orient, str(start), str(end), info[6], 'Unknown', str(results[0]), str('%.2f' % results[1]), gene_left[:-1], gene_left[-1], gene_left_dist, gene_right[:-1], gene_right[-1], gene_right_dist]
                    else:
                        results['region_' + str(region)] = [orient, str(start), str(end), info[6], 'Unknown', 'no hit', 'no hit', gene_left[:-1], gene_left[-1], gene_left_dist, gene_right[:-1], gene_right[-1], gene_right_dist]
                region += 1
            else:
                #this is something else altogether - either the gap is really large or something, place it in removed_results
                removed_results['region_' + str(region)] = line.strip() + '\tclosest.bed\n'
                region += 1

    #sort regions into the correct order
    table_keys = []
    for key in results:
        table_keys.append(key)
    region_indexes = []
    for region in table_keys:
        region_indexes.append(region.split('region_')[1])
    arr = np.vstack((table_keys, region_indexes)).transpose()
    sorted_keys = arr[arr[:,1].astype('int').argsort()]

    #write out the found hits to file
    output = open(args.output, 'w')
    output.write('\t'.join(header) + '\n')
    for key in sorted_keys[:,0]:
        output.write(key + '\t' + '\t'.join(str(i) for i in results[key]) + '\n')
    output.close()

    #write out hits that were removed for whatever reason to file
    if len(removed_results) != 0:
        output_removed = open(args.output + '_removedHits.txt', 'w')
        for region in removed_results:
            output_removed.write(removed_results[region])
        output_removed.close()

    #return(lines, len(removed_results))

if __name__ == "__main__":
    main()