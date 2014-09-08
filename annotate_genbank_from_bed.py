from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqFeature
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from argparse import (ArgumentParser, FileType)
import os

def parse_args():

    parser = ArgumentParser(description="annotate a genbank with either BLAST or BED hits")

    # required qsub options
    parser.add_argument('--intersect_bed', type=str, required=False, help='intersect bed file where blocks will be use to annotate.')
    parser.add_argument('--closest_bed', type=str, required=True, help='closest bed file where blocks will be used to annotate.')
    parser.add_argument('--insertion_seq', type=str, required=True, help='insertion sequence reference in fasta format')
    parser.add_argument('--genbank', type=str, required=False, help='original genbank file that features are going to be added to.')
    parser.add_argument('--newfile', type=str, required=True, help='new filename for genbank having features added to.')

    return parser.parse_args()

def createFeature(hits):

    x_L = hits[0]
    y_L = hits[1]
    x_R = hits[2]
    y_R = hits[3]
    quals = {}

    left_location = SeqFeature.FeatureLocation(x_L, y_L)
    right_location = SeqFeature.FeatureLocation(x_R, y_R)
    if x_L < x_R and y_L < y_R:
        #then in forward orientation, set colour to be red
        quals['colour'] = '2'
        quals['orientation'] = 'forward'
    elif x_L > x_R and y_L > y_R:
        #then in reverse orientation, set colour to be yellow
        quals['colour'] = '7'
        quals['orientation'] = 'reverse'

    left_feature = SeqFeature.SeqFeature(left_location, type='left_end', qualifiers=quals)
    right_feature = SeqFeature.SeqFeature(right_location, type='right_end', qualifiers=quals)
    
    return left_feature, right_feature

def parse_bed(bed_file, file_type, seq_length):
    '''
    Parses's bed file
    '''


    line = 0
    hit_no = 1
    hits = {}
    with open(bed_file) as summary:
        for line in summary:
            print line
            info = line.strip().split('\t')
            if file_type == 'closest':
                if info[3] == -1:
                    return hits
                elif int(info[6]) == 0:
                    #this is an overlap, so will be in the intersect file
                    pass
                elif int(info[6]) <= 10 or float(info[6]) / seq_length >= 0.8:
                    hits['hit_' + str(hit_no)] = [int(info[1]), int(info[2]), int(info[4]), int(info[5])]
                    hit_no += 1
    return hits

def insertion_length(insertion):

    sequence = SeqIO.read(insertion, "fasta")
    length = len(sequence.seq)

    return length

def main():

    args = parse_args()
    genbank = SeqIO.read(args.genbank, 'genbank')
    feature_count = 0
    length_insertion = insertion_length(args.insertion_seq)

    if os.stat(args.intersect_bed)[6] != 0:
        blocks_intersect = parse_bed(args.intersect_bed, 'intersect', length_insertion)
        print blocks_intersect
        for region in blocks_intersect:
            left_end, right_end = createFeature(blocks_intersect[region])
            genbank.features.append(left_end)
            genbank.features.append(right_end)
            feature_count += 2
    if os.stat(args.closest_bed)[6] != 0:
        blocks_closest = parse_bed(args.closest_bed, 'closest', length_insertion)
        print blocks_closest
        for region in blocks_closest:
            left_end, right_end = createFeature(blocks_closest[region])
            genbank.features.append(left_end)
            genbank.features.append(right_end)
            feature_count += 2

    SeqIO.write(genbank, args.newfile, 'genbank')
    print("Added " + str(feature_count) + " features to " + args.newfile)

if __name__ == '__main__':
    main()