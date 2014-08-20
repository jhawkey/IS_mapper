from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqFeature
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from argparse import (ArgumentParser, FileType)

def parse_args():

    parser = ArgumentParser(description="annotate a genbank with either BLAST or BED hits")

    # required qsub options
    parser.add_argument('--bed', type=str, required=False, help='intersect bed file where blocks will be use to annotate.')
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
    elif x_L > x_R and y_L > y_R:
        #then in reverse orientation, set colour to be yellow
        quals['colour'] = '7'

    left_feature = SeqFeature.SeqFeature(left_location, type='left_end', qualifiers=quals)
    right_feature = SeqFeature.SeqFeature(right_location, type='right_end', qualifiers=quals)
    
    return left_feature, right_feature

def parse_bed(bed_file):

    line = 0
    hit_no = 1
    hits = {}
    with open(bed_file) as summary:
        for line in summary:
            info = line.strip().split('\t')
            hits['hit_' + str(hit_no)] = [int(info[1]), int(info[2]), int(info[4]), int(info[5])]
            hit_no += 1
    return hits

def main():

    args = parse_args()

    blocks = parse_bed(args.bed)
    feature_count = 0

    genbank = SeqIO.read(args.genbank, 'genbank')
    for region in blocks:
        left_end, right_end = createFeature(blocks[region])
        genbank.features.append(left_end)
        genbank.features.append(right_end)
        feature_count += 2

    SeqIO.write(genbank, args.newfile, 'genbank')
    print("Added " + str(feature_count) + " features to " + args.newfile)

if __name__ == '__main__':
    main()