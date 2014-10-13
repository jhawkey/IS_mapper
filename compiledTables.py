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
from gbkToFasta import gbk_to_fasta

def parse_args():

    parser = ArgumentParser(description="create a table of features for the is mapping pipeline")
    parser.add_argument('--tables', nargs='+', type=str, required=True, help='tables to compile')
    parser.add_argument('--reference_gbk', type=str, required=True, help='gbk file of reference to report closest genes')
    parser.add_argument('--seq', type=str, required=True, help='fasta file for insertion sequence looking for in reference')
    parser.add_argument('--gap', type=int, required=False, default=0, help='distance between regions to call overlapping')
    parser.add_argument('--cds', type=str, required=False, default='gene,product', help='qualifiers to look for in reference genbank for CDS features')
    parser.add_argument('--trna', type=str, required=False, default='product', help='qualifiers to look for in reference genbank for tRNA features')
    parser.add_argument('--rrna', type=str, required=False, default='product', help='qualifiers to look for in reference genbank for rRNA features')
    parser.add_argument('--output', type=str, required=True, help='name of output file')

    return parser.parse_args()

def check_ranges(ranges, range_to_check, gap, orientation):

    start = range_to_check[0]
    stop = range_to_check[1]

    range_list = ranges.keys()

    for i in range(0, len(range_list)):
        if orientation == ranges[(range_list[i][0], range_list[i][1])]:
            x = range_list[i][0]
            y = range_list[i][1]
            if orientation == 'F':
                if x in range(start - gap, stop + 1) or x in range(start, stop + gap + 1):
                    #these ranges overlap
                    new_start = min(x, start)
                    new_end = max(y, stop)
                    return range_list[i], (new_start, new_end), orientation
                elif y in range(start - gap, stop + 1) or y in range(start, stop + gap + 1):
                    #these ranges also overlap
                    new_start = min(x, start)
                    new_end = min(y, stop)
                    return range_list[i], (new_start, new_end), orientation
            elif orientation == 'R':
                if x in range(start - gap, stop + 1) or x in range(start, stop + gap + 1):
                    #these ranges overlap
                    new_start = min(x, start)
                    new_end = max(y, stop)
                    return range_list[i], (new_start, new_end), orientation
                elif y in range(start - gap, stop + 1) or y in range(start, stop + gap + 1):
                    #these ranges also overlap
                    new_start = min(x, start)
                    new_end = min(y, stop)
                    return range_list[i], (new_start, new_end), orientation
    return False, False, False

def get_ref_positions(reference, is_query, positions_dict, orientation_dict):

    is_name = os.path.split(is_query)[1]
    ref_name = os.path.split(reference)[1]
    blast_output = os.getcwd() + '/' + is_name + '_' + ref_name + '.tmp'

    if not os.path.exists(reference):
        os.system('makeblastdb -in ' + reference + ' -dbtype nucl')
    blastn_cline = NcbiblastnCommandline(query=is_query, db=reference, outfmt="'6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs'", out=blast_output)
    stdout, stderr = blastn_cline()
    with open(blast_output) as out:
        for line in out:
            info = line.strip('\n').split('\t')
            if float(info[3]) >= 90 and float(info[4])/float(info[1]) * 100 >= 95:
                positions_dict[(int(info[6]), int(info[7]))][ref_name] = '+'
                if int(info[6]) > int(info[7]):
                    orientation_dict[(int(info[6]), int(info[7]))] = 'R'
                else:
                    orientation_dict[(int(info[6]), int(info[7]))] = 'F'

    return positions_dict, orientation_dict, ref_name

def get_qualifiers(cds_qualifiers, trna_qualifiers, rrna_qualifiers, feature):
    '''
    Takes a list of possible qualifier IDs and attempts
    to find them in the feature given.
    If the qualifier is present, appends to a list, otherwise
    just keeps going.
    '''
    
    return_quals = []
    if feature.type == 'CDS':
        qualifier_list = cds_qualifiers
    elif feature.type == 'tRNA':
        qualifier_list = trna_qualifiers
    elif feature.type == 'rRNA':
        qualifier_list = rrna_qualifiers
    for qual in qualifier_list:
        try:
            return_quals.append(feature.qualifiers[qual][0])
        except KeyError:
            pass
    return return_quals

def get_flanking_genes(reference, left, right, cds_quals, trna_quals, rrna_quals):

    gb = SeqIO.read(reference, 'genbank')
    pos_gene_left = []
    pos_gene_right = []
    distance_with_left = {}
    distance_with_right = {}

    cds_features = cds_quals.split(',')
    trna_features = trna_quals.split(',')
    rrna_features = rrna_quals.split(',')

    for feature in gb.features:
        if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
            values = get_qualifiers(cds_features, trna_features, rrna_features, feature)
            values.append(feature.strand)
            #first check to see if both coordinates fit into the feature
            if left in feature.location and right in feature.location:
                #we want the absolute value because a value with no sign in the compiled table
                #indicates that the gene is interrupted
                gene = [feature.qualifiers[cds_features[0]][0], str(abs(feature.location.start - left)), values]
                pos_gene_left = gene
                pos_gene_right = gene
                return pos_gene_left, pos_gene_right
            elif left in feature.location and right not in feature.location:
                if feature.location.start - left > 0:
                    dist = '-' + str(feature.location.start - left)
                else:
                    dist = '+' + str(abs(feature.location.start - left))
                closest_to_left_gene = [cds_features[0][0], dist, values]
                pos_gene_left = closest_to_left_gene
                other_gene = get_other_gene(reference, right, "right", cds_features, trna_features, rrna_features)
                pos_gene_right = other_gene
                return pos_gene_left, pos_gene_right
            elif left not in feature.location and right in feature.location:
                if feature.location.start - right > 0:
                    dist = '-' + str(feature.location.start - right)
                else:
                    dist = '+' + str(abs(feature.location.start - right))
                closest_to_right_gene = [feature.qualifiers[cds_features[0]][0], dist, values]
                pos_gene_right = closest_to_right_gene
                other_gene = get_other_gene(reference, left, "left", cds_features, trna_features, rrna_features)
                pos_gene_left = other_gene
                return pos_gene_left, pos_gene_right
            else:
                #the positions aren't in the middle of gene, so now need to see how close we are
                #to the current feature
                if feature.location.start - left > 0:
                    dist = '-' + str(feature.location.start - left)
                else:
                    dist = '+' + str(abs(feature.location.start - left))
                distance_with_left[abs(feature.location.start - left)] = [feature.qualifiers[cds_features[0]][0], dist, values]
                if feature.location.start - right > 0:
                    dist = '-' + str(feature.location.start - right)
                else:
                    dist = '+' + str(abs(feature.location.start - right))
                distance_with_right[abs(feature.location.start - right)] = [feature.qualifiers[cds_features[0]][0], dist, values]
            
    #we never broke out of the function, so it must mean that the insertion site
    #is intergenic                                                          
    distance_lkeys = list(OrderedDict.fromkeys(distance_with_left))
    closest_to_left_gene = distance_with_left[min(distance_lkeys)]
    pos_gene_left = closest_to_left_gene
    distance_rkeys = list(OrderedDict.fromkeys(distance_with_right))
    closest_to_right_gene = distance_with_right[min(distance_rkeys)]
    pos_gene_right = closest_to_right_gene
    #we already know that the gene isn't interrupted
    if closest_to_left_gene[0] == closest_to_right_gene[0]:
        if closest_to_left_gene[1] > closest_to_right_gene[1]:
            direction = "left"
            other_gene = get_other_gene(reference, left, direction, cds_features, trna_features, rrna_features)
            return other_gene, pos_gene_right
        elif closest_to_right_gene > closest_to_left_gene:
            direction = "right"
            other_gene = get_other_gene(reference, right, direction, cds_features, trna_features, rrna_features)
            return pos_gene_left, other_gene
    return pos_gene_left, pos_gene_right

def get_other_gene(reference, pos, direction, cds_features, trna_features, rrna_features):
    gb = SeqIO.read(reference, "genbank")
    distance = {}
    for feature in gb.features:
        #only want to look for genes that are to the left of the gene that has
        #already been found
        if feature.type == "CDS" or feature.type == "tRNA" or feature.type == "rRNA":
            values = get_qualifiers(cds_features, trna_features, rrna_features, feature)
            values.append(feature.strand)
            if direction == "left":
                #for this to be true, the position we're looking at must be
                #larger than the gene start and end (if the position is not
                #in the gene)
                if pos >= feature.location.start and feature.location.end:
                    #always want to refer to the start codon
                    if feature.location.start - pos > 0:
                        dist = '-' + str(feature.location.start - pos)
                    else:
                        dist = '+' + str(abs(feature.location.start - pos))

                    distance[abs(feature.location.start - pos)] = [feature.qualifiers[cds_features[0]][0], dist, values]
            elif direction == "right":
                if pos <= feature.location.start and feature.location.end:
                    if feature.location.start - pos > 0:
                        dist = '-' + str(feature.location.start - pos)
                    else:
                        dist = '+' + str(abs(feature.location.start - pos))
                    distance[abs(feature.location.start - pos)] = [feature.qualifiers[cds_features[0]][0], dist, values]
                    
    distance_keys = list(OrderedDict.fromkeys(distance))
    closest_gene = distance[min(distance_keys)]
    return closest_gene

def blast_db(fasta):
    
    if not os.path.exists(fasta + '.nin'):
        os.system('makeblastdb -in ' + fasta + ' -dbtype nucl')

def main():

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))
    list_of_isolates = []

    list_of_positions = collections.defaultdict(dict) # key1 = pos, key2 = isolate, value = +/-
    position_orientation = {}

    reference_fasta = args.reference_gbk.split('.g')[0]
    gbk_to_fasta(args.genbank, reference_fasta)

    blast_db(reference_fasta)

    list_of_positions, position_orientation, ref_name = get_ref_positions(reference_fasta, args.seq, list_of_positions, position_orientation)

    for result_file in unique_results_files:
        isolate = result_file.split('__')[0]
        list_of_isolates.append(isolate)
        header = 0
        with open(result_file) as file_open:
            for line in file_open:
                if header == 0:
                    header += 1
                elif 'No hits found' not in line and line != '':
                    #print isolate
                    info = line.strip('\n').split('\t')
                    #print info
                    orientation = info[1]
                    is_start = int(info[2])
                    #print is_start
                    is_end = int(info[3])
                    #print is_end
                    if (is_start, is_end) not in list_of_positions:
                        if list_of_positions.keys() != []:
                            old_range, new_range, new_orientation = check_ranges(position_orientation, (is_start, is_end), args.gap, orientation)
                            #print old_range, new_range, new_orientation
                            if old_range != False:
                                store_values = list_of_positions[old_range]
                                del list_of_positions[old_range]
                                list_of_positions[new_range] = store_values
                                list_of_positions[new_range][isolate] = '+'
                                del position_orientation[old_range]
                                position_orientation[new_range] = new_orientation
                            else:
                                list_of_positions[(is_start, is_end)][isolate] = '+'
                                position_orientation[(is_start, is_end)] = orientation
                        else:
                            list_of_positions[(is_start, is_end)][isolate] = '+'
                            position_orientation[(is_start, is_end)] = orientation
                    elif (is_start, is_end) in list_of_positions:
                        list_of_positions[(is_start, is_end)][isolate] = '+'
                    #print position_orientation

    # ordering positions from smallest to largest for final table output
    order_position_list = list(OrderedDict.fromkeys(list_of_positions.keys()))
    order_position_list.sort()
    #print order_position_list

    # create header of table
    with open(args.output, 'w') as out:
        header = ['isolate']
        for position in order_position_list:
            header.append(str(position[0]) + '-' + str(position[1]))
        #print '\t'.join(header)
        out.write('\t'.join(header) + '\n')

        row = [ref_name]
        for position in order_position_list:
            if ref_name in list_of_positions[position]:
                row.append(list_of_positions[position][ref_name])
            else:
                row.append('-')
        #print '\t'.join(row)
        out.write('\t'.join(row) + '\n')
        
        # create each row
        for isolate in list_of_isolates:
            row = [isolate]
            for position in order_position_list:
                if isolate in list_of_positions[position]:
                    row.append(list_of_positions[position][isolate])
                else:
                    row.append('-')
            #row.append('\n')
            #print '\t'.join(row)
            out.write('\t'.join(row) + '\n')

        row_l_locus = ['left locus tag']
        row_r_locus = ['right locus tag']
        row_l_dist = ['left distance']
        row_r_dist = ['right distance']
        row_l_prod = ['left product']
        row_r_prod = ['right product']

        for position in order_position_list:
            genes_before, genes_after = get_flanking_genes(args.reference_gbk, position[0], position[1], args.cds, args.trna, args.rrna)
            print genes_before
            print genes_after
            row_l_locus.append(genes_before[0])
            row_r_locus.append(genes_after[0])
            if genes_before[0] == genes_after[0]:
                row_l_dist.append(genes_before[1])
                row_r_dist.append(genes_before[1])
            else:
                row_l_dist.append(genes_before[1])
                row_r_dist.append(genes_after[1])
            row_l_prod.append(genes_before[2][:-1])
            row_r_prod.append(genes_after[2][:-1])
        out.write('\t'.join(row_l_locus) + '\n')
        out.write('\t'.join(row_l_dist) + '\n')
        out.write('\t.'.join(str(i) for i in row_l_prod) + '\n')
        out.write('\t'.join(row_r_locus) + '\n')
        out.write('\t'.join(row_r_dist) + '\n')
        out.write('\t'.join(str(i) for i in row_r_prod) + '\n')

if __name__ == "__main__":
    main()
