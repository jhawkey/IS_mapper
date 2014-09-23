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
    parser.add_argument('--reference_fasta', type=str, required=True, help='fasta file of reference to determine known positions')
    parser.add_argument('--reference_gbk', type=str, required=True, help='gbk file of reference to report closest genes')
    parser.add_argument('--seq', type=str, required=True, help='fasta file for insertion sequence looking for in reference')
    parser.add_argument('--gap', type=int, required=False, default=400, help='distance between regions to call overlapping')
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

def get_flanking_genes(reference, positions):

    gb = SeqIO.read(reference, 'genbank')
    pos_gene_start = {}
    pos_gene_end = {}
    for pos in positions:
        #print pos
        x = pos[0]
        y = pos[1]
        distance_start = {}
        distance_end = {}
        for feature in gb.features:
            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                #print feature
                if feature.type == 'CDS':
                    distance_start[abs(feature.location.start - x)] = feature.qualifiers['locus_tag'][0]
                elif feature.type == 'tRNA' or feature.type == 'rRNA':
                    distance_start[abs(feature.location.start - x)] = feature.qualifiers['product'][0]
        distance_skeys = list(OrderedDict.fromkeys(distance_start))
        gene = distance_start[min(distance_skeys)]
        pos_gene_start[pos] = gene
        for feature in gb.features:
            if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
                if feature.type == 'CDS':
                    distance_end[abs(feature.location.end - y)] = feature.qualifiers['locus_tag'][0]
                elif feature.type == 'tRNA' or feature.type == 'rRNA':
                    distance_end[abs(feature.location.end - y)] = feature.qualifiers['product'][0]
        distance_ekeys = list(OrderedDict.fromkeys(distance_end))
        gene2 = distance_end[min(distance_ekeys)]
        pos_gene_end[pos] = gene2

    pos_check = {}
    '''for pos in pos_gene_start:
        if pos_gene_start[pos] == pos_gene_end[pos]:
            x = pos[0]
            y = pos[1]
            gene_test = pos_gene_start[pos]
            if x < y:
                for feature in gb.features:
                    if feature.qualifiers['locus_tag'][0] == gene_test and feature.strand == 1:
                        distance_x = abs(feature.location.start - x)
                        distance_y = abs(eature.location.start - y)
                        if distance_x < distance_y:
                            pos_check[(x,y)] = 'y+20'
                        else:
                            pos_check[(x,y)] = 'x-20'
                    elif feature.qualifiers['locus_tag'][0] == gene_test and feature.strand == -1:
                        distance_x = abs(feature.location.end - x)
                        distance_y = abs(feature.location.end - y)
                        if distance_x < distance_y:
                            pos_check[(x, y)] = 'y+20'
                        else:
                            pos_check[(x, y)] = 'x-20'
            else:
                for feature in gb.features:
                    if features.qualifiers['locus_tag'][0] == gene_test and feature.strand == 1:
                        distance_x = abs(feature.location.start - x)
                        distance_y = abs(eature.location.start - y)
                        if distance_y < distance_x:
                            pos_check[(x,y)] = 'x+20'
                        else:
                            pos_check[(x,y)] = 'y-20'
                    elif feature.qualifiers['locus_tag'][0] == gene_test and feature.strand == -1:
                        distance_x = abs(feature.location.end - x)
                        distance_y = abs(feature.location.end - y)
                        if distance_y < distance_x:
                            pos_check[(x, y)] = 'x+20'
                        else:
                            pos_check[(x, y)] = 'y-20'

    # gotta fix this!
    for position in pos_check:
        if pos_check[position] == 'y+20':
            pass'''

    return pos_gene_start, pos_gene_end

def blast_db(fasta):
    
    if not os.path.exists(fasta + '.nin'):
        os.system('makeblastdb -in ' + fasta + ' -dbtype nucl')

def main():

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))
    list_of_isolates = []

    list_of_positions = collections.defaultdict(dict) # key1 = pos, key2 = isolate, value = +/-
    position_orientation = {}

    blast_db(args.reference_fasta)

    list_of_positions, position_orientation, ref_name = get_ref_positions(args.reference_fasta, args.seq, list_of_positions, position_orientation)

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

        genes_before, genes_after = get_flanking_genes(args.reference_gbk, order_position_list)

        row = ['flanking genes']
        for position in order_position_list:
            row.append(genes_before[position])
        #print '\t'.join(row)
        out.write('\t'.join(row) + '\n')
        row = ['flanking genes']
        for position in order_position_list:
            row.append(genes_after[position])
        #print '\t'.join(row)
        out.write('\t'.join(row) + '\n')


if __name__ == "__main__":
    main()
