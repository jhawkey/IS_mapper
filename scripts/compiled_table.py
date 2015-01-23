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
from ismap import gbk_to_fasta
import time

def parse_args():

    parser = ArgumentParser(description="Create a table of IS hits in all isolates for ISMapper")
    # Inputs
    parser.add_argument('--tables', nargs='+', type=str, required=True, help='tables to compile')
    parser.add_argument('--reference_gbk', type=str, required=True, help='gbk file of reference to report closest genes')
    parser.add_argument('--seq', type=str, required=True, help='fasta file for insertion sequence looking for in reference')
    # Parameters for hits
    parser.add_argument('--gap', type=int, required=False, default=0, help='distance between regions to call overlapping')
    parser.add_argument('--cds', nargs='+', type=str, required=False, default=['locus_tag', 'gene', 'product'], help='qualifiers to look for in reference genbank for CDS features')
    parser.add_argument('--trna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for tRNA features')
    parser.add_argument('--rrna', nargs='+', type=str, required=False, default=['locus_tag', 'product'], help='qualifiers to look for in reference genbank for rRNA features')
    # Output parameters
    parser.add_argument('--output', type=str, required=True, help='name of output file')

    return parser.parse_args()

def check_ranges(ranges, range_to_check, gap, orientation):
    '''
    Takes a list of tuples with currently known ranges, and a new range
    to check against these to see if it overlaps.
    Also takes gap variable, indicating that the ranges may be gap distance
    apart and still merged.
    Also takes orientation, as only want to merge ranges that are the same 
    orientation.

    Returns the old range (to be replaced) and the new, merged range (to replace
        the old range with) and the orientation.
    If the range_to_check can't be merged with any of the known ranges, then
    just return False, False, False.
    '''
    # Get start and end coordinates
    start = range_to_check[0]
    stop = range_to_check[1]
    # Get current list of ranges
    range_list = ranges.keys()

    # For every range, check the orientation
    for i in range(0, len(range_list)):
        # Only merge hits that have the same orientation
        if orientation == ranges[(range_list[i][0], range_list[i][1])]:
            # Get coordinates of the current range to check
            x = range_list[i][0]
            y = range_list[i][1]
            # Forward orientations have certain rules
            if orientation == 'F':
                # The x value must lie between start and stop in test range 
                # taking into account gap
                if x in range(start - gap, stop + 1) or x in range(start, stop + gap + 1):
                    # If so, then these ranges overlap
                    new_start = min(x, start)
                    new_end = max(y, stop)
                    return range_list[i], (new_start, new_end), orientation
                # Otherwise the y value must lie between the start and stop in the test range
                # taking into account the gap
                elif y in range(start - gap, stop + 1) or y in range(start, stop + gap + 1):
                    # If so, then these ranges overlap
                    new_start = min(x, start)
                    new_end = min(y, stop)
                    return range_list[i], (new_start, new_end), orientation
            # Same goes for ranges that are in reverse orientation
            elif orientation == 'R':
                if x in range(start - gap, stop + 1) or x in range(start, stop + gap + 1):
                    new_start = min(x, start)
                    new_end = max(y, stop)
                    return range_list[i], (new_start, new_end), orientation
                elif y in range(start - gap, stop + 1) or y in range(start, stop + gap + 1):
                    new_start = min(x, start)
                    new_end = min(y, stop)
                    return range_list[i], (new_start, new_end), orientation
    # If the range doesn't overlap any currently know range, then return False
    return False, False, False

def get_ref_positions(reference, is_query, positions_dict, orientation_dict):
    '''
    Get the coordinates of known IS sites in the reference.

    Takes the reference genbank, the IS query and the dictionary to add_argument
    IS query positions into, as well as a dictionary to add orientations
    of each of these positions.
    Returns these positions and orientations, as well as the reference name
    for file naming.
    '''
    # Get the name of the IS query to create temp file
    is_name = os.path.split(is_query)[1]
    ref_name = os.path.split(reference)[1]
    blast_output = os.getcwd() + '/' + is_name + '_' + ref_name + '.tmp'

    # Create a BLAST database of the reference if there isn't one already
    # Should probably use this function - blast_db(reference)
    if not os.path.exists(reference):
        os.system('makeblastdb -in ' + reference + ' -dbtype nucl')
    # Do the BLAST
    blastn_cline = NcbiblastnCommandline(query=is_query, db=reference, outfmt="'6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs'", out=blast_output)
    stdout, stderr = blastn_cline()
    # Open the BLAST output and get IS query sites
    with open(blast_output) as out:
        for line in out:
            info = line.strip('\n').split('\t')
            # To be a known site, hast to match query at least 90% with coverage of 95
            if float(info[3]) >= 90 and float(info[4])/float(info[1]) * 100 >= 95:
                positions_dict[(int(info[6]), int(info[7]))][ref_name] = '+'
                # Get orientation of known site for merging purposes
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
    skips and keeps going.
    Returns a list of qualfiers found in that feature.
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

def get_main_gene_id(qualifier_list, feature):
    '''
    Takes a list of qualifiers and a genbank feature.
    Returns the name of the qualifier that contains
    the gene id.
    '''

    for qual in qualifier_list:
        try:
            id_name = feature.qualifiers[qual][0]
            return id_name
        except KeyError:
            pass

def get_flanking_genes(reference, left, right, cds_features, trna_features, rrna_features):
    '''
    Takes the reference genbank, left and right coordinates, as well as 
    qualifiers for the different possible features in the genbank.
    Looks at each feature in the genbank and measures how close it is to 
    the left and right coordinates.

    Returns the genes closest to the left (closest_left) and right (closest_right)
    coordinates.
    '''

    gb = SeqIO.read(reference, 'genbank')
    pos_gene_left = []
    pos_gene_right = []
    distance_with_left = {}
    distance_with_right = {}
    #print left
    #print right

    for feature in gb.features:
        if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
            values = get_qualifiers(cds_features, trna_features, rrna_features, feature)
            values.append(feature.strand)
            if feature.strand == 1:
                pos = feature.location.start
            else:
                pos = feature.location.end
            # First check to see if both coordinates fit into the feature
            if left in feature.location and right in feature.location:
                # We want the absolute value because a value with no sign in the final table
                # indicates that the gene is interrupted
                gene_id = get_main_gene_id(cds_features, feature)
                gene = [gene_id, str(abs(pos - left)), values]
                pos_gene_left = gene
                pos_gene_right = [gene_id, str(abs(pos - right)), values]
                return pos_gene_left, pos_gene_right
            elif left in feature.location and right not in feature.location:
                gene_id = get_main_gene_id(cds_features, feature)
                if pos - left > 0:
                    dist = '-' + str(pos - left)
                else:
                    dist = '+' + str(abs(pos - left))
                closest_to_left_gene = [gene_id, dist, values]
                pos_gene_left = closest_to_left_gene
                other_gene = get_other_gene(reference, right, "right", cds_features, trna_features, rrna_features)
                pos_gene_right = other_gene
                return pos_gene_left, pos_gene_right
            elif left not in feature.location and right in feature.location:
                gene_id = get_main_gene_id(cds_features, feature)
                if pos - right > 0:
                    dist = '-' + str(pos - right)
                else:
                    dist = '+' + str(abs(pos - right))
                closest_to_right_gene = [gene_id, dist, values]
                pos_gene_right = closest_to_right_gene
                other_gene = get_other_gene(reference, left, "left", cds_features, trna_features, rrna_features)
                pos_gene_left = other_gene
                return pos_gene_left, pos_gene_right
            else:
                #the positions aren't in the middle of gene, so now need to see how close we are
                #to the current feature
                gene_id = get_main_gene_id(cds_features, feature)
                if pos - left > 0:
                    dist = '-' + str(pos - left)
                else:
                    dist = '+' + str(abs(pos - left))
                distance_with_left[abs(pos - left)] = [gene_id, dist, values]
                if pos - right > 0:
                    dist = '-' + str(pos - right)
                else:
                    dist = '+' + str(abs(pos - right))
                distance_with_right[abs(pos - right)] = [gene_id, dist, values]
            
    #we never broke out of the function, so it must mean that the insertion site
    #is intergenic                                                          
    distance_lkeys = list(OrderedDict.fromkeys(distance_with_left))
    closest_to_left_gene = distance_with_left[min(distance_lkeys)]
    pos_gene_left = closest_to_left_gene
    distance_rkeys = list(OrderedDict.fromkeys(distance_with_right))
    closest_to_right_gene = distance_with_right[min(distance_rkeys)]
    pos_gene_right = closest_to_right_gene
    #print closest_to_left_gene
    #print closest_to_right_gene

    #we already know that the gene isn't interrupted
    if closest_to_left_gene[0] == closest_to_right_gene[0]:
        #print 'same gene'
        if closest_to_left_gene[1] > closest_to_right_gene[1]:
            #print 'we look left'
            direction = "left"
            other_gene = get_other_gene(reference, left, direction, cds_features, trna_features, rrna_features)
            return other_gene, pos_gene_right
        elif closest_to_left_gene[1] < closest_to_right_gene[1]:
            #print 'we look right'
            direction = "right"
            other_gene = get_other_gene(reference, right, direction, cds_features, trna_features, rrna_features)
            return pos_gene_left, other_gene
    return pos_gene_left, pos_gene_right

def get_other_gene(reference, pos, direction, cds_features, trna_features, rrna_features, known=False):
    '''
    Takes reference genbank to look for features in, a coordinate (pos), 
    the direction to look in (upstream or downstream), and qualifiers for
    each of the feature types. If known is set to True, 
    this tells the function that the position could be inside a gene that 
    is flanking a known IS site in the reference.

    Measure distance compared to the start codon of the feature.
    (could be feature.location.start or .end depending on strand)
    - value before the distance indicates we are downstream
    + value before the distance indicates we are upstream

    Returns the gene closest to the coordinate given (pos).
    '''

    gb = SeqIO.read(reference, "genbank")
    distance = {}

    for feature in gb.features:
        # Only want to look for genes that are to the left of the gene that has
        # already been found
        if feature.type == "CDS" or feature.type == "tRNA" or feature.type == "rRNA":
            values = get_qualifiers(cds_features, trna_features, rrna_features, feature)
            values.append(feature.strand)
            if feature.strand == 1:
                # Foward strand, so start and end are simple
                feature_start = feature.location.start
                feature_end = feature.location.end
            else:
                # Reverse strand, so start of gene is actually
                # the end of the location for the feature
                feature_start = feature.location.end
                feature_end = feature.location.start
            if direction == "left":
                # For this to be true, the position we're looking at must be
                # larger than the gene start and end (if the position is not
                # in the gene)
                if pos in range(min(feature_start, feature_end), max(feature_start, feature_end)) and known == True:
                    # We're inside a gene, and this is a known hit, so the 
                    # flanking region could be inside the gene, but the gene
                    # is not necessarily interrupted by the known site.
                    gene_id = get_main_gene_id(cds_features, feature)
                    if feature_start - pos > 0:
                        dist = '-' + str(feature_start - pos)
                    else:
                        dist = '+' + str(abs(feature_start - pos))
                    closest_gene = [gene_id, dist, values]
                    # Return this gene as it must be the answer if we're inside
                    # this feature
                    return closest_gene
                elif pos > feature_start and feature_end:
                    # Otherwise just check to see how close the
                    # pos is to this gene
                    gene_id = get_main_gene_id(cds_features, feature)
                    # Always want to refer to the start codon
                    if feature_start - pos > 0:
                        dist = '-' + str(feature_start - pos)
                    else:
                        dist = '+' + str(abs(feature_start - pos))
                    distance[abs(feature_start - pos)] = [gene_id, dist, values]
            elif direction == "right":
                if pos in range(min(feature_start, feature_end), max(feature_start, feature_end)) and known == True:
                    # We're inside a gene, and this is a known hit, so the 
                    # flanking region could be inside the gene, but the gene
                    # is not necessarily interrupted by the known site.
                    gene_id = get_main_gene_id(cds_features, feature)
                    if feature_start - pos > 0:
                        dist = '-' + str(feature_start - pos)
                    else:
                        dist = '+' + str(abs(feature_start - pos))
                    closest_gene = [gene_id, dist, values]
                    # Return this gene as it must be the answer if we're inside
                    # this feature
                    return closest_gene
                elif pos < feature_start and feature_end:
                    # Otherwise just check to see how close the
                    # pos is to this gene
                    gene_id = get_main_gene_id(cds_features, feature)
                    # Always want to refer to the start codon
                    if feature_start - pos > 0:
                        dist = '-' + str(feature_start - pos)
                    else:
                        dist = '+' + str(abs(feature_start - pos))
                    distance[abs(feature_start - pos)] = [gene_id, dist, values]
                    
    # Get all the distances and order them
    distance_keys = list(OrderedDict.fromkeys(distance))
    # If an empty list is returned, then we must be either at the very
    # beginning or the very end of the genbank
    if len(distance_keys) == 0:
        # If the distance to the end of the genbank is
        # smaller than the distance from the start of the
        # genbank, we're at the end
        if abs(pos - len(gb)) < abs(1 - pos):
            # We need the first gene
            feature_no = 0
            for feature in gb.features:
                if feature_no == 0 and (feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA'):
                    gene_id = get_main_gene_id(cds_features, feature)
                    values = get_qualifiers(cds_features, trna_features, rrna_features, feature)
                    values.append(feature.strand)
                    closest_gene = [gene_id, 'start of genbank', values]
                    feature_no += 1
                    return closest_gene
                else:
                    pass
        # Otherwise we're closest to the start
        else:
            # We need the last gene
            closest_gene = []
            feature_no = 1
            while closest_gene == []:
                index_no = '-' + str(feature_no)
                if gb.features[int(index_no)].type == 'CDS' or gb.features[int(index_no)].type == 'tRNA' or gb.features[int(index_no)].type == 'rRNA':
                    gene_id = get_main_gene_id(cds_features, feature)
                    values = get_qualifiers(cds_features, trna_features, rrna_features, feature)
                    values.append(feature.strand)
                    closest_gene = [gene_id, 'end of genbank', values]
                    return closest_gene
                else:
                    feature_no += 1

    # The closest gene is the one with the smallest distance
    closest_gene = distance[min(distance_keys)]
    return closest_gene

def blast_db(fasta):
    '''
    Takes a fasta file and creates a BLAST database 
    if one doesn't exist already.
    '''
    
    if not os.path.exists(fasta + '.nin'):
        os.system('makeblastdb -in ' + fasta + ' -dbtype nucl')

def main():
    
    start_time = time.time()

    args = parse_args()

    unique_results_files = list(OrderedDict.fromkeys(args.tables))
    list_of_isolates = []

    # key1 = (start, end), key2 = isolate, value = +/*/?
    list_of_positions = collections.defaultdict(dict)
    # key1 = (start, end), key2 = ref, value = +
    list_of_ref_positions = collections.defaultdict(dict)
    # key = (start, end), value = orientation (F/R)
    position_orientation = {}

    reference_fasta = args.reference_gbk.split('.g')[0]
    # Create a fasta file of the reference for BLAST
    print 'Creating fasta file and database of reference ...'
    gbk_to_fasta(args.reference_gbk, reference_fasta)
    # Make a BLAST database
    blast_db(reference_fasta)
    # Get the reference positions and orientations for this IS query
    print '\nGetting query positions in reference ...'
    list_of_ref_positions, ref_position_orientation, ref_name = get_ref_positions(reference_fasta, args.seq, list_of_ref_positions, position_orientation)

    elapsed_time = time.time() - start_time
    print 'Time taken: ' + str(elapsed_time)
    # Loop through each table give to --tables
    print 'Collating results files ...'
    for result_file in unique_results_files:
        # Get isolate name
        isolate = result_file.split('_table.txt')[0]
        list_of_isolates.append(isolate)
        # Skip the header
        header = 0
        with open(result_file) as file_open:
            for line in file_open:
                # Skip header
                if header == 0:
                    header += 1
                # Check to make sure there were actually hits
                elif 'No hits found' not in line and line != '':
                    info = line.strip('\n').split('\t')
                    # Get orientation for hit and start/end coordinates
                    orientation = info[1]
                    is_start = int(info[2])
                    is_end = int(info[3])
                    # Note whether call is Known, Novel or Possible related IS
                    call = info[5]
                    # If the position Know and therefore in the reference, 
                    # compare the hit to the known reference positions
                    if call == 'Known' or call == 'Known?':
                        if (is_start, is_end) not in list_of_ref_positions:
                            # Looking for hits that are with 100 bp of the Known reference hits
                            old_range, new_range, new_orientation = check_ranges(ref_position_orientation, (is_start, is_end), 100, orientation)
                            # If we can merge ranges
                            if old_range != False:
                                store_values = list_of_ref_positions[old_range]
                                # Remove the old range, and add the new range
                                del list_of_ref_positions[old_range]
                                list_of_ref_positions[new_range] = store_values
                                # Note whether the hit is uncertain (?)
                                # or confident (+)
                                if '?' in call:
                                    list_of_ref_positions[new_range][isolate] = '?'
                                else:
                                    list_of_ref_positions[new_range][isolate] = '+'
                                # Remove the old range from the reference positions
                                # and add the new merged range
                                del ref_position_orientation[old_range]
                                ref_position_orientation[new_range] = new_orientation
                            # If we can't merge with a known position
                            else:
                                # Mark as uncertain if ? in call
                                if '?' in call:
                                    list_of_positions[(is_start, is_end)][isolate] = '?'
                                # Otherwise just append it as a new reference position
                                else:
                                    list_of_ref_positions[(is_start, is_end)][isolate] = '+'
                                # Note the orientation of the hit
                                position_orientation[(is_start, is_end)] = orientation
                    # Otherwise try and merge with positions that are novel
                    elif (is_start, is_end) not in list_of_positions:
                        # If the list of positions isn't empty, then there are ranges to check against
                        if list_of_positions.keys() != []:
                            old_range, new_range, new_orientation = check_ranges(position_orientation, (is_start, is_end), args.gap, orientation)
                            # So the current range overlaps with a range we already have
                            if old_range != False:
                                # Remove the old range and add the new one
                                store_values = list_of_positions[old_range]
                                del list_of_positions[old_range]
                                list_of_positions[new_range] = store_values
                                # Mark as ? if uncertain, * if imprecise
                                # or + if confident
                                if '?' in call:
                                    list_of_positions[new_range][isolate] = '?'
                                elif '*' in call:
                                    list_of_positions[new_range][isolate] = '*'
                                else:
                                    list_of_positions[new_range][isolate] = '+'
                                # Remove the old range from the orientations 
                                # and add the new one
                                del position_orientation[old_range]
                                position_orientation[new_range] = new_orientation
                            # Otherwise this range hasn't been seen before, so all values are False
                            else:
                                if '?' in call:
                                    list_of_positions[(is_start, is_end)][isolate] = '?'
                                elif '*' in call:
                                    list_of_positions[(is_start, is_end)][isolate] = '*'
                                else:
                                    list_of_positions[(is_start, is_end)][isolate] = '+'
                                position_orientation[(is_start, is_end)] = orientation
                        # Otherwise the position list is empty, so there are no ranges to check against
                        else:
                            if '?' in call:
                                    list_of_positions[(is_start, is_end)][isolate] = '?'
                            elif '*' in call:
                                list_of_positions[(is_start, is_end)][isolate] = '*'
                            else:
                                list_of_positions[(is_start, is_end)][isolate] = '+'
                            position_orientation[(is_start, is_end)] = orientation
                    # This position is already in the list, so just append the new isolate
                    elif (is_start, is_end) in list_of_positions:
                        if '?' in call:
                            list_of_positions[(is_start, is_end)][isolate] = '?'
                        elif '*' in call:
                            list_of_positions[(is_start, is_end)][isolate] = '*'
                        else:
                            list_of_positions[(is_start, is_end)][isolate] = '+'

    elapsed_time = time.time() - start_time
    print 'Time taken: ' + str(elapsed_time)

    # key = (start, end), valye = [left_gene, right_gene]
    position_genes = {}
    # Get the flanking genes for each know position
    print 'Getting flanking genes for each position (this step is the longest and could take some time) ...'
    if len(list_of_ref_positions.keys()) != 0:
        for position in list_of_ref_positions.keys():
            left_pos = min(position[0], position[1])
            right_pos = max(position[0], position[1])
            flanking_left = get_other_gene(args.reference_gbk, left_pos, "left", args.cds, args.trna, args.rrna)
            flanking_right = get_other_gene(args.reference_gbk, right_pos, "right", args.cds, args.trna, args.rrna)
            position_genes[(position[0], position[1])] = [flanking_left, flanking_right]
    # Get flanking genes for novel positions
    for position in list_of_positions.keys():
        genes_before, genes_after = get_flanking_genes(args.reference_gbk, position[0], position[1], args.cds, args.trna, args.rrna)
        position_genes[(position[0], position[1])] = [genes_before, genes_after]

    # Order positions from smallest to largest for final table output
    order_position_list = list(OrderedDict.fromkeys(list_of_positions.keys())) + list(OrderedDict.fromkeys(list_of_ref_positions.keys()))
    order_position_list.sort()

    elapsed_time = time.time() - start_time
    print 'Time taken: ' + str(elapsed_time)

    # Create header of table
    print 'Writing output table to ' + args.output + ' ...'
    with open(args.output, 'w') as out:
        header = ['isolate']
        for position in order_position_list:
            header.append(str(position[0]) + '-' + str(position[1]))
        out.write('\t'.join(header) + '\n')
        # Add the values for the reference positions
        row = [ref_name]
        for position in order_position_list:
            if position in list_of_ref_positions:
                if ref_name in list_of_ref_positions[position]: 
                    row.append(list_of_ref_positions[position][ref_name])
            else:
                row.append('-')
        out.write('\t'.join(row) + '\n')
        
        # Loop through each isoalte
        # and create each row
        for isolate in list_of_isolates:
            row = [isolate]
            for position in order_position_list:
                if position in list_of_ref_positions:
                    if isolate in (list_of_ref_positions[position]):
                        row.append(list_of_ref_positions[position][isolate])
                    else:
                        row.append('-')
                elif position in list_of_positions:
                    if isolate in list_of_positions[position]:
                        row.append(list_of_positions[position][isolate])
                    else:
                        row.append('-')
            out.write('\t'.join(row) + '\n')
        # Set up flanking genes
        row_l_locus = ['left locus tag']
        row_r_locus = ['right locus tag']
        row_l_dist = ['left distance']
        row_r_dist = ['right distance']
        row_l_prod = ['left product']
        row_r_prod = ['right product']

        # Print flanking genes for each position
        for position in order_position_list:
            #   genes_before, genes_after = get_flanking_genes(args.reference_gbk, position[0], position[1], args.cds, args.trna, args.rrna)
            if position in position_genes:
                row_l_locus.append(position_genes[position][0][0])
                row_r_locus.append(position_genes[position][1][0])
                row_l_dist.append(position_genes[position][0][1])
                row_r_dist.append(position_genes[position][1][1])
                row_l_prod.append(position_genes[position][0][2][:-1])
                row_r_prod.append(position_genes[position][1][2][:-1])
        out.write('\t'.join(row_l_locus) + '\n')
        out.write('\t'.join(row_l_dist) + '\n')
        out.write('\t.'.join(str(i) for i in row_l_prod) + '\n')
        out.write('\t'.join(row_r_locus) + '\n')
        out.write('\t'.join(row_r_dist) + '\n')
        out.write('\t'.join(str(i) for i in row_r_prod) + '\n')

    elapsed_time = time.time() - start_time
    print 'Table compilation finished in ' + str(elapsed_time) 

if __name__ == "__main__":
    main()
