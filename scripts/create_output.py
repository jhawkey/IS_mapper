#!/usr/bin/env python

import os, logging
import operator
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from run_commands import run_command

class ISHit(object):
    def __init__(self, left_pos, right_pos):
        # hit type will either be novel or known
        self.hit_type = None
        # left and right features are the actual feature surrounding the hit
        # 0: locus tag, 1: distance from start site, 2: all other qualifiers
        self.left_feature = None
        self.right_feature = None
        # These will be integers of either 1 or -1 to tell you the strand of the
        # genes surrounding the hit
        self.left_strand = None
        self.right_strand = None
        # These will be integers detailing the distance from the IS hit to
        # either the left or right genes
        self.left_distance = None
        self.right_distance = None

        # these are the gene IDs (locus_tags) closest to the IS hit
        self.gene_left = None
        self.gene_right = None
        # interrupted only becomes true if gene left and right are the same
        self.interrupted = False
        # confidence level is either confident, imprecise (*) or unpaired (?)
        self.confidence_level = None
        # x and y coordinates of the hit
        # x is always the smaller and y the larger, REGARDLESS of orientation!
        self.x = left_pos
        self.y = right_pos
        # set this to true when the hit has come from the intersect file
        self.overlap = False
        # will be either forward or reverse
        self.orientation = None
        # the distance between the coordinates
        self.gap = None

        # store the blast values if its a known hit
        self.per_id = ''
        self.coverage = ''

    def get_gap_distance(self):

        self.gap = self.y - self.x

        if self.overlap:
            self.gap = -self.gap

        return self

    def determine_interrupted(self):

        if self.gene_left == self.gene_right:
            self.interrupted = True

        return self



    def get_flanking_genes(self, genbank_obj, feature_list):
        '''
        :param genbank_obj: The parsed genbank object
        :param feature_list: The parsed features list from the genbank file for easy searching
        :param is_hit: The IS hit object
        :return: The modified IS hit object containing the features of the left and right genes, their strand, and
        the distance from the IS hit to each gene
        '''

        # Find the correct indexes
        left_feature_index = self.binary_search(feature_list, 'L')
        right_feature_index = self.binary_search(feature_list, 'R')
        # Print out information if returning one of the error codes
        if type(left_feature_index) != int or type(right_feature_index) != int:
            print('left index')
            print(left_feature_index)
            print('right index')
            print(right_feature_index)
            print('left position: ' + str(self.left_pos))
            print('right position: ' + str(self.right_pos))
        # Extract the SeqFeature object that corresponds to that index
        left_feature = genbank_obj.features[left_feature_index]
        right_feature = genbank_obj.features[right_feature_index]

        # Add the SeqFeatures to the IS hit object, and strand information
        self.left_feature = left_feature
        self.right_feature = right_feature
        self.left_strand = left_feature.strand
        self.right_strand = right_feature.strand

        # Get the distances from the IS hit to the left and right genes

        # The distance to the left gene is the endmost position of the feature - the left IS coord
        left_dist = abs(max(left_feature.location.start, left_feature.location.end) - self.x)
        # The distance to the right gene is the startmost position of the feature - the right IS coord
        right_dist = abs(min(right_feature.location.start, right_feature.location.end) - self.y)

        # Here is probably a good place to check if we've got a position that wraps around from start
        # to end of the reference
        # Get the size of the reference genome
        genome_size = len(genbank_obj.seq)
        # If we've got a distance that is close to the size of the reference, then we know we need to
        # alter it
        if left_dist in range(int(round(genome_size * 0.9)), int(round(genome_size * 1.1))):
            # The the left hand feature is at the end of the genome
            # Distance from IS position to start of the genome is the
            # position itself
            dist_to_start = self.x
            # Distance from the end of the final gene to the end of the
            # genome
            dist_to_end = abs(left_feature.location.end - genome_size)
            # So the total distance is those two added together
            left_dist = dist_to_start + dist_to_end

        elif right_dist in range(int(round(genome_size * 0.9)), int(round(genome_size * 1.1))):
            # Then the right hand feature is at the start of the genome
            # Distance from the IS position to the end of the genome
            dist_to_end = abs(genome_size - self.y)
            # Distance from the start of the genome to the start of the first feature
            # is the start position of the first feature
            dist_to_feature = right_feature.location.start
            # So the total distance is those two added together
            right_dist = dist_to_end + dist_to_feature

        self.left_distance = left_dist
        self.right_distance = right_dist

        # add the gene names
        self.gene_left = self.left_feature.qualifiers['locus_tag'][0]
        self.gene_right = self.right_feature.qualifiers['locus_tag'][0]

        return self

    def binary_search(self, features, direction):
        if direction == 'R':
            isPosition = self.y
        else:
            isPosition = self.x
        min = 0
        max = len(features) - 1

        while True:

            # If the min has exceeded the max, then the IS position is not
            # inside a feature, and m will now be pointing to a
            # feature next to the IS position.
            if max < min:
                if direction == 'R':
                    return self.findFeatureAfterPosition(features, isPosition, m)
                else:
                    return self.findFeatureBeforePosition(features, isPosition, m)

            # Find the midpoint and save the feature attributes
            m = (min + max) // 2
            featureStart = features[m][0]
            featureEnd = features[m][1]
            featureIndex = features[m][2]

            # If the IS position is after the feature, move the minimum to
            # be after the feature.
            if featureEnd < isPosition:
                min = m + 1

            # If the IS position is before the feature, move the maximum to
            # be before the feature.
            elif featureStart > isPosition:
                max = m - 1

            # If the IS position is inside the feature, return only that feature
            elif isPosition >= featureStart and isPosition <= featureEnd:
                return featureIndex

            else:
                return "1 - THIS SHOULDN'T HAPPEN!"

    def findFeatureBeforePosition(self, features, isPosition, m):
        # If we are looking for the feature to the left of the
        # IS position, then either m-1 or m is our answer

        # If the start of the m feature is after the IS position,
        # then m is after the IS and m-1 is the correct feature
        if features[m][0] > isPosition:
            return features[m - 1][2]

        # If both m and m+1 features are before the IS position,
        # then m will be closer to the IS and is the correct feature
        elif features[m - 1][1] < isPosition and features[m][1] < isPosition:
            return features[m][2]

        else:
            return "2 - THIS SHOULDN'T HAPPEN!"

    def findFeatureAfterPosition(self, features, isPosition, m):
        # If we are looking for the feature to the right of the
        # IS position, then either m or m+1 is our answer

        # an index error will occur if m is the final feature, so just check that the first part is true
        # and return m
        try:
            features[m + 1]
        except IndexError:
            if features[m][0] > isPosition:
                return features[m][2]
            # otherwise we must be after the final position, so need to
            # return the start position of the very first feature
            else:
                return features[0][2]
        # If the end of the m feature is before the IS position,
        # then m is before the IS and m+1 is the correct feature
        if features[m][1] < isPosition:
            index = m + 1
            if index >= len(features):
                return features[0][2]
            return features[m + 1][2]

        # If both m and m+1 features are after the IS position,
        # then m will be closer to the IS and is the correct feature
        elif features[m][0] > isPosition and features[m + 1][0] > isPosition:
            return features[m][2]
        else:
            return "3 - THIS SHOULDN'T HAPPEN!"

def get_features(genbank_object):

    feature_list = []
    feature_count_list = 0
    feature_types = ["CDS", "tRNA", "rRNA"]

    for feature in genbank_object.features:
        if feature.type in feature_types:
            feature_list.append([int(feature.location.start), int(feature.location.end), feature_count_list])
            feature_count_list += 1
        else:
            feature_count_list += 1

    feature_list = sorted(feature_list, key=operator.itemgetter(0))
    return(feature_list)

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

def get_orientation(left_coords, right_coords):
    '''
    :param left_coords: list of coordinates for left end of hit
    :param right_coords: list of coordinates for right end of hit
    :return: return ISHit object, intialised with orienation and left/right positions
    '''

    # x must always be the smallest position, and y the largest position
    # regardless of orientation
    if left_coords[0] < right_coords[0] or left_coords[1] < right_coords[1]:
        smallest = min(right_coords[0], left_coords[1])
        biggest = max(right_coords[0], left_coords[1])
        new_hit = ISHit(smallest, biggest)
        # we are in forward orientation
        new_hit.orientation = 'F'
    else:
        smallest = min(left_coords[0], right_coords[1])
        biggest = max(left_coords[0], right_coords[1])
        new_hit = ISHit(smallest, biggest)
        # we are in reverse orientation
        new_hit.orientation = 'R'

    return new_hit

def doBlast(blast_input, blast_output, database):
    '''
    Perform a BLAST using the NCBI command line tools
    in BioPython.
    '''
    run_command(['makeblastdb', '-dbtype nucl', '-in', database], shell=True)
    #print(' '.join(['blastn', '-query', blast_input, '-db', database, '-outfmt "6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs"', '>', blast_output]))
    run_command(['blastn', '-query', blast_input, '-db', database, '-outfmt "6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs"', '>', blast_output], shell=True)

def check_seq_between(genbank_seq, insertion, start, end, name, temp):
    '''
    Check the sequence between two ends to see
    if it matches the IS query or not, and what
    the coverage and %ID to the query.

    :param genbank_seq: Whole sequence from genbank file
    :param insertion: IS query object to BLAST against
    :param start: Smallest coordinate, to extract sequence
    :param end: Largest coordinate, to extract sequence
    :param name: prefix for the file of this sequence
    :param temp: folder for the file of this sequence to go to
    :return: If there is a BLAST hit, return a dictionary with the 'coverage' and 'per_id' values, else return
            an empty dict
    '''

    # Get sequence between left and right ends
    seq_between = genbank_seq[start:end]
    # Turn the sequence into a fasta file
    seq_between = SeqRecord(Seq(str(seq_between), generic_dna), id=name)
    out_seq_between = os.path.join(temp, name + '.fasta')
    out_insertion = os.path.join(temp, name + 'ISseq.fasta')
    SeqIO.write(seq_between, out_seq_between, 'fasta')
    SeqIO.write(insertion, out_insertion, 'fasta')
    # Perform the BLAST
    doBlast(temp + name + '.fasta', temp + name + '_out.txt', temp + name + 'ISseq.fasta')
    # Only want the top hit, so set count variable to 0
    first_result = 0
    # Open the BLAST output file
    with open(temp + name + '_out.txt') as summary:
        for line in summary:
            # Get coverage and % ID for top hit
            if first_result == 0:
                info = line.strip().split('\t')
                hit = {'coverage': float(info[-1]), 'per_id': float(info[3])}
                first_result += 1
            return hit
    # If there is no hit, just return an empty dict
    return {}

def check_unpaired_hits(line_check, ref_gbk_obj, ref_feature_list, is_query_obj, min_range, max_range,
                        tmp_output_folder):

    # intialise a list of all the hits found in this file
    IS_hit_list = []
    # get length of IS
    is_query_length = len(is_query_obj.seq)

    # go through each line
    for info in line_check:
        # get the distance between the hits
        gap = int(info[6])
        # separate out info on the left and right sides of the hit
        intersect_left = [int(info[1]), int(info[2])]
        intersect_right = [int(info[4]), int(info[5])]
        # TODO: check_hit_within_hit
        # get the orientation and the IS hit object
        new_hit = get_orientation(intersect_left, intersect_right)
        # if the gap is small, it's a novel hit
        if gap <= 15:
            new_hit.hit_type = 'novel'
            new_hit.confidence_level = 'unpaired'
            new_hit.get_flanking_genes(ref_gbk_obj, ref_feature_list)
            IS_hit_list.append(new_hit)
        # if the gap is big enough, could be the IS itself, so do a BLAST check
        elif float(gap) / is_query_length >= min_range and float(gap) / is_query_length <= max_range:
            new_hit = get_orientation(intersect_left, intersect_right)
            seq_check_results = check_seq_between(ref_gbk_obj.seq, is_query_obj, new_hit.x,
                                                  new_hit.y, 'tmp_seq', tmp_output_folder)
            # if it's a good hit, add it
            if len(seq_check_results) != 0 and seq_check_results['per_id'] >= 80 and seq_check_results[
                'coverage'] >= 80:
                # get the flanking genes
                new_hit.get_flanking_genes(ref_gbk_obj, ref_feature_list)
                # make sure its a confident, novel hit
                new_hit.hit_type = 'known'
                new_hit.confidence_level = 'unpaired'
                new_hit.per_id = str(seq_check_results['per_id'])
                new_hit.coverage = str(seq_check_results['coverage'])
                # add it to the list
                IS_hit_list.append(new_hit)
            # if the thresholds are low, then mark it as a possible related IS
            elif len(seq_check_results) != 0 and seq_check_results['per_id'] >= 50 and seq_check_results[
                'coverage'] >= 50:
                new_hit.get_flanking_genes(ref_gbk_obj, ref_feature_list)
                # mark it as a possible related IS, but confident
                new_hit.hit_type = 'possible related IS'
                new_hit.confidence_level = 'unpaired'
                new_hit.per_id = str(seq_check_results['per_id'])
                new_hit.coverage = str(seq_check_results['coverage'])
                # add it to the list
                IS_hit_list.append(new_hit)
            # otherwise this is a spurious result, remove
            else:
                # TODO: remove this hit
                pass
        # the gap is too small to be the IS, but larger than a novel hit
        elif float(gap) / is_query_length <= min_range and float(gap) / is_query_length < max_range:
            new_hit = get_orientation(intersect_left, intersect_right)
            # add the relevant information
            new_hit.hit_type = 'novel'
            new_hit.confidence_level = 'unpaired'
            new_hit.get_flanking_genes(ref_gbk_obj, ref_feature_list)
            # add it to the list
            IS_hit_list.append(new_hit)
        # otherwise remove!
        else:
            # TODO: remove this hit
            pass

    return IS_hit_list

def write_typing_output(IShits, output_table):

    with open(output_table, 'w') as out:

        # set the header and write it to the output file
        header = ["region", "orientation", "x", "y", "gap", "call", "percent_ID", "percent_cov", "left_gene", "left_description", "left_strand",
                  "left_distance", "right_gene", "right_description", "right_strand", "right_distance", "gene_interruption"]
        out.write('\t'.join(header) + '\n')

        # if there are no hits, record this and exit the function
        if len(IShits) == 0:
            out.write('No hits found')
            out.close()
            return

        # sort IS hits by left position, ascending order
        IShits.sort(key=lambda x: x.x)

        # loop through each hit
        region = 1
        for IShit in IShits:
            region_num = 'region_%s' % region
            region += 1
            call_type = IShit.hit_type
            if IShit.confidence_level == 'imprecise':
                call_type = call_type + '*'
            elif IShit.confidence_level == 'unpaired':
                call_type = call_type + '?'
            # calculate gap distance
            IShit.get_gap_distance()
            # determine if gene is interrupted or not
            IShit.determine_interrupted()
            # get qualifiers for left and right genes
            # TODO: make sure this qualifier call is robust
            left_description = IShit.left_feature.qualifiers['product'][0]
            right_description = IShit.right_feature.qualifiers['product'][0]
            # put together row
            line_list = [region_num, IShit.orientation, str(IShit.x), str(IShit.y), str(IShit.gap),
                         call_type, IShit.per_id, IShit.coverage, IShit.gene_left, left_description,
                         str(IShit.left_strand), str(IShit.left_distance), IShit.gene_right,
                         right_description, str(IShit.right_strand), str(IShit.right_distance),
                         str(IShit.interrupted)]
            # write out the information
            out.write('\t'.join(line_list) + '\n')

        # close the file and exit the function
        out.close()
        return

def create_typing_output(filenames, ref_gbk_obj, is_query_obj, min_range, max_range, tmp_output_folder):

    # first we need all the input files so we can match hits up
    intersect_file = filenames['intersect']
    closest_file = filenames['closest']
    left_unp = filenames['left_unpaired']
    right_unp = filenames['right_unpaired']
    #left_bedfile = filenames['left_merged_bed']
    #right_bedfile = filenames['right_merged_bed']

    # we also need to know the name of the table file where we'll write the final output to
    final_table_file = filenames['table']

    # final list of IS hit objects to make into a table at the end
    IS_hits = []

    # If both intersect and closest bed files are empty, there are no hits
    # write out an empty file and record this in the log file
    if os.stat(intersect_file)[6] == 0 and os.stat(closest_file)[6] == 0:
        write_typing_output(IS_hits, final_table_file)
        # TODO: add sample name here
        logging.info('No hits found for sample XX')
        return

    # If there are hits, read in the genbank file we're mapping to,
    # and create feature list for searching
    ref_feature_list = get_features(ref_gbk_obj)

    all_intersect_left = []
    all_intersect_right = []
    all_closest_left = []
    all_closest_right = []

    # Start with the intersect file (novel hits)
    if os.stat(intersect_file)[6] != 0:
        with open(intersect_file) as bed_intersect:
            # loop through each set of hits
            for line in bed_intersect:
                # extract the information
                info = line.strip().split('\t')
                # separate out info on the left and right sides of the hit
                intersect_left = [int(info[1]), int(info[2])]
                intersect_right = [int(info[4]), int(info[5])]
                # add this information to the master lists, for checking against unpaired hits later
                all_intersect_left.append(intersect_left)
                all_intersect_right.append(intersect_right)
                # get the gap between the hits, as determined by bedtools
                gap = int(info[6])
                # if the gap is small, then lets process this hit
                if gap <= 15:
                    # check if one hit is actually within the other hit
                    # if it is we need to remove it
                    # TODO: make this a function (check_hit_within_hit)
                    left_range = range(min(intersect_left), max(intersect_left))
                    right_range = range(min(intersect_right), max(intersect_right))
                    if (intersect_left[0] in right_range and intersect_left[1] in right_range) or (intersect_right[0] in left_range and intersect_right[1] in left_range):
                        # TODO: remove this hit
                        # set 'removed' variable to True
                        pass
                    # otherwise we need to process the hit
                    else:
                        # determine orientation and coordinates of hit
                        # process hit
                        new_hit = get_orientation(intersect_left, intersect_right)
                        # add the relevant information to the hit that we already know
                        new_hit.hit_type = 'novel'
                        new_hit.confidence_level = 'confident'
                        # make sure we note that this is overlapping because it's the intersect file
                        new_hit.overlap = True

                        # determine the features flanking the hit, and add the details to the hit object
                        new_hit.get_flanking_genes(ref_gbk_obj, ref_feature_list)
                        # TODO: determine gap, write a specific method for the class
                        #new_hit.get_gap_distance()
                        # append the hit to our list
                        IS_hits.append(new_hit)
                # If the gap is too big, we need to remove this hit
                else:
                    #TODO: remove this hit
                    pass
    # For the next section, grab the IS query length
    is_query_length = len(is_query_obj.seq)
    # Move on to the hits found in the closest bed file (known or imprecise hits)
    if os.stat(closest_file)[6] != 0:
        with open(closest_file) as bed_closest:
            # loop through each line
            for line in bed_closest:
                # extract all the information
                info = line.strip().split('\t')
                # if the fourth column contains a -1, there are no hits in this file
                if info[3] == '-1':
                    # exit the file and do not process further
                    break
                # get the distance between the hits
                gap = int(info[6])
                # separate out info on the left and right sides of the hit
                intersect_left = [int(info[1]), int(info[2])]
                intersect_right = [int(info[4]), int(info[5])]
                # add this information to the master lists, for checking against unpaired hits later
                all_closest_left.append(intersect_left)
                all_closest_right.append(intersect_right)
                # If the gap distance is 0, then this hit will be in the intersect file, so ignore
                if gap == 0:
                    pass
                # If the gap distance is small, this is likely a novel hit where no overlap was detected
                # TODO: make this gap size changeable in IS Mapper parameters
                elif gap <= 10:
                    new_hit = get_orientation(intersect_left, intersect_right)
                    # add the relevant information to the hit that we already know
                    new_hit.hit_type = 'novel'
                    new_hit.confidence_level = 'confident'
                    # determine the features flanking the hit, and add the details to the hit object
                    new_hit.get_flanking_genes(ref_gbk_obj, ref_feature_list)
                    # TODO: determine gap, write a specific method for the class
                    # new_hit.get_gap_distance()
                    # append the hit to our list
                    IS_hits.append(new_hit)
                # The gap size is within the range of the actual IS query size, and so probably indicates a known hit
                # Need to BLAST the sequence here to check it matches the IS
                elif float(gap) / is_query_length >= min_range and float(gap) / is_query_length <= max_range:
                    new_hit = get_orientation(intersect_left, intersect_right)
                    #genbank_seq, insertion, start, end, name, temp
                    seq_check_results = check_seq_between(ref_gbk_obj.seq, is_query_obj, new_hit.x,
                                                          new_hit.y, 'tmp_seq', tmp_output_folder)
                    # if it's a good hit, add it
                    if len(seq_check_results) != 0 and seq_check_results['per_id'] >= 80 and seq_check_results['coverage'] >= 80:
                        # get the flanking genes
                        new_hit.get_flanking_genes(ref_gbk_obj, ref_feature_list)
                        # make sure its a confident, novel hit
                        new_hit.hit_type = 'known'
                        new_hit.confidence_level = 'confident'
                        new_hit.per_id = str(seq_check_results['per_id'])
                        new_hit.coverage = str(seq_check_results['coverage'])
                        # add it to the list
                        IS_hits.append(new_hit)
                    # if the thresholds are low, then mark it as a possible related IS
                    elif len(seq_check_results) != 0 and seq_check_results['per_id'] >=50 and seq_check_results['coverage'] >= 50:
                        new_hit.get_flanking_genes(ref_gbk_obj, ref_feature_list)
                        # mark it as a possible related IS, but confident
                        new_hit.hit_type = 'possible related IS'
                        new_hit.confidence_level = 'confident'
                        new_hit.per_id = str(seq_check_results['per_id'])
                        new_hit.coverage = str(seq_check_results['coverage'])
                        # add it to the list
                        IS_hits.append(new_hit)
                    # otherwise this is a spurious result, remove
                    else:
                        # TODO: remove this hit
                        pass
                # The gap size here is smaller than the actual IS query, but larger than expected for a novel hit
                # This is an imprecise hit
                elif float(gap) / is_query_length <= min_range and float(gap) / is_query_length < max_range:
                    new_hit = get_orientation(intersect_left, intersect_right)
                    # add the relevant information
                    new_hit.hit_type = 'novel'
                    new_hit.confidence_level = 'imprecise'
                    new_hit.get_flanking_genes(ref_gbk_obj, ref_feature_list)
                    # add it to the list
                    IS_hits.append(new_hit)
                # This hit is way too big and doesn't fit any of the other criteria, so needs to be recorded as removed
                else:
                    # TODO: remove this hit
                    pass

    # Looking for unpaired hits which are not in the merged/closest bed files
    # Possibly unpaired because the pair is low coverage and didn't pass
    # depth cutoff

    # Start with the left hand unpaired file
    # FIRST, remove any positions which match regions we have ALREADY processed, by comparing to the left hand
    # master lists
    line_check = []
    with open(left_unp) as left_bed:
        for line in left_bed:
            info = line.strip().split('\t')
            left_coords = [int(info[1]), int(info[2])]
            if left_coords not in all_intersect_left and left_coords not in all_closest_left:
                line_check.append(line.strip().split('\t'))
    if len(line_check) != 0:
        all_new_hits = check_unpaired_hits(line_check, ref_gbk_obj, ref_feature_list, is_query_obj,
                                           min_range, max_range, tmp_output_folder)
        # add them to our current list
        IS_hits = IS_hits + all_new_hits

    # Then check the right hand unpaired file, again removing positions already processed
    line_check = []
    with open(right_unp) as right_bed:
        for line in right_bed:
            info = line.strip().split('\t')
            right_coords = [int(info[4]), int(info[5])]
            if right_coords not in all_intersect_right and right_coords not in all_closest_right:
                line_check.append(line.strip().split('\t'))
    if len(line_check) != 0:
        all_new_hits = check_unpaired_hits(line_check, ref_gbk_obj, ref_feature_list, is_query_obj,
                                           min_range, max_range, tmp_output_folder)

        # add them to our current list
        IS_hits = IS_hits + all_new_hits

    write_typing_output(IS_hits, final_table_file)

    return IS_hits
