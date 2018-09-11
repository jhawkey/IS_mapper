#!/usr/bin/env python

import os, logging
import operator


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
        # left and right coordinates of the hit
        self.left_pos = left_pos
        self.right_pos = right_pos
        # will be either forward or reverse
        self.orientation = None

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

def get_flanking_genes(genbank_obj, feature_list, is_hit):
    '''
    :param genbank_obj: The parsed genbank object
    :param feature_list: The parsed features list from the genbank file for easy searching
    :param is_hit: The IS hit object
    :return: The modified IS hit object containing the features of the left and right genes, their strand, and
    the distance from the IS hit to each gene
    '''

    # Find the correct indexes
    left_feature_index = binary_search(feature_list, is_hit.left_pos, 'L')
    right_feature_index = binary_search(feature_list, is_hit.right_pos, 'R')
    # Print out information if returning one of the error codes
    if type(left_feature_index) != int or type(right_feature_index) != int:
        print('left index')
        print(left_feature_index)
        print('right index')
        print(right_feature_index)
        print('left position: ' + str(is_hit.left_pos))
        print('right position: ' + str(is_hit.right_pos))
    # Extract the SeqFeature object that corresponds to that index
    left_feature = genbank_obj.features[left_feature_index]
    right_feature = genbank_obj.features[right_feature_index]

    # Add the SeqFeatures to the IS hit object, and strand information
    is_hit.left_feature = left_feature
    is_hit.right_feature = right_feature
    is_hit.left_strand = left_feature.strand
    is_hit.right_strand = right_feature.strand

    # The info we require is:
    # [geneid, distance, [locus_tag, (gene), product, strand]]
    #left_values = get_qualifiers(cds_features, trna_features, rrna_features, left_feature)
    #right_values = get_qualifiers(cds_features, trna_features, rrna_features, right_feature)

    # Get the distances from the IS hit to the left and right genes

    # The distance to the left gene is the endmost position of the feature - the left IS coord
    left_dist = abs(max(left_feature.location.start, left_feature.location.end) - is_hit.left_pos)
    # The distance to the right gene is the startmost position of the feature - the right IS coord
    right_dist = abs(min(right_feature.location.start, right_feature.location.end) - is_hit.right_pos)

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
        dist_to_start = is_hit.left_pos
        # Distance from the end of the final gene to the end of the
        # genome
        dist_to_end = abs(left_feature.location.end - genome_size)
        # So the total distance is those two added together
        left_dist = dist_to_start + dist_to_end

    elif right_dist in range(int(round(genome_size * 0.9)), int(round(genome_size * 1.1))):
        # Then the right hand feature is at the start of the genome
        # Distance from the IS position to the end of the genome
        dist_to_end = abs(genome_size - is_hit.right_pos)
        # Distance from the start of the genome to the start of the first feature
        # is the start position of the first feature
        dist_to_feature = right_feature.location.start
        # So the total distance is those two added together
        right_dist = dist_to_end + dist_to_feature

    is_hit.left_distance = left_dist
    is_hit.right_distance = right_dist

    return is_hit

def binary_search(features, isPosition, direction):
    min = 0
    max = len(features) - 1

    while True:

        # If the min has exceeded the max, then the IS position is not
        # inside a feature, and m will now be pointing to a
        # feature next to the IS position.
        if max < min:
            if direction == 'R':
                return findFeatureAfterPosition(features, isPosition, m)
            else:
                return findFeatureBeforePosition(features, isPosition, m)

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

def findFeatureBeforePosition(features, isPosition, m):
    # If we are looking for the feature to the left of the
    # IS position, then either m-1 or m is our answer

    # If the start of the m feature is after the IS position,
    # then m is after the IS and m-1 is the correct feature
    if features[m][0] > isPosition:
        return features[m-1][2]

    # If both m and m+1 features are before the IS position,
    # then m will be closer to the IS and is the correct feature
    elif features[m-1][1] < isPosition and features[m][1] < isPosition:
        return features[m][2]

    else:
        return "2 - THIS SHOULDN'T HAPPEN!"

def findFeatureAfterPosition(features, isPosition, m):
    # If we are looking for the feature to the right of the
    # IS position, then either m or m+1 is our answer

    # an index error will occur if m is the final feature, so just check that the first part is true
    # and return m
    try:
        features[m+1]
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
        return features[m+1][2]

    # If both m and m+1 features are after the IS position,
    # then m will be closer to the IS and is the correct feature
    elif features[m][0] > isPosition and features[m+1][0] > isPosition:
        return features[m][2]
    else:
        return "3 - THIS SHOULDN'T HAPPEN!"

def novel_hit(is_hit_obj, genbank_features):
    '''

    :param is_hit_obj:
    :param genbank_features:
    :return:
    '''

    # The hit is novel



def create_typing_output(filenames, ref_gbk_obj, is_query_obj, min_range, max_range):

    # first we need all the input files so we can match hits up
    intersect_file = filenames['intersect']
    closest_file = filenames['closest']
    #left_unp = filenames['left_unpaired']
    #right_unp = filenames['right_unpaired']
    #left_bedfile = filenames['left_merged_bed']
    #right_bedfile = filenames['right_merged_bed']

    # we also need to know the name of the table file where we'll write the final output to
    #final_table_file = filenames['table']

    # this is the dictionary of results
    results = {}
    # these are the results which are being excluded based on some criteria
    removed_results = {}
    # we're starting at the first region, line 0
    region = 1
    lines = 0
    # this is the header of the table output
    header = ["region", "orientation", "x", "y", "gap", "call", "Percent_ID", "Percent_Cov", "left_gene", "left_strand", "left_distance", "right_gene", "right_strand", "right_distance", "functional_prediction"]

    # If both intersect and closest bed files are empty, there are no hits
    if os.stat(intersect_file)[6] == 0 and os.stat(closest_file)[6] == 0:
        output = open(final_table_file, 'w')
        output.write('\t'.join(header) + '\n')
        output.write('No hits found')
        output.close()
        # TODO: add sample name here
        logging.info('No hits found for sample XX')
        return('No hits')

    # Now we need to read in the genbank file we're mapping to,
    # and create feature list for searching
    ref_feature_list = get_features(ref_gbk_obj)

    # Initialise feature count
    feature_count = 0

    intersect_left = []
    intersect_right = []
    closest_left = []
    closest_right = []

    # final list of IS hit objects
    IS_hits = []

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
                # get the gap between the hits, as determined by bedtools
                gap = int(info[6])
                # if the gap is small, then lets process this hit
                if gap <= 15:
                    # check if one hit is actually within the other hit
                    # if it is we need to remove it
                    left_range = range(min(intersect_left), max(intersect_left))
                    right_range = range(min(intersect_right), max(intersect_right))
                    if (intersect_left[0] in right_range and intersect_left[1] in right_range) or (intersect_right[0] in left_range and intersect_right[1] in left_range):
                        # TODO: remove this hit
                        pass
                    # otherwise we need to process the hit
                    else:
                        # determine orientation and coordinates of hit
                        # process hit
                        if intersect_left[0] < intersect_right[0] or intersect_left[1] < intersect_right[1]:
                            orientation = 'F'
                            new_hit = ISHit(intersect_left[1], intersect_right[0])
                        else:
                            orientation = 'R'
                            new_hit = ISHit(intersect_right[1], intersect_left[0])
                        # add the relevant information to the hit that we already know
                        new_hit.hit_type = 'novel'
                        new_hit.confidence_level = 'confident'
                        new_hit.orientation = orientation

                        # determine the features flanking the hit, and add the details to the hit object
                        new_hit = get_flanking_genes(ref_gbk_obj, ref_feature_list, new_hit)
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
                # If the gap distance is 0, then this hit will be in the intersect file, so ignore
                if gap == 0:
                    pass
                # If the gap distance is small, this is likely a novel hit where no overlap was detected
                # TODO: make this gap size changeable in IS Mapper parameters
                elif gap <= 10:
                    pass
                # The gap size is within the range of the actual IS query size, and so probably indicates a known hit
                # Need to BLAST the sequence here to check it matches the IS
                elif float(gap) / is_query_length >= min_range and float(gap) / is_query_length <= max_range:
                    pass
                # The gap size here is smaller than the actual IS query, but larger than expected for a novel hit
                # This is an imprecise hit
                elif float(gap) / is_query_length <= min_range and float(gap) / is_query_length < max_range:
                    pass
                # This hit is way too big and doesn't fit any of the other criteria, so needs to be recorded as removed
                else:
                    # TODO: remove this hit
                    pass

    return is_hits
