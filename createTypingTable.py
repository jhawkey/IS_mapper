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
    parser.add_argument('--genbank', type=str, required=True, help='genbank file to look for features in')
    parser.add_argument('--insertion', type=str, required=True, help='path to insertion sequence fasta file for BLAST hits')
    parser.add_argument('--temp', type=str, required=False, help='path to temp folder for storing temporary BLAST files if needed')
    parser.add_argument('--blast_gap', type=float, required=False, default=300, help='number of bases between two hits to BLAST the middle and see what it is')
    parser.add_argument('--reference_genbank', type=str, required=True, help='reference genbank file to find flanking genes of regions')
    parser.add_argument('--cds', type=str, required=False, default='locus_tag,gene,product', help='qualifiers to look for in reference genbank for CDS features')
    parser.add_argument('--trna', type=str, required=False, default='locus_tag,product', help='qualifiers to look for in reference genbank for tRNA features')
    parser.add_argument('--rrna', type=str, required=False, default='locus_tag,product', help='qualifiers to look for in reference genbank for rRNA features')
    parser.add_argument('--output', type=str, required=True, help='name for output file (extension .txt already given)')
    return parser.parse_args()

def extractFeatures(genbank, feature_name):

    lines = []
    ranges = []

    #read in the genbank
    record = SeqIO.read(genbank, "genbank")

    #for each feature, if the feature type matches the feature_name, then process
    for region in record.features:
        if region.type in feature_name:
            #get start and end
            start = region.location.nofuzzy_start
            end = region.location.nofuzzy_end

            #make sure start and end are in correct order (lowest number first), don't care about strandedness
            new_start = min(start, end)
            new_end = max(start, end)

            #add this to the list of lines
            lines.append([int(new_start),int(new_end),region.type])
            tuplerange = (int(new_start), int(new_end))
            ranges.append(tuplerange)

    #sort those features from lowest base to highest base
    ranges.sort()

    return ranges

def collapseRanges(ranges, gap):
    (start, stop) = ranges.pop(0)
    new_ranges = [(start,stop)] # initialise with first tuple (start, stop)
    for (start,stop) in ranges:
        overlap_starts = []
        overlap_stops = []
        toDrop = []
        for i in range(0,len(new_ranges)):
            (x,y) = new_ranges[i]
            if start in range(x-gap,y+1) or stop in range(x,y+1+gap):
                overlap_starts.append(min(x,start)) # overlaps with this feature
                overlap_stops.append(max(y,stop))
                toDrop.append(i) # record overlapping features to remove later
        if len(overlap_starts) > 0:
            toDrop.reverse()
            for i in toDrop:
                del new_ranges[i] # delete overlapping features from the ranges list
            new_start = min(overlap_starts+overlap_stops)
            new_stop = max(overlap_starts+overlap_stops)
            new_ranges.append((new_start,new_stop)) # add the new merged feature
        else:
            new_ranges.append((start,stop)) # if no overlap, add as a new range
    return new_ranges

def parseBLAST(hits):

    #opens file which contains locations, reads each line and saves that in a new variable
    summary = open(hits, "r")
    summary_list = summary.readlines()
    summary.close()

    hits = {}

    #lists for qualifiers
    start = []
    end = []
    percentID = []
    node = []
    blast_score = []
    record_name = []
    query_length = []
    hit_length = []

    #append information to appropriate lists
    for columns in (raw.strip().split() for raw in summary_list):
        start.append(columns[6])
        end.append(columns[7])
        node.append(columns[0])
        percentID.append(columns[3])
        blast_score.append(columns[9])
        record_name.append(columns[2])
        query_length.append(columns[1])
        hit_length.append(columns[4])

    #add values to dicionary
    for i in range(0, len(record_name)):
        try:
            hits[node[i]] = [start[i], end[i], percentID[i], record_name[i], blast_score[i], query_length[i], hit_length[i]]
        except KeyError:
            pass
    return hits
    
def insertionLength(insertion):

    sequence = SeqIO.read(insertion, "fasta")
    length = len(sequence.seq)

    return length

def pairHits(five_ranges, three_ranges, seqLength, genbank, blast_gap, output_file):

    record = SeqIO.read(genbank, 'genbank')
    output = open(output_file, 'w')
    paired_hits = {}
    table_keys = []
    found_threes = []
    count = 1

    for i in range(0, len(five_ranges)):
        distances = []
        for l in range(0, len(three_ranges)):
            if five_ranges[i][0] < three_ranges[l][0]:
                distance = three_ranges[l][0] - five_ranges[i][1]
                distances.append(distance)
            else:
                distance = five_ranges[i][0] - three_ranges[l][1]
                distances.append(distance)
        if min(distances) < (2 * seqLength):
            correct_index = distances.index(min(distances))
            if five_ranges[i][0] > three_ranges[correct_index][0]:
                orientation = "R"
                gap = min(distances)
                x = min(five_ranges[i][0], three_ranges[correct_index][1])
                y = max(five_ranges[i][0], three_ranges[correct_index][1])
                paired_hits['region_' + str(count)] = [orientation, str(x), str(y), gap]
            elif five_ranges[i][0] < three_ranges[correct_index][0]:
                orientation = "F"
                gap = min(distances)
                x = min(five_ranges[i][1], three_ranges[correct_index][0])
                y = max(five_ranges[i][1], three_ranges[correct_index][0])
                paired_hits['region_' + str(count)] = [orientation, str(x), str(y), gap]
            # append the 3' hit into this list so all unpaired 3's can be found later
            found_threes.append(three_ranges[correct_index])
            count += 1
        else:
            # an unpaired 5'
            # in this case gap (z) is 0, and x is the first coord and y is the second coord
            paired_hits['region_' + str(count)] = ["5' unpaired", str(five_ranges[i][0]), str(five_ranges[i][1]), '-']
            seq_before = record[five_ranges[i][0] - seqLength:five_ranges[i][0]]
            seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id='region_' + str(count) + '_before')
            seq_after = record[five_ranges[i][1]:five_ranges[i][1] + seqLength]
            seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id='region_' + str(count) + '_after')
            SeqIO.write(seq_before, output, 'fasta')
            SeqIO.write(seq_after, output, 'fasta')
            count += 1

    for value in three_ranges:
        if value not in found_threes:
            # an unpaired 3'
            # in this case gap (z) is again 0, and x is the first coord and y is the second coord
            paired_hits['region_' + str(count)] = ["3' unpaired", str(value[0]), str(value[1]), '-']
            seq_before = record[value[0] - seqLength:value[0]]
            seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id='region_' + str(count) + '_before')
            seq_after = record[value[1]:value[1] + seqLength]
            seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id='region_' + str(count) + '_after')
            SeqIO.write(seq_before, output, 'fasta')
            SeqIO.write(seq_after, output, 'fasta')
            count += 1

    for key in paired_hits:
        if 'F' in paired_hits[key] or 'R' in paired_hits[key]:
            seq_between = record.seq[int(paired_hits[key][1]):int(paired_hits[key][2])]
            seq_between = SeqRecord(Seq(str(seq_between), generic_dna), id=key)
            if len(seq_between) >= blast_gap:
                SeqIO.write(seq_between, output, 'fasta')

    return paired_hits

def unpairedHits(ranges, seqLength, genbank, output_file, orientation):

    record = SeqIO.read(genbank, 'genbank')
    output = open(output_file, 'w')
    count = 1
    hits = {}
    for i in range(0, len(ranges)):
        hits['region_' + str(count)] = [orientation, str(ranges[i][0]), str(ranges[i][1]), '']
        seq_before = record[ranges[i][0] - seqLength:ranges[i][1]]
        seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id = 'region_' + str(count) + '_before')
        seq_after = record[ranges[i][1]:ranges[i][1] + seqLength]
        seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id='region_' + str(count) + '_after')
        SeqIO.write(seq_before, output, 'fasta')
        SeqIO.write(seq_after, output, 'fasta')
        count += 1
    return hits

def createTable(table, blast_results, insertionSeqLength, blast_gap, output, reference, cds_quals, trna_quals, rrna_quals):

    # add the percent ID and query coverage for the BLAST hits to the table
    if blast_results != 0:
        for i in blast_results:
            # caclulate percent ID and coverage
            percentID = float(blast_results[i][2])
            queryCoverage = (float(blast_results[i][6])/float(blast_results[i][5])) * 100
            hitLength = blast_results[i][5]
            # this is for all the between hits
            if i in table:
                # if the hits are right next to each other or overlapping, report as novel
                #if int(table[i][3]) <= 0 or int(table[i][3]) <= blast_gap:
                    #table[i].extend(['Novel', '', ''])

                # if percent ID and cov meet a certain threshold, annotate as known, otherwise annotate as unknown
                if percentID >= 80 and queryCoverage >= 60:
                    table[i].extend(['Known', str(percentID), str('%.2f' % queryCoverage)])
                else:
                    table[i].extend(['Unknown', str(percentID), str('%.2f' % queryCoverage)])
            # this is for the unpaired hits where the before or after sequence has been taken
            if 'before' in i:
                # get the region number
                region_no = i.split('_')[1]
                # if the BLAST hit is above a certain id and cov, then anntoate as known
                if percentID >= 80 and queryCoverage >= 60:
                    table['region_' + region_no].extend(['Known before', str(percentID), str('%.2f' % queryCoverage)])
                # otherwise unknown, but give location of lower quality BLAST hit
                else:
                    table['region_' + region_no].extend(['Unknown: positioned before BLAST hit', str(percentID), str('%.2f' % queryCoverage)])

            elif 'after' in i:
                # get the region number
                region_no = i.split('_')[1]
                # if the BLAST hit meets the threshold requirements annotate as known
                if percentID >= 80 and queryCoverage >= 60:
                    table['region_' + region_no].extend(['Known after', str(percentID), str('%.2f' % queryCoverage)])
                # otherwise unknown, but give location of lower quality BLAST hit
                else:
                    table['region_' + region_no].extend(['Unknown: positioned after BLAST hit', str(percentID), str('%.2f' % queryCoverage)])
            else:
                pass
            
    # go through the keys and find these in table so it's printed in order
    table_keys = []
    for key in table:
        table_keys.append(key)
    region_indexes = []
    for region in table_keys:
        region_indexes.append(region.split('region_')[1])
    arr = np.vstack((table_keys, region_indexes)).transpose()
    sorted_keys = arr[arr[:,1].astype('int').argsort()]
    
    #header = ["region", "orientation", "hit start", "IS start", "IS end", "hit end", "length of IS region", "percent ID to IS", "coverage of region to IS", "call"]
    header = ["region", "orientation", "x", "y", "gap", "call", "%id", "%cov", "left_gene", "left_strand", "left_distance", "right_gene", "right_strand", "right_distance", "functional_prediction"]
    print "\t".join(header)
    output.write('\t'.join(header) + '\n')
    
    for key in sorted_keys[:,0]:
        #print key
        #print type(key)
        #print table[key]
        # if region already has a call, write it out
        gene_left, gene_left_dist, gene_right, gene_right_dist = get_flanking_genes(reference, int(table[key][1]), int(table[key][2]), cds_quals, trna_quals, rrna_quals)
        if len(table[key]) > 4:
            table[key].extend([gene_left[:-1], gene_left[-1], gene_left_dist, gene_right[:-1], gene_right[-1], gene_right_dist])
            print key + '\t' + '\t'.join(str(i) for i in table[key])
            output.write(key + '\t' + '\t'.join(str(i) for i in table[key]) + '\n')
        elif  'F' in table[key] or 'R' in table[key]:
            if table[key][3] <= blast_gap:
                table[key].extend(['Novel', '', ''])
                table[key].extend([gene_left[:-1], gene_left[-1], gene_left_dist, gene_right[:-1], gene_right[-1], gene_right_dist])
                print key + '\t' + '\t'.join(str(i) for i in table[key])
                output.write(key + '\t' + '\t'.join(str(i) for i in table[key]) + '\n')
            else:
                table[key].extend(['Unknown (no BLAST hit)', '', ''])
                table[key].extend([gene_left[:-1], gene_left[-1], gene_left_dist, gene_right[:-1], gene_right[-1], gene_right_dist])
                print key + '\t' + '\t'.join(str(i) for i in table[key])
                output.write(key + '\t' + '\t'.join(str(i) for i in table[key]) + '\n')
        elif "5' unpaired" in table[key] or "3' unpaired" in table[key]:
            table[key].extend(['Unknown (no BLAST hit before or after)', '', ''])
            table[key].extend([gene_left[:-1], gene_left[-1], gene_left_dist, gene_right[:-1], gene_right[-1], gene_right_dist])
            print key + '\t' + '\t'.join(str(i) for i in table[key])
            output.write(key + '\t' + '\t'.join(str(i) for i in table[key]) + '\n')

def get_other_gene(features, pos, cds_features, trna_features, rrna_features, direction):

    distance = {}
    for feature in features:
        if feature.type == 'CDS' or feature.type == 'tRNA' or feature.type == 'rRNA':
            if pos in feature.location:
                values = []
                if feature.type == 'CDS':
                    for qual in cds_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                elif feature.type == 'tRNA':
                    for qual in trna_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                elif feature.type == 'rRNA':
                   for qual in rrna_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                # get the strand
                values.append(feature.strand)
                gene_distance = feature.location.start - pos
                return values, gene_distance
            else:
                dist = feature.location.start - pos
                if dist > 0 and direction == 'downstream':
                    values = []
                    if feature.type == 'CDS':
                        for qual in cds_features:
                            try:
                                values.append(feature.qualifiers[qual][0])
                            except KeyError:
                                pass
                    elif feature.type == 'tRNA':
                        for qual in trna_features:
                            try:
                                values.append(feature.qualifiers[qual][0])
                            except KeyError:
                                pass
                    elif feature.type == 'rRNA':
                       for qual in rrna_features:
                            try:
                                values.append(feature.qualifiers[qual][0])
                            except KeyError:
                                pass
                    # get the strand
                    values.append(feature.strand)
                    distance[dist] = values
                elif dist < 0 and direction == 'upstream':
                    values = []
                    if feature.type == 'CDS':
                        for qual in cds_features:
                            try:
                                values.append(feature.qualifiers[qual][0])
                            except KeyError:
                                pass
                    elif feature.type == 'tRNA':
                        for qual in trna_features:
                            try:
                                values.append(feature.qualifiers[qual][0])
                            except KeyError:
                                pass
                    elif feature.type == 'rRNA':
                       for qual in rrna_features:
                            try:
                                values.append(feature.qualifiers[qual][0])
                            except KeyError:
                                pass
                    # get the strand
                    values.append(feature.strand)
                    distance[dist] = values

    distance_keys = list(OrderedDict.fromkeys(distance))
    gene = distance[min(distance_keys)]
    gene_distance = min(distance_keys)
    return gene, gene_distance



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
                values = []
                # get qualifiers of interest for this feature
                if feature.type == 'CDS':
                    for qual in cds_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                elif feature.type == 'tRNA':
                    for qual in trna_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                elif feature.type == 'rRNA':
                   for qual in rrna_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                # get the strand
                values.append(feature.strand) 
                gene_left = values
                gene_left_distance = feature.location.start - pos_x
                gene_right = values
                gene_right_distance = feature.location.end - pos_y

                # exit the function 
                return gene_left, gene_left_distance, gene_right, gene_right_distance
            
            # if x inside gene but y is not, need to report gene x is in gene further down genome from y
            elif pos_x in feature.location and pos_y not in feature.location:
                values = []
                if feature.type == 'CDS':
                    for qual in cds_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                elif feature.type == 'tRNA':
                    for qual in trna_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                elif feature.type == 'rRNA':
                   for qual in rrna_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                # get the strand
                values.append(feature.strand) 
                gene_left = values
                gene_left_distance = feature.location.start - pos_x

                # then go and find the next closest gene that y is next to
                gene_right, gene_right_distance = get_other_gene(gb.features, pos_y, cds_features, trna_features, rrna_features, 'downstream')

                # exit function
                return gene_left, gene_left_distance, gene_right, gene_right_distance

            elif pos_y in feature.location and pos_x not in feature.location:
                values = []
                if feature.type == 'CDS':
                    for qual in cds_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                elif feature.type == 'tRNA':
                    for qual in trna_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
                elif feature.type == 'rRNA':
                   for qual in rrna_features:
                        try:
                            values.append(feature.qualifiers[qual][0])
                        except KeyError:
                            pass
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
            values = []
            # get qualifiers of interest for feature
            if feature.type == 'CDS':
                for qual in cds_features:
                    try:
                        values.append(feature.qualifiers[qual][0])
                    except KeyError:
                        pass
            elif feature.type == 'tRNA':
                for qual in trna_features:
                    try:
                        values.append(feature.qualifiers[qual][0])
                    except KeyError:
                        pass
            elif feature.type == 'rRNA':
                for qual in rrna_features:
                    try:
                        values.append(feature.qualifiers[qual][0])
                    except KeyError:
                        pass
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
    gene_left_dist = min(distance_lkeys)
    gene_right_dist = min(distance_rkeys)

    return gene_left, gene_left_dist, gene_right, gene_right_dist

def doBlast(blast_input, blast_output, database):
    #perform BLAST
    blastn_cline = NcbiblastnCommandline(query=blast_input, db=database, outfmt="'6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs'", out=blast_output)
    stdout, stderr = blastn_cline()

def main():

    args = parse_args() 

    #find the length of the insertion sequence
    insertionSeqLength = insertionLength(args.insertion)

    #create the prefix of the file which will contain sequences for blast and then the blast output
    if args.temp == None:
        region_blast_fasta = os.path.split(args.genbank)[1].split('.gbk')[0]
    else:
        region_blast_fasta = args.temp + os.path.split(args.genbank)[1].split('.gbk')[0]

    output_file = open(args.output + '.txt', 'w')
    #get the features from the genbank
    five_ranges = extractFeatures(args.genbank, '5_prime_end')
    three_ranges = extractFeatures(args.genbank, '3_prime_end')

    #combine hits next to each other into one hit
    if five_ranges == [] and three_ranges != []:
        three_rangesNew = collapseRanges(three_ranges, 400)
        unpaired_three = unpairedHits(three_rangesNew, insertionSeqLength, args.genbank, region_blast_fasta + '_3only.fasta', "3' unpaired")
        doBlast(region_blast_fasta + '_3only.fasta', region_blast_fasta + '_3only.txt', args.insertion)
        blast_results = parseBLAST(region_blast_fasta + '_3only.txt')
        createTable(unpaired_three, blast_results, insertionSeqLength, args.blast_gap, output_file)
    elif three_ranges == [] and five_ranges != []:
        five_rangesNew = collapseRanges(five_ranges, 400)
        unpaired_five = unpairedHits(five_ranges, insertionSeqLength, args.genbank, region_blast_fasta + '_5only.fasta', "5' unpaired")
        doBlast(region_blast_fasta + '_5only.fasta', region_blast_fasta + '_5only.txt', args.insertion)
        blast_results = parseBLAST(region_blast_fasta + '_5only.txt')
        createTable(unpaired_five, blast_results, insertionSeqLength, args.blast_gap, output_file)
    elif five_ranges != [] and three_ranges != []:
        five_rangesNew = collapseRanges(five_ranges, 400)
        three_rangesNew = collapseRanges(three_ranges, 400)
        # work out which hits pair together and return the correct indexes and the group that have the most number of hits (both even if all paired)
        table = pairHits(five_rangesNew, three_rangesNew, insertionSeqLength, args.genbank, args.blast_gap, region_blast_fasta + '.fasta')
        if os.stat(region_blast_fasta + '.fasta')[6] != 0:
            doBlast(region_blast_fasta + '.fasta', region_blast_fasta + '.txt', args.insertion)
            # parse the BLAST output 
            blast_results = parseBLAST(region_blast_fasta + '.txt')
        else:
            blast_results = 0
        createTable(table, blast_results, insertionSeqLength, args.blast_gap, output_file, args.reference_genbank, args.cds, args.trna, args.rrna)
    else:
        print "\t".join(["region", "orientation", "hit start", "IS start", "IS end", "hit end", "length of IS region", "percent ID to IS", "coverage of region to IS", "call"])
        output_file.write('No hits found')
        print "No hits found"
    output_file.close()

if __name__ == "__main__":
    main()