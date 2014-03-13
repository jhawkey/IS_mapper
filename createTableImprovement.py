# read in genbank file, print out coordinates & strand of features
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

def parse_args():

	parser = ArgumentParser(description="create a table of features for the is mapping pipeline")
	parser.add_argument('--genbank', type=str, required=True, help='genbank file to look for features in')
	parser.add_argument('--output', type=str, required=True, help='file to store table in')
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

def extractContigs(genbank):

	contig_dict = {}

	record = SeqIO.read(genbank, "genbank")

	for region in record.features:
		if region.type == "fasta_record":
			start = region.location.nofuzzy_start
			end = region.location.nofuzzy_end
			contig_name = region.qualifiers['note'][0]
			contig_dict[contig_name] = [start, end]

	return contig_dict


def collapseRanges(ranges, gap, prime_end):
	(start, stop) = ranges.pop(0)
	new_ranges = [(start,stop)] # initialise with first tuple (start, stop)
	ranges_dict = {}
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
	
	for i in new_ranges:
		ranges_dict[i[0]] = [i[1], prime_end]
					
	return ranges_dict

def withinContig(ranges, contigs):

	lines_dict = {}
	count = 1
	lines_list = []

	for start in ranges:
		for key in contigs:
			if start >= contigs[key][0] and ranges[start][0] <= contigs[key][1]:
				lines_dict["region_" + str(count)] = [start, ranges[start][0], key, ranges[start][1]]
				count = count + 1

	for key in lines_dict:
		lines_list.append([key, lines_dict[key][0], lines_dict[key][1], lines_dict[key][2], lines_dict[key][3]])

	lines_list.sort(key=itemgetter(1))

	return lines_list
	
def main():
	args = parse_args()

	#get the features from the genbank
	five_ranges = extractFeatures(args.genbank, "5_prime_end")
	three_ranges = extractFeatures(args.genbank, "3_prime_end")
	contigs = extractContigs(args.genbank)

	#combine hits next to each other into one hit
	five_rangesNew = collapseRanges(five_ranges, 500, "5 prime")
	three_rangesNew = collapseRanges(three_ranges, 500, "3 prime")

	new_ranges_dict = dict(three_rangesNew.items() + five_rangesNew.items())
	lines_list = withinContig(new_ranges_dict, contigs)
	count = 1	

	output = open(args.output, "w")	

	output.write("\t".join(["region", "orientation", "start", "end", "contig", "\n"]))
	for i in lines_list:
		output.write("\t".join(["region_" + str(count), i[4], str(i[1]), str(i[2]), str(i[3]) + "\n"]))
		count = count + 1

	output.close()

if __name__ == "__main__":
	main()