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
	parser.add_argument('--genbank', type=str, required=False, help='genbank file to look for features in')
	parser.add_argument('--insertion', type=str, required=False, help='path to insertion sequence fasta file for BLAST hits')
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

def pairHits(first_ranges, second_ranges):

	correct_indexes = []

	if len(first_ranges) > len(second_ranges):
		longest_ranges = first_ranges
		shorter_ranges = second_ranges
	else:
		longest_ranges = second_ranges
		shorter_ranges = first_ranges

	#for each 5' hit
	for i in shorter_ranges:
		#only look at each hit in the 3' hits once
		values = len(longest_ranges)
		count = 0
		#track distances between each 5' and 3' possible pair
		distances = []

		while count < values:
			#calculate the distance between, taking the absolute value (and getting orientation)
			if longest_ranges[count][0] > i[1]: #therefore 5' to 3' orientation
				distance = abs(longest_ranges[count][0] - i[1])
			else:
				distance = abs(i[0] - longest_ranges[count][1]) #threfore 3' to 5' orientation
			#append the distance
			distances.append(distance)
			count = count + 1

		#the correct pair is the one with the smallest distance between them
		correct_indexes.append([distances.index(min(distances)), shorter_ranges.index(i)])

	return correct_indexes, longest_ranges

def createTableLines(five_ranges, three_ranges, paired_indexes, genbank, insertion, output_file):

	count = 1
	table = {}
	table_keys = []

	output = open(output_file, "w")

	paired_hits = len(paired_indexes)

	insertionSeqLength = insertionLength(insertion)

	record = SeqIO.read(genbank, "genbank")

	if len(five_ranges) > len(three_ranges):
		longest_ranges = five_ranges
		shorter_ranges = three_ranges
		normal_orient = "3' to 5'"
		reverse_orient = "5' to 3'"
	else:
		longest_ranges = three_ranges
		shorter_ranges = five_ranges
		normal_orient = "5' to 3'"
		reverse_orient = "3' to 5'"

	while count <= paired_hits:
		for i in paired_indexes:
			if longest_ranges[i[0]][0] > shorter_ranges[i[1]][1]:
				start = shorter_ranges[i[1]][1]
				end = longest_ranges[i[0]][0]
				table["region_" + str(count)] = [normal_orient, str(shorter_ranges[i[1]][0]), str(start), str(end), str(longest_ranges[i[0]][1])]
				table_keys.append("region_" + str(count))
			else:
				start = longest_ranges[i[0]][1]
				end = shorter_ranges[i[1]][0]
				table["region_" + str(count)] = [reverse_orient, str(longest_ranges[i[0]][0]), str(start), str(end), str(shorter_ranges[i[1]][1])]
				table_keys.append("region_" + str(count))

			seq_between = record.seq[start:end]
			seq_between = SeqRecord(Seq(str(seq_between), generic_dna), id="region_" + str(count))
			if len(seq_between) > 0:
				SeqIO.write(seq_between, output, "fasta")
			table["region_" + str(count)].append(str(len(seq_between)))
			#table["region_" + str(count)].append("")
			#table["region_" + str(count)].append("")

			count = count + 1

	new_line = []

	region_no = count
	count = count - 1  


	if count < len(five_ranges):
		while count > paired_hits and count <= len(five_ranges):
			for i in five_ranges:
				boolean = []
				for values in table:
					boolean.append(str(i[0]) in table[values] or str(i[1]) in table[values])
				if True not in boolean:
					start = i[0]
					end = i[1]
					new_line.append(["5' unpaired", str(start), str(end), "region_" + str(region_no)])
					table_keys.append("region_" + str(region_no))
					
					seq_before = record[start - insertionSeqLength:start]
					seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id="region_" + str(region_no) + "_before")
					SeqIO.write(seq_before, output, "fasta")

					#extract seq after for blasting
					seq_after = record[end:end + insertionSeqLength]
					seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id="region_" + str(region_no) + "_after")
					SeqIO.write(seq_after, output, "fasta")
					region_no = region_no + 1
					count = count + 1

	if count < len(three_ranges):
		while count >= paired_hits and count <= len(three_ranges):
			for i in three_ranges:
				boolean = []
				for values in table:
					boolean.append(str(i[0]) in table[values] or str(i[1]) in table[values])
				if True not in boolean:
					start = i[0]
					end = i[1]
					new_line.append(["3' unpaired", str(start), str(end), "region_" + str(region_no)])
					table_keys.append("region_" + str(region_no))
					
					seq_before = record[start - insertionSeqLength:start]
					seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id="region_" + str(region_no) + "_before")
					SeqIO.write(seq_before, output, "fasta")

					#extract seq after for blasting
					seq_after = record[end:end + insertionSeqLength]
					seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id="region_" + str(region_no) + "_after")
					SeqIO.write(seq_after, output, "fasta")
					region_no = region_no + 1
				
				count = count + 1

	for i in range(0, len(new_line)):
		table[new_line[i][3]] = [new_line[i][0], new_line[i][1], new_line[i][2], "", "", str(insertionSeqLength)]

	output.close()
	return table, table_keys

def createUnpairedTableLine(hit_ranges, count, paired_hits, region_no, insertion_length, table, table_keys, fasta_file):
	#CURRENTLY NOT A FUNCTIONAL FUNCTION

	new_line = []

	output = open(fasta_file, "rU")

	while count >= paired_hits and count <= len(hit_ranges):
		for i in hit_ranges:
			boolean = []
			for values in table:
				boolean.append(str(i[0]) in table[values] or str(i[1]) in table[values])
			if True not in boolean:
				start = i[0]
				end = i[1]
				new_line.append(["3' unpaired", str(start), str(end), "region_" + str(region_no)])
				table_keys.append("region_" + str(region_no))
				
				seq_before = record[start - insertion_length:start]
				seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id="region_" + str(region_no) + "_bef")
				SeqIO.write(seq_before, output, "fasta")

				#extract seq after for blasting
				seq_after = record[end:end + insertion_length]
				seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id="region_" + str(region_no) + "_aft")
				SeqIO.write(seq_after, output, "fasta")
				region_no = region_no + 1
				
			count = count + 1


def main():

	args = parse_args()	

	#get the features from the genbank
	five_ranges = extractFeatures(args.genbank, "5_prime_end")
	three_ranges = extractFeatures(args.genbank, "3_prime_end")

	#combine hits next to each other into one hit
	five_rangesNew = collapseRanges(five_ranges, 300)
	three_rangesNew = collapseRanges(three_ranges, 300)

	fasta_region = open("test_regions.fasta", "w")

	indexes, longest_ranges = pairHits(five_rangesNew, three_rangesNew)

	table, table_keys = createTableLines(five_rangesNew, three_rangesNew, indexes, args.genbank, args.insertion, "test_regions.fasta")

	#perform BLAST
	blastn_cline = NcbiblastnCommandline(query="test_regions.fasta", db=args.insertion, outfmt="'6 qseqid qlen sacc pident length slen sstart send evalue bitscore'", out="test_regions.txt")
	stdout, stderr = blastn_cline()

	dictionary = (parseBLAST("test_regions.txt"))

	insertionSeqLength = insertionLength(args.insertion)

	#add the percent ID and query coverage for the blast hits to the table
	for i in dictionary:
		percentID = float(dictionary[i][2])
		queryCoverage = (float(dictionary[i][6])/float(dictionary[i][5])) * 100
		hitLength = dictionary[i][5]
		if i in table or i[:-4] in table:
			table[i].append(str(percentID))
			table[i].append(str(queryCoverage))

	#go through the keys and find these in table so it's printed in order
	header = ["region", "orientation", "hit start", "IS start", "IS end", "hit end", "length of IS region", "percent ID to IS", "coverage of region to IS", "call"]
	print "\t".join(header)
	for key in table_keys:
		try:
			if float(table[key][5]) <= 100 and insertionSeqLength >= 100:
				print key + "\t", "\t".join(table[key]) + "\tNovel insertion site"
			#if the hit is a good percentage and good coverage, count it as known
			elif float(table[key][6]) >= 90 and float(table[key][7]) >= 90:
				print key + "\t", "\t".join(table[key]) + "\tKnown insertion site"
			else:
				print key + "\t", "\t".join(table[key]) + "\tUnknown"

		#this is for all the unpaired hits
		except IndexError:
			print key + "\t", "\t".join(table[key]) + "\t \t \tUnknown"


if __name__ == "__main__":
	main()

	



