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
	parser.add_argument('--temp', type=str, required=False, help='path to temp folder for storing temporary BLAST files if needed')
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

def pairHits(first_ranges, second_ranges, seqLength):

	correct_indexes = []

	if len(first_ranges) >= len(second_ranges):
		longest_ranges = first_ranges
		shorter_ranges = second_ranges
	else:
		longest_ranges = second_ranges
		shorter_ranges = first_ranges

	#for each 5' hit
	for i in shorter_ranges:
		#only look at each hit in the shortest set once
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
		if min(distances) < 2 * seqLength:
			correct_indexes.append([distances.index(min(distances)), shorter_ranges.index(i)])
		else:
			pass

	return correct_indexes, longest_ranges

def createTableLines(five_ranges, three_ranges, paired_indexes, genbank, insertion, output_file):

	#set up count (becomes region number keeps track of what has already been looked at), and empty dict and list for the table info
	count = 1
	table = {}
	table_keys = []

	#open the fasta file for writing sequences into that will be used in the BLAST step
	output = open(output_file, "w")

	#get the number of paired hits
	paired_hits = len(paired_indexes)

	#get the length of the insertion sequence
	insertionSeqLength = insertionLength(insertion)

	#open the genbank so sequences between or on the ends of hits can be pulled out
	record = SeqIO.read(genbank, "genbank")

	#set up orientaitons, 'normal' and 'reverse' orientations differ depending on which set of ranges is the longest
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

	#loop through the paired hits
	while count <= paired_hits:
		#look at each set of indexes
		for i in paired_indexes:
			#print "Count value " + str(count)
			#print "Number of paired hits: " + str(paired_hits)
			
			#work out orientaiton based on the start and end coordinates of the hits
			#give each region it's own unique identifier that is also saved in the keys list so this can be iterated through later
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

			#pull out the sequence between the hits and save it in the fasta file for later blasting
			seq_between = record.seq[start:end]
			seq_between = SeqRecord(Seq(str(seq_between), generic_dna), id="region_" + str(count))
			if len(seq_between) > 0:
				SeqIO.write(seq_between, output, "fasta")
			table["region_" + str(count)].append(str(len(seq_between)))
			#table["region_" + str(count)].append("")
			#table["region_" + str(count)].append("")

			#update the count value
			count = count + 1

	#setup empty list for unpaired hits
	new_line = []

	#set the region number to the same value as the count variable
	region_no = count  

	#print "Count value now: " + str(count)
	#print "Length of five_ranges: " + str(len(five_ranges))

	#go through the five end hits and compare the start and end coordinates to everything in the dict. If it's always false,
	#then the hit is not paired, and not already in the table so must be an unpaired hit
	for i in five_ranges:
		boolean = []
		for values in table:
			boolean.append(str(i[0]) in table[values] or str(i[1]) in table[values])
		if True not in boolean:
			start = i[0]
			end = i[1]
			new_line.append(["5' unpaired", str(start), str(end), "region_" + str(region_no)])
			table_keys.append("region_" + str(region_no))
			
			#extract the sequence before and after the hit and add it to the fasta file (IS could be on either end)
			seq_before = record[start - insertionSeqLength:start]
			seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id="region_" + str(region_no) + "_before")
			SeqIO.write(seq_before, output, "fasta")

			seq_after = record[end:end + insertionSeqLength]
			seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id="region_" + str(region_no) + "_after")
			SeqIO.write(seq_after, output, "fasta")
			region_no = region_no + 1
			count = count + 1

	#do the same for the three prime hits
	for i in three_ranges:
		boolean = []
		for values in table:
			boolean.append(str(i[0]) in table[values] or str(i[1]) in table[values])
		if True not in boolean:
			start = i[0]
			end = i[1]
			new_line.append(["3' unpaired", str(start), str(end), "region_" + str(region_no)])
			table_keys.append("region_" + str(region_no))
			
			#extract the sequence before and after the hit and add it to the fasta file (IS could be on either end)
			seq_before = record[start - insertionSeqLength:start]
			seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id="region_" + str(region_no) + "_before")
			SeqIO.write(seq_before, output, "fasta")

			seq_after = record[end:end + insertionSeqLength]
			seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id="region_" + str(region_no) + "_after")
			SeqIO.write(seq_after, output, "fasta")
			region_no = region_no + 1
			count = count + 1

	#add all these new lines to the table
	for i in range(0, len(new_line)):
		table[new_line[i][3]] = [new_line[i][0], new_line[i][1], new_line[i][2], "", "", str(insertionSeqLength)]

	output.close()

	return table, table_keys

def main():

	args = parse_args()	

	#find the length of the insertion sequence
	insertionSeqLength = insertionLength(args.insertion)

	#get the features from the genbank
	five_ranges = extractFeatures(args.genbank, "5_prime_end")
	three_ranges = extractFeatures(args.genbank, "3_prime_end")
	#print 'five_ranges'
	#print five_ranges
	#print 'three_ranges'
	#print three_ranges

	#combine hits next to each other into one hit
	five_rangesNew = collapseRanges(five_ranges, 300)
	three_rangesNew = collapseRanges(three_ranges, 300)

	#create the prefix of the file which will contain sequences for blast and then the blast output
	if args.temp == "":
		region_blast_fasta = os.path.split(args.genbank)[1].split('.gbk')[0]
	else:
		region_blast_fasta = args.temp + os.path.split(args.genbank)[1].split('.gbk')[0]

	#work out which hits pair together and return the correct indexes and the group that have the most number of hits (both even if all paired)
	indexes, longest_ranges = pairHits(five_rangesNew, three_rangesNew, insertionSeqLength)

	#return a dictionary with all the information for each region and a list that gives you the keys used in that dictionary
	table, table_keys = createTableLines(five_rangesNew, three_rangesNew, indexes, args.genbank, args.insertion, region_blast_fasta + ".fasta")

	#perform BLAST
	blastn_cline = NcbiblastnCommandline(query=region_blast_fasta + ".fasta", db=args.insertion, outfmt="'6 qseqid qlen sacc pident length slen sstart send evalue bitscore qcovs'", out=region_blast_fasta + ".txt")
	stdout, stderr = blastn_cline()

	#parse the BLAST output 
	blast_results = (parseBLAST(region_blast_fasta + ".txt"))

	#add the percent ID and query coverage for the blast hits to the table
	for i in blast_results:

		#caclulate percent ID and coverage
		percentID = float(blast_results[i][2])
		queryCoverage = (float(blast_results[i][6])/float(blast_results[i][5])) * 100
		hitLength = blast_results[i][5]
		#this is for all the between hits
		if i in table or i[:-4] in table:
			table[i].append(str(percentID))
			table[i].append(str("%.2f" % queryCoverage))
		#this is for the unpaired hits where the before or after sequence has been taken
		if "before" in i or "after" in i:
			region_no = i.split('_')[1]
			table["region_" + region_no].append(str(percentID))
			table["region_" + region_no].append(str("%.2f" % queryCoverage))
			if i.split('_')[2] == "before":
				table["region_" + region_no].append("before")
			else:
				table["region_" + region_no].append("after")

	#go through the keys and find these in table so it's printed in order
	header = ["region", "orientation", "hit start", "IS start", "IS end", "hit end", "length of IS region", "percent ID to IS", "coverage of region to IS", "call"]
	print "\t".join(header)
	for key in table_keys:
		if "before" not in table[key] and "after" not in table[key]:
			try:
				#if the hits are right next to each other and there is little or no sequence in between, report it as novel	
				if float(table[key][5]) == 0 or float(table[key][5]) <= (insertionSeqLength * 0.25):
					print key + "\t", "\t".join(table[key]) + "\t \t \tNovel insertion site"
				
				#if the sequence between has good ID and coverage to the IS in question, report it as known
				elif float(table[key][6]) >= 80 and float(table[key][7]) >= 60:
					print key + "\t", "\t".join(table[key]) + "\tKnown insertion site"
					print(table[key][6])
				#Otherwise report as unknown
				else:
					print key + "\t", "\t".join(table[key]) + "\tUnknown"
			except (KeyError, IndexError):
				print key + "\t", "\t".join(table[key]) + "\t \t \tUnknown: no BLAST hit before or after"
		if "after" in table[key] or "before" in table[key]:
			try:
				if float(table[key][6]) >= 80 and float(table[key][7]) >= 60:
					print key + "\t", "\t".join(table[key][:-1]) + "\tKnown insertion site " + table[key][8]
				else:
					print key + "\t", "\t".join(table[key][:-1]) + "\tUnknown: positioned " + table[key][8] + " BLAST hit"
			except IndexError:
				print key + "\t", "\t".join(table[key][:-1]) + "\tUnknown: no BLAST hit before or after"

if __name__ == "__main__":
	main()