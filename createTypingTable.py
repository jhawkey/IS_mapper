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
import numpy as np

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

def pairHits(five_ranges, three_ranges, seqLength, genbank, output_file):

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
			if five_ranges[i][1] > three_ranges[correct_index][0]:
				orientation = "3' to 5'"
				paired_hits['region_' + str(count)] = [orientation, str(three_ranges[correct_index][0]), str(three_ranges[correct_index][1]), str(five_ranges[i][0]), str(five_ranges[i][1]), str(abs(min(distances)))]
				#if five_ranges[i] in five_ranges:
				#	five_ranges.remove(five_ranges[i])
				#if three_ranges[correct_index] in three_ranges:
				#	three_ranges.remove(three_ranges[correct_index])
			else:
				orientation = "5' to 3'"
				paired_hits['region_' + str(count)] = [orientation, str(five_ranges[i][0]), str(five_ranges[i][1]), str(three_ranges[correct_index][0]), str(three_ranges[correct_index][1]), str(abs(min(distances)))]
				#if five_ranges[i] in five_ranges:
				#	five_ranges.remove(five_ranges[i])
				#if three_ranges[correct_index] in three_ranges:
				#	three_ranges.remove(three_ranges[correct_index])
			# append the 3' hit into this list so all unpaired 3's can be found later
			found_threes.append(three_ranges[correct_index])
			count += 1
		else:
			#an unpaired 5'
			paired_hits['region_' + str(count)] = ["5' unpaired", str(five_ranges[i][0]), str(five_ranges[i][1]), '', '', '']
			seq_before = record[five_ranges[i][0] - seqLength:five_ranges[i][0]]
			seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id='region_' + str(count) + '_before')
			seq_after = record[five_ranges[i][1]:five_ranges[i][1] + seqLength]
			seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id='region_' + str(count) + '_after')
			SeqIO.write(seq_before, output, 'fasta')
			SeqIO.write(seq_before, output, 'fasta')
			count += 1

	for value in three_ranges:
		if value not in found_threes:
			paired_hits['region_' + str(count)] = ["3' unpaired", str(value[0]), str(value[1]), '', '', '']
			seq_before = record[value[0] - seqLength:value[0]]
			seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id='region_' + str(count) + '_before')
			seq_after = record[value[1]:value[1] + seqLength]
			seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id='region_' + str(count) + '_after')
			SeqIO.write(seq_before, output, 'fasta')
			SeqIO.write(seq_before, output, 'fasta')
			count += 1

	for key in paired_hits:
		if "3' to 5" in paired_hits[key] or "3' to 5'" in paired_hits[key]:
			seq_between = record.seq[int(paired_hits[key][2]):int(paired_hits[key][3])]
			seq_between = SeqRecord(Seq(str(seq_between), generic_dna), id=key)
			if len(seq_between) >= (seqLength * 0.25):
				SeqIO.write(seq_between, output, 'fasta')

	return paired_hits
def unpairedHits(ranges, seqLength, genbank, output_file, orientation):
	record = SeqIO.read(genbank, 'genbank')
	output = open(output_file, 'w')
	count = 1
	hits = {}
	for i in range(0, len(ranges)):
		hits['region_' + str(count)] = [orientation, str(ranges[i][0]), str(ranges[i][1]), '', '', '']
		seq_before = record[ranges[i][0] - seqLength:ranges[i][1]]
		seq_before = SeqRecord(Seq(str(seq_before.seq), generic_dna), id = 'region_' + str(count) + '_before')
		seq_after = record[ranges[i][1]:ranges[i][1] + seqLength]
		seq_after = SeqRecord(Seq(str(seq_after.seq), generic_dna), id='region_' + str(count) + '_after')
		SeqIO.write(seq_before, output, 'fasta')
		SeqIO.write(seq_after, output, 'fasta')
		count += 1
	return hits
def createTable(table, blast_results, insertionSeqLength):

	#add the percent ID and query coverage for the blast hits to the table
	if blast_results != 0:
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
			if "before" in i:
				if percentID > 80 and queryCoverage > 60:
					region_no = i.split('_')[1]
					table["region_" + region_no].append(str(percentID))
					table["region_" + region_no].append(str("%.2f" % queryCoverage))
					table["region_" + region_no].append("before")
			elif "after" in i:
				if percentID > 80 and queryCoverage > 60:
					region_no = i.split('_')[1]
					table["region_" + region_no].append(str(percentID))
					table["region_" + region_no].append(str("%.2f" % queryCoverage))
					table["region_" + region_no].append("after")
			else:
				pass
			
	table_keys = []
	for key in table:
		table_keys.append(key)
	region_indexes = []
	for region in table_keys:
		region_indexes.append(region.split('region_')[1])
	arr = np.vstack((table_keys, region_indexes)).transpose()
	sorted_keys = arr[arr[:,1].astype('int').argsort()]
	#go through the keys and find these in table so it's printed in order
	header = ["region", "orientation", "hit start", "IS start", "IS end", "hit end", "length of IS region", "percent ID to IS", "coverage of region to IS", "call"]
	print "\t".join(header)
	for key in sorted_keys[:,0]:
		if "5' unpaired" not in table[key] and "3' unpaired" not in table[key]:
			try:
				#if the hits are right next to each other and there is little or no sequence in between, report it as novel	
				if float(table[key][5]) == 0 or float(table[key][5]) <= (insertionSeqLength * 0.25):
					print key + "\t", "\t".join(table[key]) + "\t \t \tNovel insertion site"
				#if the sequence between has good ID and coverage to the IS in question, report it as known
				elif float(table[key][6]) >= 80 and float(table[key][7]) >= 60:
					print key + "\t", "\t".join(table[key]) + "\tKnown insertion site"
				#Otherwise report as unknown
				else:
					print key + "\t", "\t".join(table[key]) + "\tUnknown"
			except (KeyError, IndexError):
				print key + "\t", "\t".join(table[key]) + "\t \t \tUnknown: no BLAST hit before or after"
		if "5' unpaired" in table[key] or "3' unpaired" in table[key]:
			try:
				if float(table[key][6]) >= 80 and float(table[key][7]) >= 60:
					print key + "\t", "\t".join(table[key][:-1]) + "\tKnown insertion site " + table[key][8]
				else:
					print key + "\t", "\t".join(table[key][:-1]) + "\tUnknown: positioned " + table[key][8] + " BLAST hit"
			except IndexError:
				print key + "\t", "\t".join(table[key][:-1]) + "\tUnknown: no BLAST hit before or after"

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

	#get the features from the genbank
	five_ranges = extractFeatures(args.genbank, '5_prime_end')
	three_ranges = extractFeatures(args.genbank, '3_prime_end')

	#combine hits next to each other into one hit
	if five_ranges == [] and three_ranges != []:
		three_rangesNew = collapseRanges(three_ranges, 400)
		unpaired_three = unpairedHits(three_rangesNew, insertionSeqLength, args.genbank, region_blast_fasta + '_3only.fasta', "3' unpaired")
		doBlast(region_blast_fasta + '_3only.fasta', region_blast_fasta + '_3only.txt', args.insertion)
		blast_results = parseBLAST(region_blast_fasta + '_3only.txt')
		createTable(unpaired_three, blast_results, insertionSeqLength)
	elif three_ranges == [] and five_ranges != []:
		five_rangesNew = collapseRanges(five_ranges, 400)
		unpaired_five = unpairedHits(five_ranges, insertionSeqLength, args.genbank, region_blast_fasta + '_5only.fasta', "5' unpaired")
		doBlast(region_blast_fasta + '_5only.fasta', region_blast_fasta + '_5only.txt', args.insertion)
		blast_results = parseBLAST(region_blast_fasta + '_5only.txt')
		createTable(unpaired_five, blast_results, insertionSeqLength)
	elif five_ranges != [] and three_ranges != []:
		five_rangesNew = collapseRanges(five_ranges, 400)
		three_rangesNew = collapseRanges(three_ranges, 400)
		#work out which hits pair together and return the correct indexes and the group that have the most number of hits (both even if all paired)
		table = pairHits(five_rangesNew, three_rangesNew, insertionSeqLength, args.genbank, region_blast_fasta + '.fasta')
		if os.stat(region_blast_fasta + '.fasta')[6] != 0:
			doBlast(region_blast_fasta + '.fasta', region_blast_fasta + '.txt', args.insertion)
			#parse the BLAST output 
			blast_results = parseBLAST(region_blast_fasta + '.txt')
		else:
			blast_results = 0
		createTable(table, blast_results, insertionSeqLength)
	else:
		print "\t".join(["region", "orientation", "hit start", "IS start", "IS end", "hit end", "length of IS region", "percent ID to IS", "coverage of region to IS", "call"])
		print "No hits found"

if __name__ == "__main__":
	main()