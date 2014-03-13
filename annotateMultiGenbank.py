from Bio import SeqIO
from Bio import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqFeature
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from optparse import OptionParser
import string, re
import os, sys, subprocess
import collections

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	# required qsub options
	parser.add_option("-s", "--summary", action="store", dest="summary", help="text file output from BLAST, needs to be of format '6 qseqid qlen sacc pident length slen sstart send evalue bitscore'", default="")
	parser.add_option("-g", "--genbank", action="store", dest="genbank", help="original genbank file that features are going to be added to, should be multi-genbank", default="")
	parser.add_option("-i", "--prefix", action="store", dest="prefix", help="prefix for genbank created from fasta with added features", default="")
	parser.add_option("-n", "--newfile", action="store", dest="newfile", help="new filename for genbank having features added to", default="")
	parser.add_option("-t", "--genbank_type", action="store", dest="genbank_type", help="type of genbank file being created - multi or single, default is multi", default="multi")
	parser.add_option("-f", "--fasta", action="store", dest="fasta", help="fasta file to create multi-genbank from", default="")
	parser.add_option("-p", "--pid", action="store", dest="pid", help="minimum percent id of blast hit to be annotated (default 95)", default="95")
	parser.add_option("-c", "--qcov", action="store", dest="qcov", help="minimum coverage of the query to the reference to be annotated (default 100)", default="100")
	parser.add_option("-l", "--limit", action="store", dest="limit", help="distance from ends of contigs for start or end of contig value to be called" , default=300)

	return parser.parse_args()

def parseBLAST(hits, genbank_type):

	#opens file which contains locations, reads each line and saves that in a new variable
	summary = open(hits, "r")
	summary_list = summary.readlines()
	summary.close()

	hits = collections.defaultdict(dict)

	#lists for qualifiers
	start = []
	end = []
	percentID = []
	query_id = []
	blast_score = []
	hit_id = []
	query_length = []
	hit_length = []

	#append information to appropriate lists
	for columns in (raw.strip().split() for raw in summary_list):
		start.append(columns[6])
		end.append(columns[7])
		query_id.append(columns[0])
		percentID.append(columns[3])
		blast_score.append(columns[9])
		hit_id.append(columns[2])
		query_length.append(columns[1])
		hit_length.append(columns[4])

	#add values to dicionary
	if genbank_type == "multi":	
		for i in range(0, len(hit_id)):
			try:
				hits[hit_id[i]][query_id[i]] = [start[i], end[i], percentID[i], blast_score[i], query_length[i], hit_length[i]]
				#hits[node[i]] = [start[i], end[i], percentID[i], record_name[i], blast_score[i], query_length[i], hit_length[i]]
			except KeyError:
				pass
		return (hits)

	#different dictionary setup if genbank is a single type
	elif genbank_type == "single":

		count = 1

		for i in range(0, len(hit_id)):
			hits[hit_id[i]][query_id[i]] = [start[i], end[i], percentID[i], blast_score[i], query_length[i], hit_length[i]]
			count += 1

		return(hits)
	
	else:
		print("Genbank type not correctly specified, must be either single or multi.")


def createFeature(hits, record_id, limit, record_length = None):

	#set feature type
	if "_5_" in options.summary:
		feature_type = "5_prime_end"
	elif "_3_" in options.summary:
		feature_type = "3_prime_end"
	else:
		feature_type = "misc_feature"

	#get query coverage for hit
	queryCoverage = (float(hits[record_id][5])/float(hits[record_id][4])) * 100

	if record_length != None:
		#find out if hit is at the beginning or end of a contig based on a specified limit
		start_range = range(0, limit)
		end_region = record_length - limit
		end_range = range(end_region, record_length + 1)
		start_set = set(start_range)
		end_set = set(end_range)
		if len(start_set.intersection(end_set)) != 0:
			hit_location_value = "contig too small to call"

		else:
			if int(hits[record_id][0]) in start_range or int(hits[record_id][1]) in start_range:
				hit_location_value = "start of contig"
			elif int(hits[record_id][0]) in end_range or int(hits[record_id][1]) in end_range:
				hit_location_value = "end of contig"
			else:
				hit_location_value = "middle of contig"


	#check percent ID is at least minimum value set by user
	if float(hits[record_id][2]) >= float(options.pid) and queryCoverage >= float(options.qcov):
		quals = {}
		quals['note'] = "Node: " + record_id + " query length: " + hits[record_id][4] + " blast score: " + hits[record_id][3] + " query coverage: " + str(queryCoverage) + " percent ID: " + hits[record_id][2]
		if record_length != None:
			quals['location'] = hit_location_value

		#set Artemis colour to represent percent ID of hit
		quals['colour'] = artemisColour(hits[record_id][2])

		#make the feature
		new_location = SeqFeature.FeatureLocation(int(hits[record_id][0]), int(hits[record_id][1]))
		new_feature = SeqFeature.SeqFeature(new_location, type = feature_type, qualifiers = quals)

		return(new_feature)
	else:
		return(0)

def artemisColour(percentID):

	if float(percentID) <= 60:
		colour = "1"
	elif float(percentID) > 60 and float(percentID) <= 70:
		colour = "17"
	elif float(percentID) > 70 and float(percentID) <= 75:
		colour = "16"
	elif float(percentID) > 75 and float(percentID) <= 80:
		colour = "7"
	elif float(percentID) > 80 and float(percentID) <= 85:
		colour = "6"
	elif float(percentID) > 85 and float(percentID) <= 90:
		colour = "4"
	elif float(percentID) > 90 and float(percentID) <= 95:
		colour = "9"
	elif float(percentID) > 95 and float(percentID) <= 99:
		colour = "8"
	elif float(percentID) >= 99:
		colour = "3"

	return colour


if __name__ == "__main__":

	(options, args) = main()

	(fastaDir, fastaFileName) = os.path.split(options.fasta)

	#check for summary file
	if options.summary=="":
		DoError("No summary file provided (-s)")

	#if no genbank already, fasta only
	if options.fasta != "":

		#create multi entry genbank with entry contig as its own entry
		print("Creating multi entry genbank for annotation...")
		count = SeqIO.convert(options.fasta, "fasta", options.prefix + ".gbk", "genbank", generic_dna)

		print("Successfully converted %i records" % count)

		hits_dictionary = parseBLAST(options.summary, options.genbank_type)
		print hits_dictionary

		annotatedGenbank = options.newfile

		#open file to write
		handle = open(annotatedGenbank, "w")

		feature_count = 0

		#read in file
		records = SeqIO.parse(options.prefix + ".gbk", "genbank")
		record_list = []
		for record in records:
			record_list.append(record)

		new_record_list = []
		
		#iterate through the contigs in the genbank file
		for record in record_list:
			
			#check to see if that contig is in the genbank file
			if record.id in hits_dictionary:
				print record.id

				for node in hits_dictionary[record.id]:
					print hits_dictionary[record.id][node]
					record_length = len(record)

					#create the feature
					new_feature = createFeature(hits_dictionary[record.id], node, options.limit, record_length = record_length)	

					#check that the feature has passed the creation step
					if new_feature != 0:

						#add feature to record	
						record.features.append(new_feature)
						feature_count = feature_count + 1
					#just append the record if there was no feature to be added
					else:
						new_record_list.append(record)

				new_record_list.append(record)
				
			#just append the record if the id was not in the list
			else:
				new_record_list.append(record)

		SeqIO.write(new_record_list, handle, "genbank")

		handle.close()	

		print("Added " + str(feature_count) + " features to " + options.newfile)

	else:

		#get variables for annotations
		hits_dictionary = parseBLAST(options.summary, options.genbank_type)
		print hits_dictionary

		feature_count = 0

		if options.genbank_type == "multi":

			new_record_list = []

			record_list = SeqIO.parse(options.genbank, "genbank")
		
			#iterate through the contigs in the genbank file
			for record in record_list:
				
				#check to see if that contig is in the genbank file
				if record.id in hits_dictionary:
					print record.id

					for node in hits_dictionary[record.id]:
						print hits_dictionary[record.id][node]
						record_length = len(record)

						#create the feature
						new_feature = createFeature(hits_dictionary[record.id], node, options.limit, record_length = record_length)	

						#check that the feature has passed the creation step
						if new_feature != 0:

							#add feature to record	
							record.features.append(new_feature)
							feature_count = feature_count + 1
						#just append the record if there was no feature to be added
						else:
							new_record_list.append(record)

					new_record_list.append(record)
					
				#just append the record if the id was not in the list
				else:
					new_record_list.append(record)

			SeqIO.write(new_record_list, options.newfile, "genbank")

			print("Added " + str(feature_count) + " features to " + options.newfile)

		if options.genbank_type == "single":

			genbank = SeqIO.read(options.genbank, "genbank")

			for key in hits_dictionary:

				for node in hits_dictionary[key]:

					new_feature = createFeature(hits_dictionary[key], node, options.limit)

					if new_feature != 0:

						genbank.features.append(new_feature)

						feature_count += 1
				
			SeqIO.write(genbank, options.newfile, "genbank")
			print("Added " + str(feature_count) + " features to " + options.newfile)






