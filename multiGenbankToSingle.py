import string, re
import os, sys, subprocess
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i", "--in", action="store", dest="input", help="input file (gbk from RAST)", default="")
	parser.add_option("-t", "--tags", action="store", dest="tags", help="prefix for locus_tag (if not specified, none is created)", default="")
	parser.add_option("-n", "--name", action="store", dest="name", help="name for gbk header", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="name for gbk output")
	
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
	
	total = 0 # total bases
	
	handle = open(options.input, "rU")
	
	records = list(SeqIO.parse(handle, "genbank"))
	
	feature_count = 0

	colour_count = 0

	if len(options.name) >= 10:
		options.name = options.name[:9]

	for r in records:
		length = len(r)
		id = r.name
		seq = r.seq
		seq.alphabet=generic_dna

		if total > 0:
			newrecord.seq = newrecord.seq + seq
		else:
			# first sequence, initialise seqrecord
			newrecord = SeqRecord(seq=r.seq,name=options.name,id=options.name)
			newrecord.seq.alphabet=generic_dna
		
		# create feature for contig
		if colour_count % 2 == 0:
			newrecord.features.append(SeqFeature(FeatureLocation(total, total + length), type="fasta_record", qualifiers = {'note' : [r.name], 'colour':'11'}))
			colour_count = colour_count + 1
		else:
			newrecord.features.append(SeqFeature(FeatureLocation(total, total + length), type="fasta_record", qualifiers = {'note' : [r.name], 'colour':'10'}))
			colour_count = colour_count + 1

		# copy CDS features
		for f in r.features:
			feature_count += 1
			f.qualifiers["locus_tag"] = options.tags + str(feature_count)
			newrecord.features.append(SeqFeature(FeatureLocation(f.location.nofuzzy_start + total, f.location.nofuzzy_end + total), strand = f.strand, type=f.type, qualifiers = f.qualifiers))

		total += length

	handle.close()
	
	SeqIO.write(newrecord, options.output, "genbank")