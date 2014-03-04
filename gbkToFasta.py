import string, re
import os, sys, subprocess
from optparse import OptionParser
from Bio import SeqIO
from Bio.Alphabet import generic_dna

def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i", "--in", action="store", dest="input", help="input gbk", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="output fasta")
	
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()

	input_handle = open(options.input, "rU")
	output_handle = open(options.output, "w")

	sequences = SeqIO.parse(input_handle, "genbank")
	SeqIO.write(sequences, output_handle, "fasta")

	input_handle.close()
	output_handle.close()

	