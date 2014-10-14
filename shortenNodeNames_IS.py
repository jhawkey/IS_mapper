#!/usr/bin/env python

import string, re
import os, sys, subprocess
from optparse import OptionParser

def shortendNode(inputfile, outputfile):

	# read in file
	a = file(inputfile, "r")
	o = file(outputfile, "w")
	
	#for each line, remove newline and break into contigs
	for line in a:
		line = line.rstrip()
		f = line.split('>')
		
		#take id line in fasta and split on _, then create new id line which contains
		#node number only
		if len(f) > 1:
			originalID = f[1]
			IDsegments = line.split('_')
			o.write(str(IDsegments[0] + '_' + IDsegments[1] + "\n"))
		else:
			o.write(line + "\n")
	a.close()
	o.close()


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i", "--input", action="store", dest="infile", help="input contig files", default="")
	parser.add_option("-o", "--output", action="store", dest="outfile", help="name of output file with new contig names", default="")
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
	
	if options.infile=="":
		DoError("No input file specified")
	if options.outfile=="":
		DoError("No output file specified")      
	
	else:
		shortendNode(options.infile, options.outfile)