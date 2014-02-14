import logging
import sys, re, os
from argparse import (ArgumentParser, FileType)

def parse_args():
	'Parse the input arguments, use -h for help'

	parser = ArgumentParser(description='IS mapper')

	# need to add verison info later
	#parser.add_argument("--version", action='version', ...)

	parser.add_argument('--reads', nargs = '', type = str, required=True, help='Paired end reads for analysing (can be gzipped)')
	parser.add_argument('--reference', nargs='', type = str, required=True, help='Fasta file for reference gene (eg: insertion sequence) that will be mapped to')
	parser.add_argument('--assemblies', nargs='', type=str, required=False, help='Contig assembly')
	parser.add_argument('--typingRef', nargs='', type=str, required=False, help='Reference genome for typing against')
	parser.add_argument('--coverage', nargs='', type=str, required=False, default='90', help='Minimum coverage for hit to be annotated (default 90)')
	parser.add_argument('--percentid', nargs='', type=str, required=False, default='90', help='Minimum percent ID for hit to be annotated (default 90')
	parser.add_argument('--log', action="store_true", required=False, help='Switch on logging to file (otherwise log to stdout')

	# Do I need this?
	parser.add_argument('--output', type=str, required=True, help='Location to store output files')