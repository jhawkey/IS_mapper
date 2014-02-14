import logging
import sys, re, os
from argparse import (ArgumentParser, FileType)

def parse_args():
	'Parse the input arguments, use -h for help'

	parser = ArgumentParser(description='IS mapper')

	# need to add verison info later
	#parser.add_argument("--version", action='version', ...)

	parser.add_argument('--reads', nargs = '+', type = str, required=True, help='Paired end reads for analysing (can be gzipped)')
	parser.add_argument('--reference', type = str, required=True, help='Fasta file for reference gene (eg: insertion sequence) that will be mapped to')
	parser.add_argument('--assemblies', nargs='+', type=str, required=False, help='Contig assemblies, one for each read set')
	parser.add_argument('--typingRef', type=str, required=False, help='Reference genome for typing against')
	parser.add_argument('--coverage', type=float, required=False, default=90.0, help='Minimum coverage for hit to be annotated (default 90.0)')
	parser.add_argument('--percentid', type=float, required=False, default=90.0, help='Minimum percent ID for hit to be annotated (default 90.0')
	parser.add_argument('--log', action="store_true", required=False, help='Switch on logging to file (otherwise log to stdout')

	# Do I need this?
	parser.add_argument('--output', type=str, required=True, help='Location to store output files')

# Exception to raise if the command we try to run fails for some reason
class CommandError(Exception):
	pass

def run_command(command, **kwargs):
	'Execute a shell command and check the exit status and any O/S exceptions'
	command_str = ' '.join(command)
	logging.info('Running: {}'.format(command_str))
	try:
		exit_status = call(command, **kwargs)
	except OSError as e:
		message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
		raise CommandError({"message": message})
	if exit_status != 0:
		message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
		raise CommandError({"message": message})


# Change this so it uses BWA
def bowtie_index(fasta_files):
	'Build a bowtie2 index from the given input fasta(s)'

	# check that both bowtie and samtools have the right versions
	check_command_version(['bowtie2', '--version'],
				'bowtie2-align version 2.1.0',
				'bowtie',
				'2.1.0')

	for fasta in fasta_files:
		built_index = fasta + '.1.bt2'
		if os.path.exists(built_index):
			logging.info('Index for {} is already built...'.format(fasta))
		else:
			logging.info('Building bowtie2 index for {}...'.format(fasta))
			run_command(['bowtie2-build', fasta, fasta])


# Check that an acceptable version of a command is installed
# Exits the program if it can't be found.
# - command_list is the command to run to determine the version.
# - version_identifier is the unique string we look for in the stdout of the program.
# - command_name is the name of the command to show in error messages.
# - required_version is the version number to show in error messages.
def check_command_version(command_list, version_identifier, command_name, required_version):
	try:
		command_stdout = check_output(command_list, stderr=STDOUT)
	except OSError as e:
		logging.error("Failed command: {}".format(' '.join(command_list)))
		logging.error(str(e))
		logging.error("Could not determine the version of {}.".format(command_name))
		logging.error("Do you have {} installed in your PATH?".format(command_name))
		exit(-1)
	except CalledProcessError as e:
		# some programs such as samtools return a non-zero exit status
		# when you ask for the version (sigh). We ignore it here.
		command_stdout = e.output

	if version_identifier not in command_stdout:
		logging.error("Incorrect version of {} installed.".format(command_name))
		logging.error("{} version {} is required by SRST2.".format(command_name, required_version))
		exit(-1)