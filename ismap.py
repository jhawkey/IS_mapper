import logging
import sys, re, os
from argparse import (ArgumentParser, FileType)
from subprocess import call, check_output, CalledProcessError, STDOUT
from Bio.Blast.Applications import NcbiblastnCommandline

def parse_args():
    '''
    Parse the input arguments, use -h for help.
    '''

    parser = ArgumentParser(description='IS mapper')

    # need to add verison info later
    #parser.add_argument("--version", action='version', ...)

    parser.add_argument('--reads', nargs = '+', type = str, required=True, help='Paired end reads for analysing (can be gzipped)')
    parser.add_argument('--forward', type = str, required=False, default = '_1', help = 'Identifier for forward reads if not in MiSeq format (default _1)')
    parser.add_argument('--reverse', type=str, required=False, default='_2', help='Identifier for reverse reads if not in MiSeq format (default _2)')
    parser.add_argument('--reference', type = str, required=True, help='Fasta file for reference gene (eg: insertion sequence) that will be mapped to')
    parser.add_argument('--assemblies', nargs='+', type=str, required=False, help='Contig assemblies, one for each read set')
    parser.add_argument('--assemblyid', type=str, required=False, help='Identifier for assemblies eg: sampleName_contigs (specify _contigs) or sampleName_assembly (specify _assembly). If there is no extension leave blank. Do not specify the . extension, eg .fasta or .fa')
    parser.add_argument('--type', type=str, required=False, default='fasta', help='Indicator for contig assembly type, genbank or fasta (default fasta)')
    parser.add_argument('--extension', type=str, required=False, default='_contigs', help='Identifier for assemblies (default _contigs')
    parser.add_argument('--typingRef', type=str, required=False, help='Reference genome for typing against')
    parser.add_argument('--coverage', type=float, required=False, default=90.0, help='Minimum coverage for hit to be annotated (default 90.0)')
    parser.add_argument('--percentid', type=float, required=False, default=90.0, help='Minimum percent ID for hit to be annotated (default 90.0')
    parser.add_argument('--log', action="store_true", required=False, help='Switch on logging to file (otherwise log to stdout')
    parser.add_argument('--runtype', type=str, required=True, help='"typing" or "improvement"')

    # Do I need this?
    parser.add_argument('--output', type=str, required=True, help='Path to location for output files')

    return parser.parse_args()

# Exception to raise if the command we try to run fails for some reason
class CommandError(Exception):
    pass

def run_command(command, **kwargs):
    '''
    Execute a shell command and check the exit status and any O/S exceptions.
    '''
    
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

def bwa_index(fasta):
    '''
    Check to see if bwa index for given input fasta exists.
    If it doesn't, build an index from the given input fasta.
    '''
    
    #check_command_version('bwa')
    built_index = fasta + '.bwt'
    print built_index
    if os.path.exists(built_index):
        logging.info('Index for {} is already built...'.format(fasta))
    else:
        logging.info('Building bwa index for {}...'.format(fasta))
        run_command(['bwa', 'index', fasta])

def check_command_version(command_list, version_identifier, command_name, required_version):
    '''
    Check that an acceptable version of a command is installed.
    Exits the program if it can't be found.
        - command_list is the command to run to determine the version.
        - version_identifier is the unique string we look for in the stdout of the program.
        - command_name is the name of the command to show in error messages.
        - required_version is the version number to show in error messages.
    '''

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
        logging.error("{} version {} is required by ISmapper.".format(command_name, required_version))
        exit(-1)

def get_readFile_components(full_file_path):

    (file_path,file_name) = os.path.split(full_file_path)
    m1 = re.match("(.*).gz",file_name)
    ext = ""
    if m1 != None:
        # gzipped
        ext = ".gz"
        file_name = m1.groups()[0]
    (file_name_before_ext,ext2) = os.path.splitext(file_name)
    full_ext = ext2+ext

    return(file_path,file_name_before_ext,full_ext)

def read_file_sets(args):   

    fileSets = {} # key = id, value = list of files for that sample
    num_paired_readsets = 0

    # paired end
    forward_reads = {} # key = sample, value = full path to file
    reverse_reads = {} # key = sample, value = full path to file
    assemblies = {}
    num_paired_readsets = 0
    num_single_readsets = 0
    num_assemblies = 0

    for fastq in args.reads:
        (file_path,file_name_before_ext,full_ext) = get_readFile_components(fastq)
        # try to match to MiSeq format:
        m = re.match("(.*)(_S.*)(_L.*)(_R.*)(_.*)", file_name_before_ext)
        if m == None:
            # not default Illumina file naming format, expect simple/ENA format
            m = re.match("(.*)(" + args.forward + ")$",file_name_before_ext)
            if m != None:
                # store as forward read
                (baseName,read) = m.groups()
                forward_reads[baseName] = fastq
            else:
                m = re.match("(.*)(" + args.reverse + ")$",file_name_before_ext)
                if m != None:
                # store as reverse read
                    (baseName,read) = m.groups()
                    reverse_reads[baseName] = fastq
                else:
                    print "Could not determine forward/reverse read status for input file " + fastq
        else:
            # matches default Illumina file naming format, e.g. m.groups() = ('samplename', '_S1', '_L001', '_R1', '_001')
            baseName, read  = m.groups()[0], m.groups()[3]
            if read == "_R1":
                forward_reads[baseName] = fastq
            elif read == "_R2":
                reverse_reads[baseName] = fastq
            else:
                print "Could not determine forward/reverse read status for input file " + fastq
                print "  this file appears to match the MiSeq file naming convention (samplename_S1_L001_[R1]_001), but we were expecting [R1] or [R2] to designate read as forward or reverse?"
                fileSets[file_name_before_ext] = fastq
                num_single_readsets += 1

    if args.runtype == "improvement":
        for assembly in args.assemblies:
            (file_path, file_name_before_ext, full_ext) = get_readFile_components(assembly)
            identifier = file_name_before_ext.split(args.assemblyid)[0]
            assemblies[identifier] = assembly

    # store in pairs with assemblies
    for sample in forward_reads:
        if sample in reverse_reads and sample in assemblies:
            fileSets[sample] = [forward_reads[sample],reverse_reads[sample], assemblies[sample]] # store pair and assembly
            num_paired_readsets += 1
            num_assemblies += 1
        elif sample in reverse_reads:
            fileSets[sample] = [forward_reads[sample], reverse_reads[sample]] #store pair
            num_paired_readsets += 1
        else:
            fileSets[sample] = [forward_reads[sample]] # no reverse found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + forward_reads[sample])
    for sample in reverse_reads:
        if sample not in fileSets:
            fileSets[sample] = reverse_reads[sample] # no forward found
            num_single_readsets += 1
            logging.info('Warning, could not find pair for read:' + reverse_reads[sample])
    for sample in assemblies:
        if sample not in fileSets:
            fileSets[sample] = assemblies[sample]
            num_assemblies += 1
            logging.info('Warning, could not find reads for assembly:' + assemblies[sample])

    if num_paired_readsets > 0:
        logging.info('Total paired readsets found:' + str(num_paired_readsets)) 
    if num_single_readsets > 0:
        logging.info('Total single reads found:' + str(num_single_readsets))
    if num_assemblies > 0:
        logging.info('Total number of assemblies found:' + str(num_assemblies))

    return fileSets

def get_kmer_size(read):

    cmd = "gunzip -c " + read + " | head -n 400"
    info = os.popen(cmd)
    seqs = []
    count = 1

    for line in info:
        if count % 4 == 2:
            seqs.append(line)
            count = count + 1
        else:
            count = count + 1

    lens = []
    total = 0
    for i in seqs:
        lens.append(len(i.split('\n')[0]))
    for i in lens:
        total = total + i
    total = total / 100
    sKmer = total / 3
    eKmer = total / 3 * 2

    return sKmer, eKmer

def check_blast_database(fasta):

    database_path = fasta + ".nin"

    if os.path.exists(database_path):
        logging.info('Index for {} is already built...'.format(fasta))
    else:
        logging.info('Building blast index for {}...'.format(fasta))
        print ' '.join(['makeblast db -in', fasta, '-dbtype nucl'])
        run_command(['makeblast db -in', fasta, '-dbtype nucl'])

def main():

    args = parse_args()

    #check output path has a final slash
    if args.output[-1] != "/":
        output_path = args.output + "/"
    else:
        output_path = args.output
    
    #set up logfile
    if args.log is True:
        logfile = output_path + "log.log"
    else:
        logfile = None
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')
    logging.info('program started')
    logging.info('command line: {0}'.format(' '.join(sys.argv)))

    fileSets = read_file_sets(args)
    print fileSets

    bwa_index(args.reference)

    for sample in fileSets:
        forward_read = fileSets[sample][0]
        reverse_read = fileSets[sample][1]
        try:
            assembly = fileSets[sample][2]
        except IndexError:
            pass

        output_sam = output_path + sample + '.sam'
        five_bam = output_path + sample + '_5.bam'
        three_bam = output_path + sample + '_3.bam'

        sKmer, eKmer = get_kmer_size(forward_read)

        VOdir_five = output_path + sample + "_VO_5"
        VOdir_three = output_path + sample + "_VO_3"
        five_assembly = output_path + sample + "_5_contigs.fasta"
        three_assembly = output_path + sample + "_3_contigs.fasta"

        VO_fiveout = VOdir_five + "/out/"
        VO_threeout = VOdir_three + "/out/"

        five_contigHits = output_path + sample + "_5_contigHits.txt"
        three_contigHits = output_path + sample + "_3_contigHits.txt"

        #map to IS reference
        #run_command(['bwa', 'mem', args.reference, forward_read, reverse_read, '>', output_sam])
        #print(' '.join(['bwa', 'mem', args.reference, forward_read, reverse_read, '>', output_sam]))

        #pull unmapped reads flanking IS
        run_command(['samtools view', '-Sb', '-f 36', output_sam, '>', five_bam])
        #print(' '.join(['samtools', 'view -Sb -f 36', output_sam, '>', five_bam]))
        run_command(['samtools view', '-Sb', '-f 4', '-F 40', output_sam, '>', three_bam])
        #print(' '.join(['samtools', 'view -Sb -f 4 -F 40', output_sam, '>', three_bam]))

        #assemble ends
        run_command(['./velvetshell.sh', VOdir_five, str(sKmer), str(eKmer), five_bam, VO_fiveout, five_assembly])
        #print(' '.join(['./velvetshell.sh', VOdir_five, str(sKmer), str(eKmer), five_bam, VO_fiveout, five_assembly]))
        run_command(['./velvetshell.sh', VOdir_three, str(sKmer), str(eKmer), three_bam, VO_threeout, three_assembly])
        #print(' '.join(['./velvetshell.sh', VOdir_three, str(sKmer), str(eKmer), three_bam, VO_threeout, three_assembly]))

        if args.runtype == "improvement":

            #check database for assemblies and create one if it doesn't already exist
            check_blast_database(assembly)

            #get prefix for output filenames
            genbank_output = sample + "_annotated.gbk"
            final_genbank = sample + "_annotatedAll.gbk"
            table_output = sample + "_table.txt"

            #blast ends against assemblies
            run_command(['blastn', '-db', args.assemblies, '-query', five_assembly, "-max_target_seqs 1 -outfmt '6 qseqid qlen sacc pident length slen sstart send evalue bitscore' >", five_contigHits])
            #print(' '.join(['blastn', '-db', assembly, '-query', five_assembly, "-max_target_seqs 1 -outfmt '6 qseqid qlen sacc pident length slen sstart send evalue bitscore' >", five_contigHits]))
            run_command(['blastn', '-db', args.assemblies, '-query', five_assembly, "-max_target_seqs 1 -outfmt '6 qseqid qlen sacc pident length slen sstart send evalue bitscore' >", three_contigHits])
            #print(' '.join(['blastn', '-db', assembly, '-query', five_assembly, "-max_target_seqs 1 -outfmt '6 qseqid qlen sacc pident length slen sstart send evalue bitscore' >", three_contigHits]))

            #annotate hits to genbank
            run_command(['python', 'annotateMultiGenbank.py', '-s', five_contigHits, '-f', args.aseemblies, '-p', '80', '-c', '80', '-i', sample, '-n', genbank_output ])
            #print(' '.join(['python', 'annotateMultiGenbank.py', '-s', five_contigHits, '-f', assembly, '-p', '80', '-c', '80', '-i', sample, '-n', genbank_output ]))
            run_command(['python', 'annotateMultiGenbank.py', '-s', three_contigHits, '-g', genbank_output, '-n', final_genbank, '-p', '80', '-c', '80'])
            #print(' '.join(['python', 'annotateMultiGenbank.py', '-s', three_contigHits, '-g', genbank_output, '-n', final_genbank, '-p', '80', '-c', '80']))

            #create output table
            run_command(['python', 'createTableImprovement.py', '--genbank', final_genbank, '>', table_output])
            #print(' '.join(['python', 'createTableImprovement.py', '--genbank', final_genbank, '>', table_output]))

        if args.runtype == "typing":
            pass

            #check database for reference genome and create if it doesn't exist

            #blast ends against reference genome

            #annotate hits to a genbank

            #create output table


if __name__ == '__main__':
    main()
