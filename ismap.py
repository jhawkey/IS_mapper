import logging
import sys, re, os
from argparse import (ArgumentParser, FileType)
import subprocess
from subprocess import call, check_output, CalledProcessError, STDOUT
from Bio.Blast.Applications import NcbiblastnCommandline

def parse_args():
    '''
    Parse the input arguments, use -h for help.
    '''

    parser = ArgumentParser(description='IS mapper')

    # need to add verison info later
    #parser.add_argument("--version", action='version', ...)

    # Inputs
    parser.add_argument('--runtype', type=str, required=True, help='"typing" or "improvement"')
    parser.add_argument('--reads', nargs = '+', type = str, required=False, help='Paired end reads for analysing (can be gzipped)')
    parser.add_argument('--forward', type = str, required=False, default = '_1', help = 'Identifier for forward reads if not in MiSeq format (default _1)')
    parser.add_argument('--reverse', type=str, required=False, default='_2', help='Identifier for reverse reads if not in MiSeq format (default _2)')
    parser.add_argument('--reference', type = str, required=True, help='Fasta file for reference gene (eg: insertion sequence) that will be mapped to')
    parser.add_argument('--assemblies', nargs = '+', type=str, required=False, help='Contig assemblies, one for each read set')
    parser.add_argument('--assemblyid', type=str, required=False, help='Identifier for assemblies eg: sampleName_contigs (specify _contigs) or sampleName_assembly (specify _assembly). Do not specify extension.')
    parser.add_argument('--extension', type=str, required=False, help='Extension for assemblies (.fasta or .gbk, default is .fasta)', default='.fasta')
    parser.add_argument('--typingRef', type=str, required=False, help='Reference genome for typing against')
    parser.add_argument('--type', type=str, required=False, default='fasta', help='Indicator for contig assembly type, genbank or fasta (default fasta)')
    parser.add_argument('--path', type=str, required=True, default='', help='Path to folder where scripts are.')

    # Cutoffs for annotation
    parser.add_argument('--coverage', type=float, required=False, default=90.0, help='Minimum coverage for hit to be annotated (default 90.0)')
    parser.add_argument('--percentid', type=float, required=False, default=90.0, help='Minimum percent ID for hit to be annotated (default 90.0')
    
    # Reporting options
    parser.add_argument('--log', action="store_true", required=False, help='Switch on logging to file (otherwise log to stdout')
    parser.add_argument('--output', type=str, required=True, help='prefix for output files')
    parser.add_argument('--temp', action="store_true", required=False, help='Switch on keeping the temp folder instead of deleting it at the end of the program')


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
        exit_status = call(command_str, **kwargs) 
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

def check_command(command_list, command_name):
    '''
    Check that the dependency is installed.
    Exits the program if it can't be found.
        - command_list is the command to run to determine the version.
        - command_name is the name of the command to show in the error message.
    '''
    try:
        command_stdout = check_output(command_list, stderr=STDOUT)
    except OSError as e:
        logging.error("Failed command: {}".format(' '.join(command_list)))
        logging.error(str(e))
        logging.error("Do you have {} installed in your PATH?".format(command_name))
        exit(-1)
    except CalledProcessError as e:
        # some programs such as samtools return a non zero exit status
        # when you ask for the version. We ignore it here.
        command_stdout = e.output

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
        os.system(' '.join(['makeblastdb -in', fasta, '-dbtype nucl']))
        #run_command(['makeblast db -in', fasta, '-dbtype nucl'])

def make_directories(dir_list):
    '''
    Makes the directories specified in the list.
    Checks to make sure they exist first.
    If the directory is a VO directory, REMOVE IT before making it again,
    as this will cause Velvet to give an error.
    '''
    for directory in dir_list:
        if "VO" not in directory:
            run_command(['mkdir', '-p', directory], shell=True)
        elif "VO" in directory:
            if os.path.exists(directory):
                run_command(['rm', '-rf', directory], shell=True)
                run_command(['mkdir', '-p', directory], shell=True)
            else:
                run_command(['mkdir', '-p', directory], shell=True)
        else:
            logging.info('Cannot make diretory {}'.format(directory))

def main():

    args = parse_args()

    if args.path[-1] != "/":
        args.path = args.path + "/"
    
    #set up logfile
    if args.log is True:
        logfile = args.output + ".log"
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

    check_command(['bwa'], 'bwa')
    check_command(['samtools'], 'samtools')
    check_command(['VelvetOptimiser.pl', '--version'], 'VelvetOptimiser')
    check_command(['makeblastdb'], 'blast')

    # checks to make sure the runtype is valid and provides an error
    # if it's not.
    if args.runtype != "improvement" and args.runtype != "typing":
        logging.info('Invalid runtype selected: {}'.format(args.runtype))
        logging.info('Runtype should be improvement or typing (see instructions for further details)')
        exit(-1)

    fileSets = read_file_sets(args)

    bwa_index(args.reference)

    for sample in fileSets:
        forward_read = fileSets[sample][0]
        reverse_read = fileSets[sample][1]
        try:
            assembly = fileSets[sample][2]
        except IndexError:
            pass

        current_dir = os.getcwd() + '/'
        temp_folder = current_dir + sample + '_temp/'
        output_sam = temp_folder + sample + '.sam'
        five_bam = temp_folder + sample + '_5.bam'
        three_bam = temp_folder + sample + '_3.bam'
        VOdir_five = temp_folder + sample + '_VO_5'
        VOdir_three = temp_folder + sample + '_VO_3'
        VO_fiveout = VOdir_five + "/out/"
        VO_threeout = VOdir_three + "/out/"  
        five_assembly = sample + "_5_contigs.fasta"
        three_assembly = sample + "_3_contigs.fasta"
        five_contigHits = sample + "_5_contigHits.txt"
        three_contigHits = sample + "_3_contigHits.txt"

        # get Velvet kmer range
        sKmer, eKmer = get_kmer_size(forward_read)

        # create all required directories
        make_directories([VOdir_five, VOdir_three])

        # map to IS reference
        run_command(['bwa', 'mem', args.reference, forward_read, reverse_read, '>', output_sam], shell=True)

        # pull unmapped reads flanking IS
        run_command(['samtools view', '-Sb', '-f 36', output_sam, '>', five_bam], shell=True)
        run_command(['samtools view', '-Sb', '-f 4', '-F 40', output_sam, '>', three_bam], shell=True)

        # assemble ends
        #run_command(["cd", VOdir_five], shell=True) 
        #run_command(["VelvetOptimiser.pl", "-s", str(sKmer), "-e", str(eKmer), "-f '-short -bam ../" + five_bam + "'"])
        #run_command(['cd ../', '&&', 'mv', VOdir_five, '/auto*/contigs.fa', five_assembly], shell=True)
        #run_command(["cd", VOdir_three], shell=True)
        #run_command(["cd", VOdir_three, "VelvetOptimiser.pl", "-s", str(sKmer), "-e", str(eKmer), "-f '-short -bam ../" + three_bam + "'"])
        #run_command(['cd ../', '&&', 'mv', VOdir_three, '/auto*/contigs.fa', three_assembly], shell=True)
        print ' '.join([args.path + 'velvetshell.sh', VOdir_five, str(sKmer), str(eKmer), five_bam, current_dir + VO_fiveout, current_dir + five_assembly])
        print ' '.join([args.path + 'velvetshell.sh', VOdir_three, str(sKmer), str(eKmer), three_bam, current_dir + VO_threeout, current_dir + three_assembly])
        run_command([args.path + 'velvetshell.sh', VOdir_five, str(sKmer), str(eKmer), five_bam, VO_fiveout, current_dir + five_assembly], shell=True)
        run_command([args.path + 'velvetshell.sh', VOdir_three, str(sKmer), str(eKmer), three_bam, VO_threeout, current_dir + three_assembly], shell=True)

        if args.runtype == "improvement":

            # get prefix for output filenames
            genbank_output = temp_folder + sample + "_annotated.gbk"
            final_genbank = sample + "_annotatedAll.gbk"
            final_genbankSingle = sample + "_annotatedAllSingle.gbk"
            table_output = sample + "_table.txt"

            if args.extension == '.gbk':
                assembly_gbk = assembly
                (file_path, file_name_before_ext, full_ext) = get_readFile_components(assembly_gbk)
                assembly_fasta = os.path.join(temp_folder, file_name_before_ext) + '.fasta'
                run_command(['python', args.path + 'gbkToFasta', '-i', assembly, '-o', assembly_fasta], shell=True)
                assembly = assembly_fasta

            # check database for assemblies and create one if it doesn't already exist
            check_blast_database(assembly)

            # blast ends against assemblies
            run_command(['blastn', '-db', assembly, '-query', five_assembly, "-max_target_seqs 1 -outfmt '6 qseqid qlen sacc pident length slen sstart send evalue bitscore' >", five_contigHits], shell=True)
            run_command(['blastn', '-db', assembly, '-query', three_assembly, "-max_target_seqs 1 -outfmt '6 qseqid qlen sacc pident length slen sstart send evalue bitscore' >", three_contigHits], shell=True)

            # annotate hits to genbank
            run_command(['python', args.path + 'annotateMultiGenbank.py', '-s', five_contigHits, '-f', assembly, '-p', str(args.percentid), '-c', str(args.coverage), '-i', sample, '-n', genbank_output ], shell=True)
            run_command(['python', args.path + 'annotateMultiGenbank.py', '-s', three_contigHits, '-g', genbank_output, '-n', final_genbank, '-p', str(args.percentid), '-c', str(args.coverage)], shell=True)

            # create single genbank and output table
            run_command(['python', args.path + 'multiGenbankToSingle.py', '-i', final_genbank, '-n', sample, '-o', final_genbankSingle], shell=True)
            run_command(['python', args.path + 'createTableImprovement.py', '--genbank', final_genbankSingle, '--output', table_output], shell=True)

            if args.temp == False:
                run_command(['rm', '-rf', temp_folder], shell=True)

        if args.runtype == "typing":
            pass

            #check database for reference genome and create if it doesn't exist
            #check_blast_database(args.typingRef)

            #blast ends against reference genome
            #run_command(['blastn', '-db', args.typingRef, '-query', five_assembly, "-max_target_seqs 1 -outfmt '6 qseqid qlen sacc pident length slen sstart send evalue bitscore' >", five_contigHits], shell=True)
            #run_command(['blastn', '-db', args.typingRef, '-query', three_assembly, "-max_target_seqs 1 -outfmt '6 qseqid qlen sacc pident length slen sstart send evalue bitscore' >" three_contigHits], shell=True)

            #turn typingRef into a fasta
            #typingRefFasta = args.typingRef + '.fasta'
            #run_command(['python', 'gbkToFasta.py', '-i', args.typingRef, '-o', typingRefFasta])

            #annotate hits to a genbank
            #run_command(['python', 'annotateMultiGenbank.py', '-s', five_contigHits, '-f', typingRefFasta, '-p', str(args.percentid), '-c', str(args.coverage), '])

            #create output table


if __name__ == '__main__':
    main()
