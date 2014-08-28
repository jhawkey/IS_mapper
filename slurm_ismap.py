import sys, re, os
from argparse import (ArgumentParser, FileType)

def parse_args():

    parser = ArgumentParser(description="Submit ISMapper jobs to SLURM")

    parser.add_argument('--walltime', type=str, required=False, help='Amount of wall time. Default 1 hr', default='0-01:00:00')
    parser.add_argument('--memory', type=str, required=False, help='Amount of memory (in MB). Default is 16gb', default='16384')
    parser.add_argument('--rundir', type=str, required=False, help='Directory to run in. Default is current directory')

    parser.add_argument('--script', type=str, required=True, help='Location of ISMapper script, ismap.py')
    parser.add_argument('--reference', type=str, required=True, help='Path to IS reference.')
    parser.add_argument('--reads', nargs='+', type=str, required=True, help='Paired end read files in fastq.gz format')
    parser.add_argument('--forward', type = str, required=False, default = '_1', help = 'Identifier for forward reads if not in MiSeq format (default _1)')
    parser.add_argument('--reverse', type=str, required=False, default='_2', help='Identifier for reverse reads if not in MiSeq format (default _2)')
    parser.add_argument('--assemblies', nargs='+', type=str, required=False, help='Contig assemblies, one for each read set (If using improvement option)')
    parser.add_argument('--assemblyid', type=str, required=False, help='Identifier for assemblies eg: sampleName_contigs (specify _contigs) or sampleName_assembly (specify _assembly). Do not specify extension.')
    parser.add_argument('--runtype', type=str, required=True, help='Runtype for the program, either improvement or typing')
    parser.add_argument('--logprefix', type=str, required=False, help='Creates a prefix for the log file (default is just sample name)', default='')
    parser.add_argument('--other_args', type=str, required=False, help='String containing all other arguments to pass to ISMapper')

    return parser.parse_args()

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
            print('Warning, could not find pair for read:' + forward_reads[sample])
    for sample in reverse_reads:
        if sample not in fileSets:
            fileSets[sample] = reverse_reads[sample] # no forward found
            num_single_readsets += 1
            print('Warning, could not find pair for read:' + reverse_reads[sample])
    for sample in assemblies:
        if sample not in fileSets:
            fileSets[sample] = assemblies[sample]
            num_assemblies += 1
            print('Warning, could not find reads for assembly:' + assemblies[sample])

    if num_paired_readsets > 0:
        print('Total paired readsets found:' + str(num_paired_readsets)) 
    if num_single_readsets > 0:
        print('Total single reads found:' + str(num_single_readsets))
    if num_assemblies > 0:
        print('Total number of assemblies found:' + str(num_assemblies))

    return fileSets

def main():

    args = parse_args()

    if not args.rundir:
        args.rundir = os.getcwd()


    fileSets = read_file_sets(args)
    print fileSets

    for sample in fileSets:

        cmd = "#!/bin/bash"
        cmd += "\n#SBATCH -p main"
        cmd += "\n#SBATCH --job-name=ismapper" + sample
        cmd += "\n#SBATCH --ntasks=1"
        cmd += "\n#SBATCH --mem-per-cpu=" + args.memory
        cmd += "\n#SBATCH --time=" + args.walltime
        cmd += "\ncd " + args.rundir
        cmd += "\nmodule load python-gcc/2.7.5"
        cmd += "\nmodule load bwa-intel/0.7.5a"
        cmd += "\nmodule load samtools-gcc/0.1.19"
        cmd += "\nmodule load blast+-gcc/2.2.25"
        cmd += "\nmodule load velvetoptimiser/2.2.5"
        cmd += "\nmodule load bamtools-intel/2.3.0"
        cmd += "\nmodule load spades-gcc/3.0.0"
        cmd += "\nmodule load bedtools-intel/2.20.1"
        cmd += "\npython " + args.script
        cmd += " --reference " + args.reference
        if args.runtype == 'typing':
            cmd += " --runtype typing --reads " + fileSets[sample][0] + " " + fileSets[sample][1]
        elif args.runtype == 'improvement':
            cmd += " --runtype improvement --reads " + fileSets[sample][0] + " " + fileSets[sample][1] + " --assemblies " + fileSets[sample][2]
        if args.assemblyid:
            cmd += " --assemblyid " + args.assemblyid
        if args.logprefix == '':
            cmd += " --log --output " + sample
        elif args.logprefix != '':
            cmd += " --log --output " + args.logprefix + "_" + sample
        cmd += " " + args.other_args

        print cmd

        os.system('echo "' + cmd + '" | sbatch')

if __name__ == '__main__':
    main() 
