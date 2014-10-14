So there are some options that don't make sense and are a little redundant, I'm in the process of cleaning those up

To use ISMapper:

For multiple read sets use the slurm_ismap.py script

optional arguments:
  -h, --help            show this help message and exit
  --walltime WALLTIME   Amount of wall time. Default 1 hr
  --memory MEMORY       Amount of memory (in MB). Default is 16gb
  --rundir RUNDIR       Directory to run in. Default is current directory
  --script SCRIPT       Location of ISMapper script, ismap.py
  --reads READS [READS ...]
                        Paired end read files in fastq.gz format
  --forward FORWARD     Identifier for forward reads if not in MiSeq format
                        (default _1)
  --reverse REVERSE     Identifier for reverse reads if not in MiSeq format
                        (default _2)
  --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        Contig assemblies, one for each read set (If using
                        improvement option)
  --assemblyid ASSEMBLYID
                        Identifier for assemblies eg: sampleName_contigs
                        (specify _contigs) or sampleName_assembly (specify
                        _assembly). Do not specify extension.
  --runtype RUNTYPE     Runtype for the program, either improvement or typing
  --logprefix LOGPREFIX
                        Creates a prefix for the log file (default is just
                        sample name)
  --other_args OTHER_ARGS
                        String containing all other arguments to pass to
                        ISMapper


--assemblies and --assemblyid options are only required if running the improvement pathway
Assemblies must be named exactly the same as the reads: eg - CP001921_1.fastq.gz, CP001921_2.fastq.gz, CP001921_assembly.fasta
The assemblyid flag is tell the program what the extension is after the sample name. In the case above, --assemblyid would be set to _assembly
DON'T put the .fasta in the --assemblyid flag, that's specified elsewhere. 

--logprefix is the prefix for the log file (I usually use the date)

--other_args should pass all the other arguments required for ismap.py

Arguments for ismap.py:

  --runtype RUNTYPE     "typing" or "improvement"
  --reads READS [READS ...]
                        Paired end reads for analysing (can be gzipped)
  --forward FORWARD     Identifier for forward reads if not in MiSeq format
                        (default _1)
  --reverse REVERSE     Identifier for reverse reads if not in MiSeq format
                        (default _2)
  --reference REFERENCE
                        Fasta file for reference gene (eg: insertion sequence)
                        that will be mapped to
  --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        Contig assemblies, one for each read set
  --assemblyid ASSEMBLYID
                        Identifier for assemblies eg: sampleName_contigs
                        (specify _contigs) or sampleName_assembly (specify
                        _assembly). Do not specify extension.
  --extension EXTENSION
                        Extension for assemblies (.fasta or .gbk, default is
                        .fasta)
  --typingRef TYPINGREF
                        Reference genome for typing against
  --type TYPE           Indicator for contig assembly type, genbank or fasta
                        (default fasta)
  --path PATH           Path to folder where scripts are.
  --coverage COVERAGE   Minimum coverage for hit to be annotated (default
                        90.0)
  --percentid PERCENTID
                        Minimum percent ID for hit to be annotated (default
                        90.0
  --log                 Switch on logging to file (otherwise log to stdout
  --output OUTPUT       prefix for output files
  --temp                Switch on keeping the temp folder instead of deleting
                        it at the end of the program

Obviously you don't need to pass reads, forward, reverse, assemblies or assemblyid.

You will need to pass:
--reference
This is the fasta file for the IS of interest
--extension (if doing improvement)
This is the .fasta or .gbk or whatever format your contigs are in
--type (if doing improvement)
This is the same as extension, not sure why it exists, don't use it!
--typingRef (only for typing)
GENBANK reference for the genome you're typing against
--path
This is the path to the folder where the scripts are. It's needed at the moment because the program isn't installed on barcoo on its own
--log
Put this if you want to switch on the log file
--output 
Prefix for output files - don't worry about setting this if you're using the slurm option
--temp
Put this if you want to keep the temp folders once the run is completed, otherwise it'll be deleted

Eg (multiple read sets):
python /vlsci/VR0082/shared/jane/IS_mapper/slurm_ismap.py --script /vlsci/VR0082/shared/jane/IS_mapper/ismap.py --reads /vlsci/VR0082/shared/acineto/fastq/vietnam2012/*fastq.gz  --forward 1 --reverse 2 --runtype typing --logprefix 290414 --other_args "--typingRef /vlsci/VR0082/shared/acineto/refs/chromosomes/CP001921.gbk --reference /vlsci/VR0082/shared/jane/IS_refs/IS26.fasta --log --path /vlsci/VR0082/shared/jane/IS_mapper/"

Eg (single read set):
python /vlsci/VR0082/shared/jane/IS_mapper/ismap.py --runtype typing --reads /vlsci/VR0082/shared/acineto/fastq/vietnam2012/8346_4#96_1.fastq.gz /vlsci/VR0082/shared/acineto/fastq/vietnam2012/8346_4#96_2.fastq.gz --forward 1 --reverse 2 --log --output 290414__8346_4#96_ --typingRef /vlsci/VR0082/shared/acineto/refs/chromosomes/CP001921.gbk --reference /vlsci/VR0082/shared/jane/IS_refs/ISAba1.fasta --log --path /vlsci/VR0082/shared/jane/IS_mapper/