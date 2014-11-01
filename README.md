# ISMapper

This program takes paired end Illumina short read sequence data, an IS query of interest and a reference genome or assembly and reports the locations of the IS query in the reference genome or the assembly.

## Dependencies
* Python v2.7.5
* BioPython v1.63 - http://biopython.org/wiki/Main_Page
* BWA v0.7.5a - http://bio-bwa.sourceforge.net/
* Samtools v0.1.19 - http://samtools.sourceforge.net/
* Bedtools v2.20.1 - http://bedtools.readthedocs.org/en/latest/content/installation.html
* BLAST+ v2.2.28 - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/

## Installation
Install Python and its dependencies first.

Download ISMapper from GitHub:
 ```
git clone https://github.com/jhawkey/IS_mapper/
```

Install with pip:

```
pip install IS_mapper/
```

Testing the installation:

```
ismap --version
compiled_table.py -h
```

## Usage

There are two possible options for running ISMapper, depending on what reference genome you would like to compare to. The typing option looks for your IS query locations in your short read data, and compares these locations to a reference genome (which may or may not have that particular location in it).
Alternately, ISMapper can also be used for assembly improvement if the improvement option is chosen. Instead of a reference genome being supplied, the assembly of the short read data can be supplied instead. ISMapper then looks for your IS query locations in the assembly, indicating which contigs may be able to be joined together based on the available IS evidence.

### Typing

Input files:
* short read data in fastq format (can be gzipped)
* IS query of interest in fasta format
* reference genome to compare to in genbank format

Basic usage for ISMapper:

Multiple read sets can be supplied one after the other, separated by spaces, after --reads. ISMapper will pair the reads together.

```
ismap --reads [isolateA_1.fastq.gz] [isolateA_2.fastq.gz] [isolateB_1.fastq.gz] [isolateB_2.fastq.gz] --query IS_query.fasta --typingRef reference_genome.gbk --runtype typing --output prefix_out
```

Once ISMapper has finished running, for each isolate there will be multiple output files, the most interesting of which is the *_table.txt file, showing each location in the reference genome where there is a copy of your IS query in your isolate.
These files can be compiled together to generate one large table showing all possible IS query locations in all isolates as well as the reference genome by using the compiled_table script.

```
compiled_table.py --tables *_table.txt --reference_gbk reference_genome.gbk --seq IS_query.fasta --output compiled_table_out.txt
```

### Improvement

Input files:
* short read sequence data in fastq format (can be gzipped)
* IS query of interest in fasta format
* assemblies for each of your isolates, either in fasta or genbank format

Basic usage:

When supplying assemblies for each isolate, it is vital that they all share the same file name structure, with only the isolate name differentiating them.
Eg:
isolateA_assembly.fasta matches to isolateA_1.fastq.gz and isolateA_2.fastq.gz
isolateB_assembly.fasta matches to isolateB_1.fastq.gz and isolateB_2.fastq.gz
NOT
isolateA_assembly.fasta matches to isolateA_1.fastq.gz and isolateA_2.fastq.gz
isoalteB_contigs.fasta matches to isolateB_1.fastq.gz and isolateB_2.fastq.gz

```
ismap --reads [isolateA_1.fastq.gz] [isolateA_2.fastq.gz] [isolateB_1.fastq.gz] [isolateB_2.fastq.gz] --query IS_query.fasta --assemblies [isolateA_assembly.fasta] [isolateB_assembly.fasta] --assemblyid _assembly --extension .fasta --type fasta --runtype improvement --output prefix_out
```

Once ISMapper has finished running, multiple output files will be generated, the most interesting of which will be the *_table.txt file. The table file contains the names of contigs that have either a left or a right end in them.

## Adavnced Options for ismap

`--forward` and `--reverse` are used to designate which reads are forward and reverse only used if NOT inMiSeq format sample_S1_L001_R1_001.fastq.gz; otherwise default is _1, i.e. expect forward reads as sample_1.fastq.gz)

`--cutoff` is used to determine the read depth at a position that will be called as a 'true' peak, and thus make it to the next round of analysis. Default is currently 6 - if you're having issues detecting locations that you feel like should be there, lowering this cutoff may assist in finding them.

`--a` and `--T` are flags that are passed to BWA. --a will turn on all alignment reporting in BWA, and --T is used to give an integer mapping score to BWA to determine what alignments are kept. These options may be useful in finding IS query positions that are next to repeated elements as BWA will report all hits for the read not just the best random hit. However, using these options may cause noise and confusion in the final output files.

`--cds`, `--trna` and `--rrna` are used to specify what qualifiers will be looked for in the reference genbank when determining genes flanking the IS query location. Defaults are locus_tag and product.

`--log` turns on the log file.

`--temp` turns on keeping the temporary files instead of deleting them once the run has completed.


## Other options for compiled_table

`--gap` determines the overlap between nearby positions to be called as the same position in the final compiled tables output. Default is 0 (so positions must be exactly the same in different isolates to compile together), however increasing this number can simplify the final output

`--cds`, `--trna` and `--rrna` are used to specify what qualifiers will be looked for in the reference genbank when determining genes flanking the IS query location. Defaults are locus_tag and product.

## Running multiple jobs on a SLURM queing system

```
python slurm_ismap.py -h
usage: slurm_ismap.py [-h] [--walltime WALLTIME] [--memory MEMORY]
                      [--rundir RUNDIR] --script SCRIPT --query QUERY --reads
                      READS [READS ...] [--forward FORWARD]
                      [--reverse REVERSE]
                      [--assemblies ASSEMBLIES [ASSEMBLIES ...]]
                      [--assemblyid ASSEMBLYID] --runtype RUNTYPE
                      [--logprefix LOGPREFIX] [--other_args OTHER_ARGS]

Submit ISMapper jobs to SLURM

optional arguments:
  -h, --help            show this help message and exit
  --walltime WALLTIME   Amount of wall time. Default 1 hr
  --memory MEMORY       Amount of memory (in MB). Default is 16gb
  --rundir RUNDIR       Directory to run in. Default is current directory
  --script SCRIPT       Location of ISMapper script, ismap.py
  --query QUERY         Path to IS query.
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
```