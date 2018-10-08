# ISMapper

This program takes paired end Illumina short read sequence data, an IS query/queries of interest and a reference genome and reports the locations of the IS query in the reference genome or the assembly.
For a more in depth description of how the program works, see **INSERT REFERENCE HERE** below.

For support, please create an issue in the GitHub issue tracker.

Read more about ISMapper here:
> [__Hawkey J, Hamidian M, Wick RR, Edwards DJ, Billman-Jacobe H, Hall RM, Holt KE.__ ISMapper: identifying transposase insertion sites in bacterial genomes from short read sequence data _BMC Genomics_ 2015.](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1860-2)

# Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
    * [Running a test case](#running-a-test-case)
* [Usage](#usage)
    * [Basic usage](#basic-usage)
    * [All ISMapper options](#all-ismapper-options)
    * [A note on version 2](#a-note-on-version-2)
* [Outputs](#outputs)
    * [Single isolate output](#single-isolate-output)
    * [Compiled output](#compiled-output)
        * [compiled_table.py usage](#compiled_tablepy-usage)
        * [compiled_table.py outputs](#compiled_tablepy-outputs)
* [Method](#method)
    
# Requirements
* [Python](https://www.python.org/) v2.7+ or v3.6+
* [BioPython](https://biopython.org/) v1.63 or later
* [BWA](http://bio-bwa.sourceforge.net/) v0.7.5a or later
* [SAMtools](http://www.htslib.org/) v0.1.19 or later
* [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html) v2.20.1 or later
* [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/) v2.2.28 or later

# Installation
Install Python and its dependencies first.

Download ISMapper from GitHub:
 ```
git clone https://github.com/jhawkey/IS_mapper/
```

Install with pip:

```
pip install IS_mapper/
```

Alternatively, you can install ISMapper directly from github with pip.

```
pip install git+https://github.com/jhawkey/IS_mapper
```

To install locally:

```
pip install --user git+https://github.com/jhawkey/IS_mapper
```

To install a specific release:

```
pip install git+https://github.com/jhawkey/IS_mapper@my_tag
```

Test ISMapper is installed:

```
ismap --version
compiled_table.py -h
```

## Running a test case

To check that all dependencies are working correctly and outputs are as expected:

Download _Streptococcus suis_ P1/7 reads from SRA

```
cd IS_mapper/test/inputs
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225612/ERR225612_*.fastq.gz
```

Run ISMapper, using the reads you've just downloaded, as well as the provided IS query and reference genbank in test/inputs

```
ismap --reads IS_mapper/test/inputs/ERR225612_*.fastq.gz --queries IS_mapper/test/inputs/ISSsu3.fasta --reference IS_mapper/test/inputs/S_suis_P17.gbk --output_dir IS_mapper/test/
```

The resulting table files should match those found in `test/example_results`.

# Usage
ISMapper requires three main inputs:
* paired end reads
* a reference genome
* an IS query

## Basic usage

`--reads` Multiple read sets can be supplied using a wildcard (eg: `/path/to/reads/*.fastq.gz`). They can also be supplied one after the other, separated by spaces. Regardless of how the reads are supplied, ISMapper will pair the reads belonging to the same isolate together.  

`--queries` Multiple IS queries can also be supplied, either by spaces, or as a multifasta file. ISMapper will run queries sequentially in the same output folder.  

`--reference` Multiple reference genomes can also be supplied, again either separated by spaces, or as a multigenbank file. 

**Note:** _All_ reference genomes **must** have a locus_tag for each CDS, tRNA or rRNA feature.

`--output_dir` Path to the directory where you would like your output files to go. A subdirectory for each isolate will be made within this folder. Within each isolate subdirectory, subdirectories will be made for each IS query supplied.

`--log` Prefix for the log file. If not supplied, the prefix will be current date and time.

`--help_all` Display the extended help and advanced options.

Example command:
```
ismap --reads [isolateA_1.fastq.gz] [isolateA_2.fastq.gz] [isolateB_1.fastq.gz] [isolateB_2.fastq.gz] --queries IS_query.fasta --reference reference_genome.gbk --output_dir
```

## All ISMapper options
```
Basic ISMapper options:
  --reads READS [READS ...]
                        Paired end reads for analysing (can be gzipped)
  --queries QUERIES [QUERIES ...]
                        Multifasta file for query gene(s) (eg: insertion
                        sequence) that will be mapped to.
  --reference REFERENCE [REFERENCE ...]
                        Reference genome for typing against in genbank format
  --output_dir OUTPUT_DIR
                        Location for all output files (default is current
                        directory).
  --log LOG             Prefix for log file. If not supplied, prefix will be
                        current date and time.
  --help_all HELP_ALL   Display extended help

Parameters for defining hits:
  --min_clip MIN_CLIP   Minimum size for softclipped region to be extracted
                        from initial mapping (default 10).
  --max_clip MAX_CLIP   Maximum size for softclipped regions to be included
                        (default 30).
  --cutoff CUTOFF       Minimum depth for mapped region to be kept in bed file
                        (default 6)
  --novel_gap_size NOVEL_GAP_SIZE
                        Distance in base pairs between left and right flanks
                        to be called a novel hit (default 15)
  --min_range MIN_RANGE
                        Minimum percent size of the gap to be called a known
                        hit (default 0.9, or 90 percent)
  --max_range MAX_RANGE
                        Maximum percent size of the gap to be called a known
                        hit (default 1.1, or 110 percent)
  --merging MERGING     Value for merging left and right hits in bed files
                        together to simply calculation of closest and
                        intersecting regions (default 100).

BWA parameters:
  --a                   Switch on all alignment reporting for bwa.
  --T T                 Mapping quality score for bwa (default 30).
  --t T                 Number of threads for bwa (default 1).
  --forward FORWARD     Identifier for forward reads if ISMapper is unable to
                        pair (default is Miseq format _1)
  --reverse REVERSE     Identifier for forward reads if ISMapper is unable to
                        pair (default is Miseq format _2)

Parameters for output table:
  --cds CDS             qualifier containing gene information (default
                        product). Also note that all CDS features MUST have a
                        locus_tag
  --trna TRNA           qualifier containing gene information (default
                        product). Also note that all tRNA features MUST have a
                        locus_tag
  --rrna RRNA           qualifier containing gene information (default
                        product). Also note that all rRNA features MUST have a
                        locus_tag

Reporting parameters:
  --temp                Switch on keeping the temp folder instead of deleting
                        it at the end of the run
  --bam                 Switch on keeping the final bam files instead of
                        deleting them at the end of the run
```

### A note on version 2
ISMapper version 2 is now exclusively for finding IS insertion sites relative to a reference genome. This version of ISMapper is unable to operate in 'assembly improvement mode', as was documented for version 1. If you need this functionality, please download the latest version 1 branch.

# Outputs
## Single isolate output
For each isolate, a folder is made for each IS query given. Within the folder for the IS query are the output files pertaining to that query. There will be the same set of output files for each reference sequence provided.

The main output file is `[isolate]__[reference]_table.txt`. This file contains one line per IS position found. A description of the column names is below:

* `region`: a numeric code given to the IS position
* `orientation`: the orientation of the IS position, either `F` (forward) or `R` (reverse)
* `x`: the left coordinate of the position (always the smallest value, regardless of orientation)
* `y`: the right coordinate of the position (always the largest value, regardless of orientation) 
* `gap`: the distance between the `x` and `y` coordinates. A negative value indicates that the flanking regions overlap. The sequence of this overlap is usually the target site duplication the IS has made when it inserted.
* `call`: determination of whether the position hit is a novel one (not in the reference) or known (in the reference). A `*` indicates that the gap size is larger than expected, so this hit is considered imprecise. A `?` indicates that either the left or right flank was below the depth cutoff, and so this hit is uncertain.
* `percent_ID`: if the IS hit is a known one, this is the percentage nucleotide identity to the IS query
* `percent_cov`: if the IS hit is a known one, this is the percentage coverage of the region to the IS query
* `left_gene`: locus tag ID of the gene closest to the left side of the IS position
* `left_description`: description of the gene feature from the left side of the IS position. By default this is the `product` of the gene listed in the reference GenBank.
* `left_strand`: DNA strand (1: 5' -> 3', -1: 3' -> 5') of the gene on the left side of the IS position
* `left_distance`: distance of the IS hit from the start codon of the gene on the left side of the IS position
* `right_gene`: locus tag ID of the gene closest to the right side of the IS position
* `right_description`: description of the gene feature from the right side of the IS position. By default this is the `product` position
* `right_strand`: DNA strand (1: 5' -> 3', -1: 3' -> 5') of the gene on the right side of the IS hit
* `right_distance`: distance of the IS hit from the start codon of the gene on the right side of the IS position
* `gene_interruption`: either `True` or `False` depending on whether the IS position has interrupted a gene

## Compiled output
Once ISMapper has completed, the individual output tables for each isolate can be combined to create a single table containing all IS positions across all isolates using `compiled_table.py`.

### compiled_table.py usage
Required parameters are:
`--tables` Each individual isolate table from a **single** IS query and **single** reference genbank.
`--reference` The reference genbank used to generate the individual tables. This must be a **single entry** genbank file.
`--query` The IS query sequence in fasta format.

### compiled_table.py outputs
This script generates two output files. 

The `_full_compiled.txt` file contains isolates as rows, and merged IS positions as columns. The first isolate row are the positions found in the reference genome. Values in the cells are either `-` for absent, `+` for present, `*` for present but an imprecise hit and `?` for present but an unconfident hit. The final nine rows in the table show orientation of the IS position, and surrounding gene information for each IS hit.

The `_binary_compiled.txt` contains the same information as the `_full_compiled.txt` table, but the presence or absence of the IS positions are encoded as either `0` or `1`. In this table, the additional rows showing orientation of the IS position and surrounding gene information is absent. By default, IS positions found to be imprecise (`*`) are encoded as present `1` and unconfident positions (`?`) are encoded as absent `0`. However, these values can be changed using the `--imprecise` and `--unconfident` flags when running `compiled_table.py`.

All options for `compiled_tabe.py` are show below.

```Create a table of IS hits in all isolates for ISMapper

optional arguments:
  -h, --help            show this help message and exit
  --tables TABLES [TABLES ...]
                        tables to compile
  --reference REFERENCE
                        gbk file of reference
  --query QUERY         fasta file for insertion sequence query for
                        compilation
  --gap GAP             distance between regions to call overlapping, default
                        is 0
  --cds CDS             qualifier containing gene information (default
                        product). Also note that all CDS features MUST have a
                        locus_tag
  --trna TRNA           qualifier containing gene information (default
                        product). Also note that all tRNA features MUST have a
                        locus_tag
  --rrna RRNA           qualifier containing gene information (default
                        product). Also note that all rRNA features MUST have a
                        locus_tag
  --imprecise IMPRECISE
                        Binary value for imprecise (*) hit (can be 1, 0 or
                        0.5), default is 1
  --unconfident UNCONFIDENT
                        Binary value for questionable (?) hit (can be 1, 0 or
                        0.5), default is 0
  --out_prefix OUT_PREFIX
                        Prefix for output file
```

# Method

ISMapper finds locations of an IS query in short read data using a series of mapping steps.

First, all reads for an isolate are mapped to the IS query using BWA. The reads we are interested in are the ones whose pairs have mapped, but are unmapped themselves (so therefore must be flanking the IS query). These reads are selected using SamTools.

These unmapped pairs are selected and placed into two groups - left end (flanking the left end of the IS query) and right end (flanking the right end of the IS query).
Each of these groups are independently mapped (with BWA) to the reference genome, which may or may not have a particular IS query location present.

From this mapping, Bedtools is used to find the depth at each position in the reference genome. We are looking for large peaks of left end and right end reads that indicate a possibly IS query location. Positions where the read depth falls below 6 (default, can be changed using the `cutoff` parameter) are eliminated from further analysis, as these may be incorrect alignments.

<p align="center"><img src="misc/image.png" alt="caption here"></p>

Overlapping, or close, flanking regions are merged using the Bedtools merge function to prevent multiple hits representing the same region in the final output file. By default, close regions are defined as being <= 100 bp apart.

Once left and right flanking regions have been identified, Bedtools is used to find pairs of left and right flanks which either:  
a) intersect or are <= 15 bp (default) apart, which indicates a novel position that is not present in the reference, or;   
b) separated by a distance which is approximately equal to the size of the IS query, which indicates a position that is probably in the reference.

<p align="center"><img src="misc/image.png" alt="caption here"></p>

All potential IS positions are then assessed to determine if they are tabulated into the final results file. Any hits which are removed are placed in the `_removedHits.txt` file, with the reason for their removal.