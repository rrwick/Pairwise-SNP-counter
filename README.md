<p align="center"><img src="logo.png" alt="Snouter" width="500"></p>

Snouter is a tool to count the number of single nucleotide differences (a.k.a. SNPs) between very similar bacterial genome assemblies. It works with many types of assemblies: fragmented draft assemblies from Illumina reads, long read assemblies with a higher error rate, and completed assemblies.

__IMPORTANT: This tool is in active development and not yet ready for general use. Check back soon!__

[![Build Status](https://travis-ci.org/rrwick/Snouter.svg?branch=master)](https://travis-ci.org/rrwick/Snouter) [![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![Coverage Status](https://coveralls.io/repos/github/rrwick/Snouter/badge.svg?branch=master)](https://coveralls.io/github/rrwick/Snouter?branch=master)



## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [For the impatient](#for-the-impatient)
* [Usage and examples](#usage-and-examples)
* [Full command line help](#full-command-line-help)
* [FAQ](#faq)
* [License](#license)



## Introduction

The naive way to count SNPs between two assemblies would be to align them to each other and count the differences. This works with flawless assemblies, but real assemblies usually have errors that result in false positive SNPs. For example, imagine two identical 5 Mbp bacterial genomes, separately sequenced and assembled. If the assemblies have an average quality of Q50 (99.999% accurate), there may be 50 errors in each. Simply counting the differences could then give us a SNP count near 100, not the zero we expect for identical genomes. One can see that for very similar genomes (e.g. 10s of genuine SNPs), false postives caused by assembly errors can be a serious problem. It is even more problematic for Nanopore-only assemblies which may have a average quality of Q30 (99.9% accurate), where false positives can number in the thousands.

To avoid this problem, Snouter masks out unreliable parts of assemblies so it can count SNPs only in reliable regions. It does this in two ways: 1) examining base-level confidence by mapping sequencing reads back to their own assembly, and 2) excluding repeat regions which are difficult to assemble correctly (especially using short reads). This means that even in draft assemblies and noisy long read assemblies, small numbers of SNP differences can be reliably counted.



## Requirements

### Input data

* Two or more assemblies in FASTA format
* For each assembly, reads in FASTA/FASTQ format - _the same reads used to make the assembly_

### Software

Snouter requires some tools to be installed and available on the command line. If in doubt, open a terminal and try running the commands listed below. If they work, you should be good to go!

| Tool | Command(s) | Read types |
| ---- | -----------| ---------- |
| [MUMmer](http://mummer.sourceforge.net/) | `nucmer`, `show-snps` | all |
| [Samtools](http://www.htslib.org/) | `samtools` | all |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) | `bowtie2-build`, `bowtie2` | short (Illumina) reads |
| [Minimap2](https://github.com/lh3/minimap2) | `minimap2` | long (Nanopore or PacBio) reads |



## Installation



## For the impatient

```
snouter mask --assembly_fp sample_1_assembly.fasta --read_fps sample_1_reads_1.fastq.gz sample_1_reads_2.fastq.gz --read_type illumina
snouter mask --assembly_fp sample_2_assembly.fasta --read_fps sample_2_reads_1.fastq.gz sample_2_reads_2.fastq.gz --read_type illumina
snouter mask --assembly_fp sample_3_assembly.fasta --read_fps sample_3_reads_1.fastq.gz sample_3_reads_2.fastq.gz --read_type illumina
snouter count --assembly_fps sample_1_assembly.fasta sample_2_assembly.fasta sample_3_assembly.fasta
```


## Usage and examples

### Snouter mask

The first stage in Snouter, `snouter mask`, creates a mask file for each assembly: a file indicating which bases are not trustworthy. This step is relatively slow, as it involves aligning reads and assessing each base in the assemblies. However, it only needs to be run once per assembly (_O_(n) complexity).

Examples:
* Masking a paired-read Illumina assembly:<br>
`snouter mask --assembly_fp assembly.fasta --read_fps reads_1.fastq.gz reads_2.fastq.gz --read_type illumina`
* Masking a long read assembly:<br>
`snouter mask --assembly_fp assembly.fasta --read_fps nanopore_reads.fastq.gz --read_type long`
* Bulk masking using an input file (specifications below):<br>
`snouter mask --bulk_mask bulk.tsv`

Bulk masking:
* This file takes the format of a tab separated values file, with the following fields as described below:<br>
`*EXCLUDE*  READ_TYPE   ASSEMBLY_FP ASSEMBLY_FPS_1 *ASSEMBLY_FPS_2* *ASSEMBLY_FPS_3*`
* Optional fields are marker with asterix. An example of bulk.tsv  is as follows:<br>
```
illumina	as1.fasta	r1_1.fastq.gz	r1_2.fastq.gz
long	as2.fasta	r2.fastq.gz
7	illumina	as3.fasta	r3_1.fastq.gz	r3_2.fastq.gz
```

### Snouter count

The second stage is `snouter count`, which conducts pairwise assembly alignment to get SNP counts. This step is much faster, but it is done on every pair of assemblies (_O_(n<sup>2</sup>) complexity) so may be the slow if you have a large number of assemblies to compare.

Examples:
* example_1
* example_2
* example_3
* example_4


## Full command line help

### Snouter mask

```
usage: snouter mask [-h] --assembly_fp ASSEMBLY_FP --read_fps READ_FPS
                    [READ_FPS ...] --read_type {illumina,long}
                    [--threads THREADS] [--exclude EXCLUDE]

optional arguments:
  -h, --help            show this help message and exit
  --assembly_fp ASSEMBLY_FP
                        Input assembly filepath
  --read_fps READ_FPS [READ_FPS ...]
                        Input read filepaths, space separated
  --read_type {illumina,long}
                        Read type of input reads. [choices: illumina, long]
  --threads THREADS     Number of threads
  --exclude EXCLUDE     Percentage of assembly bases to exclude
```


### Snouter count

```
usage: snouter count [-h] --assembly_fps ASSEMBLY_FPS [ASSEMBLY_FPS ...]
                     [--mask_fps MASK_FPS [MASK_FPS ...]]

optional arguments:
  -h, --help            show this help message and exit
  --assembly_fps ASSEMBLY_FPS [ASSEMBLY_FPS ...]
                        Input assembly filepaths, space separated
  --mask_fps MASK_FPS [MASK_FPS ...]
                        Input masking filepaths, space separated
```



## FAQ



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
