<p align="center"><img src="logo.png" alt="Snouter" width="500"></p>

Snouter is a tool to count the number of single nucleotide differences (a.k.a. SNPs) between very similar bacterial genome assemblies. It works with many types of assemblies: fragmented draft assemblies from Illumina reads, long read assemblies with a higher error rate, and completed assemblies.

__IMPORTANT: This tool is in active development and not yet ready for general use. Check back soon!__

[![Build Status](https://travis-ci.org/rrwick/Snouter.svg?branch=master)](https://travis-ci.org/rrwick/Snouter) [![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![Coverage Status](https://coveralls.io/repos/github/rrwick/Snouter/badge.svg?branch=master)](https://coveralls.io/github/rrwick/Snouter?branch=master)



## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [Example usage](#example-usage)
* [Full usage](#full-usage)
* [FAQ](#faq)
* [License](#license)



## Introduction

The naive way to count SNPs between two assemblies would be to align them to each other and count the differences. This would work great on flawless assemblies, but real assemblies usually have errors that can result in false positive SNPs. For example, imagine two identical 5 Mbp bacterial genomes, separately sequenced and assembled. If the assemblies have an average quality of Q50 (99.999% accurate), there may be 50 errors in each. Naively counting the SNPs could then give us a value near 100, not the 0 we should get for identical genomes. One can see that for very similar genomes (e.g. 10s of genuine SNPs), the false postives caused by assembly errors can be a serious problem. It is even more problematic for Nanopore-only assemblies which may have a average quality of Q30 (99.9% accurate), where the false positives can number in the thousands.

To avoid this problem, this tool masks out unreliable parts of genome assemblies, so it can count SNPs only in the reliable regions. It does this in two ways: 1) examining base-level confidence by mapping sequencing reads back to their own assembly, and 2) excluding repeat regions which are difficult to assemble correctly, especially using short reads. This means that even in draft assemblies and noisy assemblies, small numbers of SNP differences can be reliably counted.



## Requirements

### Input data

* Two or more assemblies in FASTA format
* For each assembly, reads in FASTQ format - _the same reads used to make the assembly_

### Software

Snouter requires some tools to be installed and available on the command line. If in doubt, open a terminal and try running the commands listed below. If they work, you should be good to go!

| Tool | Command(s) | Read types |
| ---- | -----------| ---------- |
| [MUMmer](http://mummer.sourceforge.net/) | `nucmer`, `show-snps` | all |
| [Samtools](http://www.htslib.org/) | `samtools` | all |
| [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) | `bowtie2-build`, `bowtie2` | short (Illumina) reads |
| [Minimap2](https://github.com/lh3/minimap2) | `minimap2` | long (Nanopore or PacBio) reads |



## Installation



## Example usage

Snouter is run into two stages:
1) Creating a mask file for each assembly.
    * This step is relatively slow, as it involves aligning reads and assessing each base in the assemblies. However, it only needs to be run once per assembly (_O_(n) complexity).
2) Conducting pairwise assembly alignment to get SNPs.
    * This step is much faster, but it is done on every pair of assemblies (_O_(n<sup>2</sup>) complexity) so may be the slow if you have a large number of assemblies to compare.



## Full usage



## FAQ



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
