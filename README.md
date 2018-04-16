# Pairwise SNP Counter
[![Build Status](https://travis-ci.org/rrwick/Pairwise-SNP-counter.svg?branch=master)](https://travis-ci.org/rrwick/Pairwise-SNP-counter) [![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![Coverage Status](https://coveralls.io/repos/github/rrwick/Pairwise-SNP-counter/badge.svg?branch=master)](https://coveralls.io/github/rrwick/Pairwise-SNP-counter?branch=master)

__IMPORTANT NOTE: This tool is in active development and not yet ready for general use. Check back soon!__

This is a tool to count the number of single nucleotide differences (a.k.a. SNPs) between very similar bacterial genome assemblies. It works with many types of assemblies: fragmented draft assemblies from Illumina reads, long read assemblies with a higher error rate, and completed assemblies.



## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [Example commands](#example-commands)
* [FAQ](#faq)
* [License](#license)



## Introduction

The challenge in counting SNP differences comes from the fact that almost all genome assemblies have errors, leading to false positive SNPs. For example, imagine two identical 5 Mbp bacterial genomes, separately sequenced and assembled. If the assemblies have a quality of Q50 (99.999% accurate), there may be 50 errors in each. Naively counting the SNPs would then give us a value of 100, not the 0 we should get for identical genomes. One can see that for very similar genomes (e.g. 10s of genuine SNPs), the false postives caused by assembly errors can be a serious problem. It is even more problematic for Nanopore-only assemblies (Q30 to Q35), where the false positives number in the thousands.

This tool avoids that problem by masking out unreliable parts of genome assemblies, so it can count SNPs only in the reliable regions. It does this in two ways: 1) examining base-level confidence by mapping sequencing reads back to their own assembly, and 2) excluding repeat regions which are difficult to assemble correctly, especially using short reads. This means that even in draft assemblies and noisy assemblies, small numbers of SNP differences can be reliably counted.



## Requirements



## Installation



## Example commands



## FAQ



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
