#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com), Stephen Watts, Alex Tokolyi
https://github.com/rrwick/Snouter

Snouter is a command line tool for counting the single-nucleotide differences (a.k.a. SNPs) between
two closely related assemblies.

TODO: add some example commands here.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program.  If not, see
<http://www.gnu.org/licenses/>.
"""

import argparse
from distutils.version import LooseVersion
import gzip
import logging
import math
import multiprocessing
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile


def get_arguments():
    parser = argparse.ArgumentParser()

    subparser = parser.add_subparsers(metavar='mask count', dest='command')

    parser_mask = subparser.add_parser('mask')
    parser_mask.add_argument('--assembly_fp', required=True, type=pathlib.Path,
                             help='Input assembly filepath')
    # TODO: --read_fps help info stating expected read order
    parser_mask.add_argument('--read_fps', required=True, nargs='+', type=pathlib.Path,
                             help='Input read filepaths, space separated')
    parser_mask.add_argument('--read_type', required=True, choices=['illumina', 'long'],
                             help='Read type of input reads. [choices: illumina, long]')
    parser_mask.add_argument('--threads', required=False, type=int, default=default_thread_count(),
                             help='Number of threads')
    parser_mask.add_argument('--exclude', required=False, type=float,
                             help='Percentage of assembly bases to exclude')

    parser_count = subparser.add_parser('count')
    parser_count.add_argument('--assembly_fps', required=True, nargs='+', type=pathlib.Path,
                              help='Input assembly filepaths, space separated')
    parser_count.add_argument('--mask_fps', nargs='+', type=pathlib.Path,
                              help='Input masking filepaths, space separated')
    # TODO: option for output format: simple counts (one pairwise count per line), count matrix or
    #       PHYLIP distance matrix (for tree-building).

    parser.add_argument('--tmp_dir', required=False, type=pathlib.Path,
                        help='If desired, input a directory to use as temporary')

    args = parser.parse_args()
    if not args.command:
        # TODO: print better help info. see samtools for an example with subcommands
        parser.print_help()
        print('\n', end='')
        parser.error('command options include mask or count')
    return args


def main():
    # Get commandline arguments and initialise
    args = get_arguments()
    initialise_logging()
    check_arguments(args)
    check_dependencies()

    # Initialise temporary directory
    if args.tmp_dir:
        if args.tmp_dir.exists():
            tmp_dir = tempfile.TemporaryDirectory(dir=args.tmp_dir)
        else:
            logging.critical(f'Specified temporary directory ({args.tmp_dir}) does not exist')
            sys.exit(1)
    else:
        tmp_dir = tempfile.TemporaryDirectory()

    # Execute requested stage
    if args.command == 'mask':
        run_mask(args, tmp_dir.name)
    elif args.command == 'count':
        run_count(args, tmp_dir.name)


def check_arguments(args):
    log_newline()
    logging.info('Checking arguments')

    # TODO: Perform additional argument parsing, checking
    if args.command == 'mask':
        check_parsed_file_exists(args.assembly_fp)
        for read_fp in args.read_fps:
            check_parsed_file_exists(read_fp)

        if args.read_type == 'illumina':
            if len(args.read_fps) > 3:
                logging.critical('--read_fps takes no more than three illumina read sets')
                sys.exit(1)
        elif args.read_type == 'long':
            if len(args.read_fps) > 1:
                logging.critical('--read_fps takes only a single long read set')
                sys.exit(1)
        if args.exclude is None:
            args.exclude = 2.0 if args.read_type == 'illumina' else 5.0
            logging.debug(f'--exclude set to {args.exclude} based on read type of '
                          f'"{args.read_type}"')

    if args.command == 'count':
        if len(args.assembly_fps) < 2:
            logging.critical('Two or more assemblies are required, got {len(args.assembly_fps}')
            sys.exit(1)
        for assembly_fp in args.assembly_fps:
            check_parsed_file_exists(assembly_fp)
        if not args.mask_fps:
            args.mask_fps = [pathlib.Path(f'{fp}.mask') for fp in args.assembly_fps]
        if len(args.assembly_fps) != len(args.mask_fps):
            logging.critical('Need the same number of masks as assemblies')
            sys.exit(1)
        for mask_fp in args.mask_fps:
            check_parsed_file_exists(mask_fp)


def check_parsed_file_exists(filepath):
    # Check that the argument has been set; is not None
    if filepath and not filepath.exists():
        logging.critical(f'Input file {filepath} does not exist')
        sys.exit(1)


def run_mask(args, dh):
    read_filetype = check_input_mask_files(args)
    # Map reads to assembly
    if args.read_type == 'illumina':
        index_fp = index_assembly(args.assembly_fp, dh)
        bam_fp = map_illumina_reads(index_fp, args.read_fps, dh, args.threads, read_filetype)
    else:
        assert args.read_type == 'long'
        bam_fp = map_long_reads(args.assembly_fp, args.read_fps, dh, args.threads)
    scores, total_size = get_base_scores_from_mpileup(args.assembly_fp, bam_fp)
    min_score_threshold = get_score_threshold(scores, args.exclude)
    logging.debug(f'Minimum score threshold set to {min_score_threshold}')
    write_mask_file(scores, min_score_threshold, args.assembly_fp, total_size)


def run_count(args, dh):
    check_input_count_files(args)
    counts = {}
    for i in range(len(args.assembly_fps)):
        assembly_1 = args.assembly_fps[i]
        mask_1 = load_mask_file(args.mask_fps[i])
        for j in range(i+1, len(args.assembly_fps)):
            assembly_2 = args.assembly_fps[j]
            mask_2 = load_mask_file(args.mask_fps[j])
            snp_count = get_pairwise_snp_count(assembly_1, mask_1, assembly_2, mask_2, dh)
            counts[(assembly_1, assembly_2)] = snp_count
            counts[(assembly_2, assembly_1)] = snp_count
    # TODO: save results to file, in user's preferred format


class ColourFormatter(logging.Formatter):
    def format(self, record):
        s = super().format(record)
        if record.levelname == 'DEBUG':
            return '\033[2m' + s + '\033[0m'  # dim
        elif record.levelname == 'CRITICAL':
            return '\033[31m' + s + '\033[0m'  # red
        else:
            return s


def initialise_logging():
    # TODO: command line arguments for logging; save to file? print to stdout?

    fmt_str = '%(asctime)s %(message)s'
    datefmt_str = '%d/%m/%Y %H:%M:%S'
    logger = logging.getLogger()

    log_filehandler = logging.FileHandler('run.log', mode='w')
    file_formatter = logging.Formatter(fmt=fmt_str, datefmt=datefmt_str)
    log_filehandler.setFormatter(file_formatter)
    logger.addHandler(log_filehandler)

    log_streamhandler = logging.StreamHandler()
    stdout_formatter = ColourFormatter(fmt=fmt_str, datefmt=datefmt_str)
    log_streamhandler.setFormatter(stdout_formatter)
    logger.addHandler(log_streamhandler)

    # TODO: expose this option to the command line
    logger.setLevel(logging.DEBUG)


def log_newline():
    log_no_format('')


def log_no_format(msg):
    saved_formatters = [x.formatter for x in logging.getLogger().handlers]
    for h in logging.getLogger().handlers:
        h.setFormatter(logging.Formatter(fmt=''))
    logging.info(msg)
    for i, formatter in enumerate(saved_formatters):
        logging.getLogger().handlers[i].formatter = formatter


def check_dependencies():
    log_newline()
    logging.info('Checking program dependencies')
    dependencies = {'samtools': {'vcommand': 'samtools 2>&1',
                                 'vregex': re.compile(r'Version: ([^ \n]+)'),
                                 'vrequired': '1.0'},
                    'bowtie2':  {'vcommand': 'bowtie2 --version',
                                 'vregex': re.compile(r'^.+version (.+)'),
                                 'vrequired': '2.0'},
                    'minimap2': {'vcommand': 'minimap2 --version',
                                 'vregex': re.compile(r'^(.+)'),
                                 'vrequired': '2.0'},
                    'bash':     {'vcommand': 'bash --version',
                                 'vregex': re.compile(r'^.+version (.+) .+'),
                                 'vrequired': '1.0'},
                    'nucmer':   {'vcommand': 'nucmer --version 2>&1',
                                 'vregex': re.compile(r'version (.+)'),
                                 'vrequired': '3.0'}}
    for dependency, version_data in dependencies.items():
        if not shutil.which(dependency):
            logging.critical(f'Could not find dependency {dependency}')
            sys.exit(1)
        result = execute_command(version_data['vcommand'], check=False)
        try:
            version = version_data['vregex'].search(result.stdout).group(1)
        except AttributeError:
            # TODO: should we have an option to skip dependency check?
            logging.critical(f'Unable to determine version for {dependency}')
            sys.exit(1)
        if LooseVersion(version) < LooseVersion(version_data['vrequired']):
            logging.critical(f'{dependency} version {version_data["vrequired"]} or high is '
                             f'required')
            sys.exit(1)
        else:
            logging.debug(f'{dependency}: version {version} (good)')


def get_sequence_filetype(filename):
    with get_open_function(filename)(filename, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                return 'FASTA'
            elif line.startswith('@'):
                return 'FASTQ'
            else:
                logging.critical(f'Could not determine type of {filename} (must be FASTA or FASTQ')
                sys.exit(1)


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    with open(filename, 'rb') as unknown_file:
        file_start = unknown_file.read(max_len)
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        logging.critical('cannot use bzip2 format - use gzip instead')
        sys.exit(1)
    if compression_type == 'zip':
        logging.critical('cannot use zip format - use gzip instead')
        sys.exit(1)
    return compression_type


def get_open_function(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def check_input_mask_files(args):
    log_newline()
    logging.info('Checking input file types')

    assembly_filetype = get_sequence_filetype(args.assembly_fp)
    logging.debug(f'{args.assembly_fp} ({assembly_filetype})')

    read_filetypes = set()
    for read_fp in args.read_fps:
        read_filetype = get_sequence_filetype(read_fp)
        logging.debug(f'{read_fp} ({read_filetype})')
        read_filetypes.add(read_filetype)
    if args.read_type == 'illumina' and len(read_filetypes) > 1:
        logging.critical('Read files must be all FASTQ or all FASTA (not a mixture of both)')
        sys.exit(1)
    assert len(read_filetypes) == 1

    read_names = list()
    if len(args.read_fps) > 1:
        for read_fp in args.read_fps:
            with get_open_function(read_fp)(read_fp, 'rt') as fh:
                read_names.append(fh.readline())
        for forward_char, reverse_char in zip(*read_names[:2]):
            if forward_char != reverse_char:
                break
        if forward_char != '1' or reverse_char != '2':
            msg = 'Read sets specified in --read_fps do not appear to be paired.'
            msg += ' Reads are expected to be passed as --read_fps <forward> <reverse> or '
            msg += ' --read_fps <forward> <reverse> <unpaired>'
            logging.critical(msg)
            sys.exit(1)
    return list(read_filetypes)[0]


def check_input_count_files(args):
    for assembly_fp in args.assembly_fps:
        if get_sequence_filetype(assembly_fp) != 'FASTA':
            logging.critical(f'Following file is not valid fasta: {assembly_fp}')
            sys.exit(1)
    for mask_fp in args.mask_fps:
        with mask_fp.open('r') as fh:
            if fh.readline().rstrip().split() != ['contig_name', 'contig_scores']:
                logging.critical(f'The file {mask_fp} does not not look like a mask file')
                sys.exit(1)


def execute_command(command, check=True):
    logging.debug(f'Running: {command}')
    # This will abort a series of piped commands if any fails but all stderr will be
    # returned. Without pipefail, even when processes fail only the last process returncode is
    # seen (which could be 0).
    command = f'set -o pipefail; {command}'
    result = subprocess.run(command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True,
                            encoding='utf-8',
                            executable=shutil.which('bash'))
    if check and result.returncode != 0:
        # TODO: are we happy with logging over multiple lines?
        logging.critical(f'Failed to run command: {result.args}')
        logging.critical(f'stdout: {result.stdout}')
        logging.critical(f'stderr: {result.stderr}')
        sys.exit(1)
    return result


def index_assembly(assembly_fp, temp_directory):
    log_newline()
    logging.info(f'Building a Bowtie2 index for {assembly_fp}')
    index_fp = pathlib.Path(temp_directory, assembly_fp)
    execute_command(f'bowtie2-build {assembly_fp} {index_fp}')
    return index_fp


def map_illumina_reads(index_fp, read_fps, temp_directory, threads, read_filetype):
    logging.info(f'Aligning Illumina reads to with Bowtie2')
    bam_fp = pathlib.Path(temp_directory, f'{index_fp.stem}.bam')
    command = f'bowtie2 --threads {threads} --sensitive -X 1000 -x {index_fp} '
    if read_filetype == 'FASTA':
        command += '-f '
    if len(read_fps) == 1:
        command += f'-U {read_fps[0]} '
    elif len(read_fps) == 2:
        command += f'-1 {read_fps[0]} -2 {read_fps[1]} '
    elif len(read_fps) == 3:
        command += f'-1 {read_fps[0]} -2 {read_fps[1]} -U {read_fps[2]}'
    command += f'| samtools sort -o {bam_fp}'
    execute_command(command)
    execute_command(f'samtools index {bam_fp}')
    return bam_fp


def map_long_reads(assembly_fp, read_fps, temp_directory, threads):
    logging.info(f'Aligning long reads to {assembly_fp} with Minimap2')
    bam_fp = pathlib.Path(temp_directory, f'{assembly_fp.stem}.bam')
    read_fps_str = ' '.join(str(rfp) for rfp in read_fps)
    command = f'minimap2 -t {threads} -a -x map-ont {assembly_fp} {read_fps_str} '
    command += f'| samtools sort -o {bam_fp}'
    execute_command(command)
    execute_command(f'samtools index {bam_fp}')
    return bam_fp


def get_base_scores_from_mpileup(assembly_fp, bam_fp):
    # TODO: split this over multiple threads, by using samtools' -r option
    log_newline()
    logging.info(f'Get base-level metrics for {assembly_fp} using Samtools mpileup')
    command = f'samtools mpileup -A -B -Q0 -vu -t INFO/AD -f {assembly_fp} {bam_fp}'
    mpileup_output = execute_command(command).stdout
    return get_base_scores_from_mpileup_output(mpileup_output)


def get_base_scores_from_mpileup_output(mpileup_output):
    """
    This function returns per-base scores as determined by the allelic depth (AD) information
    provided by samtools mpileup. Specifically, the score is what fraction of the reads at a
    position match the assembly's base(s), so higher is better.
    """
    scores = {}  # key = contig name, value = list of scores
    total_size = 0
    ad_regex = re.compile(r';AD=([\d,]+);')

    for line in mpileup_output.splitlines():
        if line.startswith('##contig='):
            contig_name = re.search(r'<ID=([^, ]+)', line).group(1)
            contig_length = int(re.search(r'length=(\d+)', line).group(1))
            total_size += contig_length
            logging.debug(f'contig {contig_name}: {contig_length} bp')
            scores[contig_name] = [None] * contig_length
        elif line.startswith('#'):
            continue
        else:
            parts = line.split('\t')
            contig_name = parts[0]
            pos = int(parts[1]) - 1  # use 0-based indexing
            length = len(parts[3])
            allele_depths = [int(x) for x in ad_regex.search(line).group(1).split(',')]
            ref_depth = allele_depths[0]
            total_depth = sum(allele_depths)
            if total_depth == 0:
                ref_frac = 0.0
            else:
                ref_frac = ref_depth / total_depth

            # If this VCF line covers multiples bases (is the case for indels), the score applies
            # to each base in the range. When multiple VCF lines cover the same base (e.g. an indel
            # covering a few and a substitution in the same range) then the lower (worst) score is
            # kept for that base.
            for i in range(length):
                if scores[contig_name][pos + i] is None:
                    scores[contig_name][pos + i] = ref_frac
                elif ref_frac < scores[contig_name][pos + i]:
                    scores[contig_name][pos + i] = ref_frac

    # Positions that didn't get a base might have no coverage, which is very bad, so they are given
    # a score of zero.
    for contig_name in list(scores.keys()):
        scores[contig_name] = [x if x is not None else 0.0 for x in scores[contig_name]]

    # TODO: maybe 'blur' the scores a bit, so we end up masking not just individual bad bases but
    # rather regions around bad bases?

    return scores, total_size


def get_score_threshold(scores, percentile):
    """
    Returns the given percentile of the scores using the nearest-rank method.
    https://en.wikipedia.org/wiki/Percentile#The_nearest-rank_method
    """
    if not scores:
        return 0.0
    sorted_scores = []
    for contig_scores in scores.values():
        sorted_scores += contig_scores
    sorted_scores = sorted(sorted_scores)
    fraction = percentile / 100.0
    rank = int(math.ceil(fraction * len(sorted_scores)))
    if rank == 0:
        return sorted_scores[0]
    return sorted_scores[rank - 1]


def write_mask_file(scores, min_score_threshold, assembly_fp, total_size):
    mask_fp = f'{assembly_fp}.mask'
    total_masked_bases = 0
    with open(mask_fp, 'wt') as mask:
        mask.write('\t'.join(('contig_name', 'contig_scores')))
        mask.write('\n')
        for contig_name, contig_scores in scores.items():
            mask.write(contig_name)
            mask.write('\t')
            masked_bases = [i for i, s in enumerate(contig_scores) if s < min_score_threshold]
            total_masked_bases += len(masked_bases)
            mask.write(','.join((str(x) for x in masked_bases)))
            mask.write('\n')
    logging.info(f'Masked {total_masked_bases:,d} out of {total_size:,d} bases in {assembly_fp}')


def load_mask_file(mask_fp):
    masked_positions = {}
    with open(mask_fp, 'rt') as mask:
        # Skip header
        mask.readline()
        for line in mask:
            parts = line.rstrip('\n').split('\t')
            contig_name = parts[0]
            if parts[1]:
                positions = set(int(x) for x in parts[1].split(','))
            else:
                positions = set()
            masked_positions[contig_name] = positions
    return masked_positions


def default_thread_count():
    return min(multiprocessing.cpu_count(), 8)


class Snp(object):
    def __init__(self, show_snps_line):
        parts = show_snps_line.strip().split('\t')
        self.a1_pos = int(parts[0])
        self.a1_base = parts[1]
        self.a2_base = parts[2]
        assert len(self.a1_base) == 1
        assert len(self.a2_base) == 1
        self.a2_pos = int(parts[3])
        self.a1_strand = int(parts[6])
        self.a2_strand = int(parts[7])
        self.dist_to_contig_end = int(parts[5])
        self.a1_contig = parts[8]
        self.a2_contig = parts[9]

    def __str__(self):
        a1_strand = '+' if self.a1_strand == 1 else '-'
        a2_strand = '+' if self.a2_strand == 1 else '-'
        return f'{self.a1_contig}{a1_strand}:{self.a1_pos}, ' \
               f'{self.a2_contig}{a2_strand}:{self.a2_pos}: ' \
               f'{self.a1_base} -> {self.a2_base}'

    def data_tuple(self):
        assert self.a1_strand == 1
        return (self.a1_contig, self.a1_pos, self.a1_strand,
                self.a2_contig, self.a2_pos, self.a2_strand,
                self.a1_base, self.a2_base)

    def __lt__(self, other):
        return self.data_tuple() < other.data_tuple()

    def __eq__(self, other):
        return self.data_tuple() == other.data_tuple()

    def __hash__(self):
        return hash(self.data_tuple())


def get_snp_count_from_nucmer(assembly_1, mask_1, assembly_2, mask_2, prefix):
    log_newline()
    logging.info(f'Aligning {assembly_1} and {assembly_2} using numcer')
    execute_command(f'nucmer --prefix={prefix} {assembly_1} {assembly_2}')
    show_snps_output = execute_command(f'show-snps -CrTH {prefix}.delta').stdout
    snps = sorted(Snp(x) for x in show_snps_output.splitlines())
    before_mask_count = len(snps)
    snps = [x for x in snps
            if x.a1_pos not in mask_1[x.a1_contig] and x.a2_pos not in mask_2[x.a2_contig]]
    for snp in snps:
        logging.debug(str(snp))
    after_mask_count = len(snps)
    mask_diff = before_mask_count - after_mask_count
    logging.info(f'Found {before_mask_count:,d} SNPs, {mask_diff:,d} of which were masked out '
                 f'leaving {after_mask_count:,d} SNPs')
    return after_mask_count


def get_pairwise_snp_count(assembly_1, mask_1, assembly_2, mask_2, temp_dir):
    snps_1_vs_2 = get_snp_count_from_nucmer(assembly_1, mask_1, assembly_2, mask_2,
                                            pathlib.Path(temp_dir, 'out1'))
    snps_2_vs_1 = get_snp_count_from_nucmer(assembly_2, mask_2, assembly_1, mask_1,
                                            pathlib.Path(temp_dir, 'out2'))
    final_count = min(snps_1_vs_2, snps_2_vs_1)
    log_newline()
    logging.info(f'Final SNP count between {assembly_1} and {assembly_2}: {final_count:,d}')
    return final_count


if __name__ == '__main__':
    main()
