#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com), Stephen Watts, Alex Tokolyi
https://github.com/rrwick/Pairwise-SNP-counter

This is a command line tool for counting the single-nucleotide differences (a.k.a. SNPs) between
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
    subparser = parser.add_subparsers(metavar='mask align', dest='command')

    parser_mask = subparser.add_parser('mask')
    parser_mask.add_argument('--assembly_fp', required=True, type=pathlib.Path,
                             help='Input assembly filepath')
    parser_mask.add_argument('--read_fps', required=True, nargs='+', type=pathlib.Path,
                             help='Input read filepaths, space separated')
    parser_mask.add_argument('--read_type', required=True, choices=['illumina', 'long'],
                             help='Read type of input reads. [choices: illumina, long]')
    parser_mask.add_argument('--threads', required=False, type=int, default=default_thread_count(),
                             help='Number of threads')
    # TODO: option to specify temp directory
    parser_mask.add_argument('--exclude', required=False, type=float, default=2.0,
                             help='Percentage of assembly bases to exclude')

    parser_align = subparser.add_parser('align')
    parser_align.add_argument('--assembly_fps', required=True, nargs='+', type=pathlib.Path,
                              help='Input assembly filepaths, space separated')
    parser_align.add_argument('--mask_fps', nargs='+', type=pathlib.Path,
                              help='Input masking filepaths, space separated')

    args = parser.parse_args()
    if not args.command:
        # TODO: print better help info. see samtools for an example with subcommands
        parser.print_help()
        print('\n', end='')
        parser.error('command options include mask or align')

    # TODO: Perform additional argument parsing, checking
    check_parsed_file_exists(args.assembly_fp, parser)
    if args.command == 'mask':
        for read_fp in args.read_fps:
            check_parsed_file_exists(read_fp, parser)

        if args.read_type == 'illumina':
            if len(args.read_fps) > 2:
                parser.error('--read_fps takes no more than two illumina read sets')
        elif args.read_type == 'long':
            if len(args.read_fps) > 1:
                parser.error('--read_fps takes only a single long read set')

    if args.command == 'align':
        if len(args.assembly_fps) < 2:
            parser.error('Two or more assemblies are required, got {len(args.assembly_fps}')
        if not args.mask_fps:
            args.mask_fps = [f'{fp}.mask' for fp in args.assembly_fps]
        if len(args.assembly_fps) != len(args.mask_fps):
            parser.error('Need the same number of masks as assemblies')
        check_parsed_file_exists(args.mask_fp, parser)

    return args


def check_parsed_file_exists(filepath, parser):
    # Check that the argument has been set; is not None
    if filepath and not filepath.exists():
        parser.error(f'Input file {filepath} does not exist')


def main():
    # Get commandline arguments and initialise
    args = get_arguments()
    initialise_logging()
    check_dependencies()

    # Execute requested stage
    if args.command == 'mask':
        run_mask(args)
    elif args.command == 'align':
        run_align(args)


def run_mask(args):
    read_filetype = check_input_mask_files(args)
    with tempfile.TemporaryDirectory() as dh:
        # Map reads to assembly
        if args.read_type == 'illumina':
            index_fp = index_assembly(args.assembly_fp, dh)
            bam_fp = map_illumina_reads(index_fp, args.read_fps, dh, args.threads, read_filetype)
        else:
            assert args.read_type == 'long'
            bam_fp = map_long_reads(args.assembly_fp, args.read_fps, dh, args.threads)
        scores = get_base_scores_from_mpileup(args.assembly_fp, bam_fp)
        min_score_threshold = get_score_threshold(scores, args.exclude)
        write_mask_file(scores, min_score_threshold, args.assembly_fp)


def run_align(args):
    check_input_align_files(args)
    counts = {}
    for i in range(len(args.assembly_fps)):
        assembly_1 = args.assembly_fps[i]
        mask_1 = load_mask_file(args.mask_fps[i])
        for j in range(i+1, len(args.assemblies)):
            assembly_2 = args.assembly_fps[j]
            mask_2 = load_mask_file(args.mask_fps[j])
            snp_count = get_pairwise_snp_count(assembly_1, mask_1, assembly_2, mask_2)
            counts[(assembly_1, assembly_2)] = snp_count
            counts[(assembly_2, assembly_1)] = snp_count


def initialise_logging():
    # TODO: command line arguments for logging; save to file? print to stdout?
    # Set up loggers
    log_filehandler = logging.FileHandler('run.log', mode='w')
    log_streamhandler = logging.StreamHandler()
    log_format = logging.Formatter(fmt='%(asctime)s %(message)s', datefmt='%d/%m/%Y %H:%M:%S')
    log_filehandler.setFormatter(log_format)
    log_streamhandler.setFormatter(log_format)

    # Add log handles to root logger and set log level
    logger = logging.getLogger()
    logger.addHandler(log_filehandler)
    logger.addHandler(log_streamhandler)
    # TODO: expose this option to the command line
    logger.setLevel(logging.DEBUG)

    # TODO: use colours in stdout
    # https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output


def log_newline():
    logging_formatter = logging.getLogger('').handlers[0].formatter
    for h in logging.getLogger('').handlers:
        h.setFormatter(logging.Formatter(fmt=''))
    logging.info('')
    for h in logging.getLogger('').handlers:
        h.formatter = logging_formatter


def check_dependencies():
    # TODO: add all other dependencies
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
                                 'vrequired': '1.0'}}
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
    logging.info(f'{args.assembly_fp} ({assembly_filetype})')

    read_filetypes = set()
    for read_fp in args.read_fps:
        read_filetype = get_sequence_filetype(read_fp)
        logging.info(f'{read_fp} ({read_filetype})')
        read_filetypes.add(read_filetype)
    if args.read_type == 'illumina' and len(read_filetypes) > 1:
        logging.critical('Read files must be all FASTQ or all FASTA (not a mixture of both)')
        sys.exit(1)
    assert len(read_filetypes) == 1
    return list(read_filetypes)[0]


def check_input_align_files(args):
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
    ad_regex = re.compile(r';AD=([\d,]+);')

    for line in mpileup_output.splitlines():
        if line.startswith('##contig='):
            contig_name = re.search(r'ID=(\w+)', line).group(1)
            contig_length = int(re.search(r'length=(\d+)', line).group(1))
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
            ref_frac = ref_depth / sum(allele_depths)

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

    return scores


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


def write_mask_file(scores, min_score_threshold, assembly_fp):
    mask_fp = f'{assembly_fp}.mask'
    with open(mask_fp, 'wt') as mask:
        mask.write('\t'.join(('contig_name', 'contig_scores')))
        mask.write('\n')
        for contig_name, contig_scores in scores.items():
            mask.write(contig_name)
            mask.write('\t')
            mask.write(','.join((str(i) for i, s in enumerate(contig_scores)
                                 if s < min_score_threshold)))
            mask.write('\n')


def load_mask_file(mask_fp):
    masked_positions = {}
    with open(mask_fp, 'rt') as mask:
        # Skip header
        mask.readline()
        for line in mask:
            parts = line.rstrip('\n').split('\t')
            contig_name = parts[0]
            if parts[1]:
                positions = [int(x) for x in parts[1].split(',')]
            else:
                positions = []
            masked_positions[contig_name] = positions
    return masked_positions


def default_thread_count():
    return min(multiprocessing.cpu_count(), 8)


def get_pairwise_snp_count(assembly_1, mask_1, assembly_2, mask_2):
    return 0  # TEMP


if __name__ == '__main__':
    main()
