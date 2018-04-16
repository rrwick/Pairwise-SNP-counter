#!/usr/bin/env python3
import argparse
from distutils.version import LooseVersion
import logging
import math
import pathlib
import re
import shutil
import subprocess
import tempfile
import multiprocessing


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
    parser_align.add_argument('--assembly_fp', required=True, type=pathlib.Path,
                              help='Input assembly filepath, space separated')
    parser_align.add_argument('--mask_fp', type=pathlib.Path,
                              help='Input masking filepath, space separated')

    args = parser.parse_args()
    if not args.command:
        # TODO: print better help info. see samtools for an example with subcommands
        parser.print_help()
        print('\n', end='')
        parser.error('command options include mask or align')

    # TODO: Perform additional argument parsing, checking
    check_parsed_file_exists(args.assembly_fp, parser)
    if args.command == 'align':
        check_parsed_file_exists(args.mask_fp, parser)

    if args.command == 'mask':
        for read_fp in args.read_fps:
            check_parsed_file_exists(read_fp, parser)

        if args.read_type == 'illumina':
            if len(args.read_fps) > 2:
                parser.error('--read_fps takes no more than two illumina read sets')
        elif args.read_type == 'long':
            if len(args.read_fps) > 1:
                parser.error('--read_fps takes only a single long read set')

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
    check_input_mask_files(args)
    with tempfile.TemporaryDirectory() as dh:
        # Map reads to assembly
        if args.read_type == 'illumina':
            index_fp = index_assembly(args.assembly_fp, dh)
            bam_fp = map_illumina_reads(index_fp, args.read_fps, dh, args.threads)
        else:
            assert args.read_type == 'long'
            bam_fp = map_long_reads(args.assembly_fp, args.read_fps, dh, args.threads)
        scores = get_base_scores_from_mpileup(args.assembly_fp, bam_fp)
        min_score_threshold = get_score_threshold(scores, args.exclude)
        write_mask_file(scores, min_score_threshold, args.assembly_fp)


def run_align(args):
    check_input_align_files(args)


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


def check_dependencies():
    # TODO: add all other dependencies
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
            exit(1)
        result = execute_command(version_data['vcommand'], check=False)
        try:
            version = version_data['vregex'].search(result.stdout).group(1)
        except AttributeError:
            # TODO: should we have an option to skip dependency check?
            logging.critical(f'Unable to determine version for {dependency}')
            exit(1)
        if LooseVersion(version) < LooseVersion(version_data['vrequired']):
            logging.critical(f'{dependency} version {version_data["vrequired"]} or high is required')
            exit(1)


def check_input_mask_files(args):
    # TODO: check that files look like FASTQ and FASTA
    # TODO: if Illumina reads, check they are in the right order (1, 2, unpaired)
    pass


def check_input_align_files(args):
    # TODO: check that files look like our mask format and FASTA
    pass


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
        exit(1)
    return result


def index_assembly(assembly_fp, temp_directory):
    index_fp = pathlib.Path(temp_directory, assembly_fp)
    execute_command(f'bowtie2-build {assembly_fp} {index_fp}')
    return index_fp


def map_illumina_reads(index_fp, read_fps, temp_directory, threads):
    bam_fp = pathlib.Path(temp_directory, f'{index_fp.stem}.bam')
    command = f'bowtie2 --threads {threads} --sensitive -X 1000 -x {index_fp} '
    if len(read_fps) == 1:
        command += f'-U {read_fps[0]} '
    elif len(read_fps) == 2:
        command += f'-1 {read_fps[0]} -2 {read_fps[1]} '
    command += f'| samtools sort -o {bam_fp}'
    execute_command(command)
    execute_command(f'samtools index {bam_fp}')
    return bam_fp


def map_long_reads(assembly_fp, read_fps, temp_directory, threads):
    bam_fp = pathlib.Path(temp_directory, f'{assembly_fp.stem}.bam')
    read_fps_str = ' '.join(str(rfp) for rfp in read_fps)
    command = f'minimap2 -t {threads} -a -x map-ont {assembly_fp} {read_fps_str} '
    command += f'| samtools sort -o {bam_fp}'
    execute_command(command)
    execute_command(f'samtools index {bam_fp}')
    return bam_fp


def get_base_scores_from_mpileup(assembly_fp, bam_fp):
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
        for contig_name, contig_scores in scores.items():
            mask.write(contig_name)
            mask.write('\t')
            mask.write(','.join((str(i) for i, s in enumerate(contig_scores)
                                 if s < min_score_threshold)))
            mask.write('\n')


def default_thread_count():
    return min(multiprocessing.cpu_count(), 8)


if __name__ == '__main__':
    main()
