#!/usr/bin/env python3
import argparse
import pathlib
import logging
import subprocess
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

    parser_align = subparser.add_parser('align')
    parser_align.add_argument('--assembly_fp', required=True, type=pathlib.Path,
            help='Input assembly filepath, space separated')
    parser_align.add_argument('--mask_fp', type=pathlib.Path,
            help='Input masking filepath, space separated')

    # TODO: Perform additional argument parsing, checking
    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        print('\n', end='')
        parser.error('command options include mask or align')

    return args


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
    # Create temp working directory
    with tempfile.TemporaryDirectory() as dh:
        index_fp = index_assembly(args.assembly_fp, dh)
        map_reads(index_fp, args.read_fps, args.read_type, dh)


def run_align(args):
    check_input_align_files(args)
    pass


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
    logger.setLevel(logging.INFO)


def check_dependencies():
    # TODO: do we need to check for version as well?
    command_template = 'which %s'
    # TODO: add dependencies
    dependencies = ['samtools', 'bowtie2']
    for dependency in dependencies:
        result = execute_command(command_template % dependency, quiet=True)
        if result.returncode != 0:
            logging.critical('Could not find dependency %s' % dependency)
            logging.critical('%s', result.stderr)


def check_input_mask_files(args):
    # TODO: check that files look like FASTQ and FASTA
    pass


def check_input_align_files(args):
    # TODO: check that files look like our mask format and FASTA
    pass


def execute_command(command, quiet=False):
    logging.info(command)
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=True, encoding='utf-8')
    if not quiet and result.returncode != 0:
        # TODO: are we happy with logging over multiple lines?
        logging.critical('Failed to run command: %s', result.args)
        logging.critical('stdout: %s', result.stdout)
        logging.critical('stderr: %s', result.stderr)

    return result


def index_assembly(assembly_fp, temp_directory):
    command_template = 'bowtie2-build %s %s'
    execute_command(command_template % (assembly_fp, temp_directory))
    return pathlib.Path(temp_directory, assembly_fp)


def map_reads(index_fp, reads_fp, read_type, temp_directory):
    sam_fp = pathlib.Path(temp_directory, '%s.sam' % index_fp.stem)
    # TODO: get correct commands
    if read_type == 'illumina':
        # TODO: is it worth refactoring to avoid verbosity?
        if len(reads_fp) == 1:
            command_template = 'bowtie2 --sensistive -X 2000 -x -U %s -S %s'
            command = command_template % (reads_fp, sam_fp)
        elif len(reads_fp) == 2:
            command_template = 'bowtie2 --sensistive -X 2000 -x -1 %s -2 %s -S %s'
            command = command_template % (*reads_fp, sam_fp)
    elif read_type == 'long':
        command_template = 'freebayes -U %s > %s'
        command = command_template % (reads_fp, sam_fp)
    execute_command(command)
    return sam_fp


if __name__ == '__main__':
    main()
