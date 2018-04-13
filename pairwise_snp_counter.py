#!/usr/bin/env python3
import argparse
import pathlib
import logging
import subprocess


def get_arguments():
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(metavar='mask align', dest='command')

    parser_mask = subparser.add_parser('mask')
    parser_mask.add_argument('--input_fps', required=True, nargs='+', type=pathlib.Path,
            help='Input file paths')

    parser_align = subparser.add_parser('align')
    parser_align.add_argument('--input_fps', required=True, nargs='+', type=pathlib.Path,
            help='Input file paths')

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

    # Execute requested stage
    if args.command == 'mask':
        run_mask(args)
    elif args.command == 'align':
        run_align(args)


def run_mask(args):
    pass


def run_align(args):
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


def execute_command(command):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
                            encoding='utf-8')
    if result.returncode != 0:
        # TODO: are we happy with logging over multiple lines?
        logging.critical('Failed to run command: %s', result.args)
        logging.critical('stdout: %s', result.stdout)
        logging.critical('stderr: %s', result.stderr)


if __name__ == '__main__':
    main()
