#!/usr/bin/env python3
import argparse
import pathlib


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
    # Get commandline arguments
    args = get_arguments()

    # Execute requested stage
    if args.command == 'mask':
        run_mask(args)
    elif args.command == 'align':
        run_align(args)


def run_mask(args):
    pass


def run_align(args):
    pass


if __name__ == '__main__':
    main()
