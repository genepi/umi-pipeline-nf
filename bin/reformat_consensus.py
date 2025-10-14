"""
This is a modified version of the code present in:
https://github.com/nanoporetech/pipeline-umi-amplicon/blob/master/lib/umi_amplicon_tools/reformat_consensus.py
"""

import argparse
import logging
import sys
import os

import pysam


def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Command line interface to telemap"
    parser = argparse.ArgumentParser(
        description=usage, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-l",
        "--log",
        dest="log",
        choices=[
            "DEBUG",
            "INFO",
            "WARNING",
            "ERROR",
            "CRITICAL",
            "debug",
            "info",
            "warning",
            "error",
            "critical",
        ],
        default="INFO",
        help="Print debug information",
    )

    parser.add_argument(
        "--consensus_fasta", dest="CONS_FASTA", required=True, type=str, help="consensus fasta file"
    )

    parser.add_argument(
        "-o", "--output", dest="OUTPUT", required=True, help="Output folder"
    )

    parser.add_argument(
        "-t", "--threads", dest="THREADS", type=int, default=1, help="Number of threads."
    )

    args = parser.parse_args(argv)

    return args

def get_read_seq(read):
    return read.name.split(";seq=")[1].split(";")[0]

def get_read_qual(read):
    return read.name.split(";qual=")[1].split(";seqs=")[0]

def get_read_name(read):
    return read.name.split(";")[0]

def parse_stdin(args):
    consensus_filename = args.CONS_FASTA
    output_folder = args.OUTPUT
    cluster_filename = os.path.join(output_folder, "final.fastq")

    with open(cluster_filename, "w") as out, pysam.FastxFile(consensus_filename) as reads:
        for read in reads:
            read_name = get_read_name(read)
            read_seq = get_read_seq(read)
            read_qual = get_read_qual(read)
            print("@{}".format(read_name), file=out)
            print(read_seq, file=out)
            print("+", file=out)
            print(read_qual, file=out)


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to telemap.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % args.log.upper())
    logging.basicConfig(level=numeric_level, format="%(message)s")

    parse_stdin(args)


if __name__ == "__main__":
    main()
