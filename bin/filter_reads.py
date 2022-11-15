import argparse
import logging
import os
import sys

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
        "-t", "--threads", dest="THREADS", type=int, default=1, help="Number of threads."
    )

    parser.add_argument(
        "--min_overlap",
        dest="MIN_OVERLAP",
        type=float,
        default=0.9,
        help="Min overlap with target region",
    )

    parser.add_argument(
        "--include_secondary_reads",
        dest="INCL_SEC",
        action="store_true",
        help="Include secondary alignments",
    )

    parser.add_argument(
        "-o", "--output", dest="OUT", type=str, required=False, help="Output folder"
    )

    parser.add_argument(
        "--tsv",
        dest="TSV",
        action="store true",
        help="Write tsv file containing filtering stats"
    )

    parser.add_argument("BED", type=str, nargs=1, help="BED file")

    parser.add_argument(
        "BAM", type=str, nargs="?", default="/dev/stdin", help="BAM file"
    )

    parser.add_argument(
        "--output_format",
        dest="OUT_FORMAT",
        type=str,
        help="Choose fastq or fasta",
        default="fasta"
    )

    args = parser.parse_args(argv)

    return args


def parse_bed(bed_regions):
    with open(bed_regions) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                logging.warning("Ignoring BED entry: {}".format(line))
                continue

            region = {
                "chr": cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
                "name": cols[3],
            }
            return region


def write_read(read, output, region, type, format):
    output_fastx = os.path.join(
        output, "{}_{}.{}".format(region["name"], type, format)
    )

    # see if appending line is no problem by running it with nextflow (Otherwise delete files before appending for the first time)
    with open(output_fastx, "a") as out_f:
        if format == "fasta":
            write_fasta(read, out_f)
        elif format == "fastq":
            write_fastq(read, out_f)
        else:
            raise RuntimeError("specified format incorrect: {}".format(format))


def write_fasta(read, out_f):
    read_strand = "-"
    if read.is_reverse:
        print(
            ">{};strand={}".format(read.query_name, read_strand), file=out_f
        )
        print(read.get_forward_sequence(), file=out_f)
    else:
        read_strand = "+"
        print(
            ">{};strand={}".format(read.query_name, read_strand), file=out_f
        )
        print(read.query_sequence, file=out_f)


def write_fastq(read, out_f):
    read_strand = "-"
    if read.is_reverse:
        print(
            "@{};strand={}".format(read.query_name, read_strand), file=out_f
        )
        print(read.get_forward_sequence(), file=out_f)
        print("+", file=out_f)
        print(pysam.qualities_to_qualitystring(
            read.get_forward_qualities()), file=out_f)
    else:
        read_strand = "+"
        print(
            "@{};strand={}".format(read.query_name, read_strand), file=out_f
        )
        print(read.query_sequence, file=out_f)
        print("+", file=out_f)
        print(pysam.qualities_to_qualitystring(
            read.query_qualities), file=out_f)


def filter_reads(args):
    bed_regions = args.BED[0]
    bam_file = args.BAM
    max_clipping = 250
    min_overlap = args.MIN_OVERLAP
    incl_sec = args.INCL_SEC
    output = args.OUT
    out_format = args.OUT_FORMAT
    tsv = args.TSV
    stats_out_filename = "umi_filter_reads_stats"

    n_non_reads = 0
    n_unmapped = 0
    n_concatamer = 0
    n_short = 0
    n_ontarget = 0
    n_reads_region = 0
    n_supplementary = 0
    n_secondary = 0
    n_total = 0

    filtered_perc = 0
    unmapped_perc = 0
    ontarget_perc = 0
    concatermer_perc = 0
    short_perc = 0

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        region = parse_bed(bed_regions)
        region_length = region["end"] - region["start"]
        logging.info("Region: {}".format(region["name"]))

        for read in bam.fetch(
            until_eof=True
        ):

            if (read.query_sequence is None):
                n_non_reads += 1
                continue

            n_total += 1
            if read.is_unmapped:
                n_unmapped += 1
                write_read(read, output, region, "unmapped", out_format)
                continue

            if read.is_secondary:
                if not incl_sec:
                    n_secondary += 1
                    write_read(read, output, region, "secondary", out_format)
                    continue

            if read.is_supplementary:
                n_supplementary += 1
                write_read(read, output, region, "supplementary", out_format)
                continue

            n_ontarget += 1
            if read.query_alignment_length < (
                read.query_length - 2 * max_clipping
            ):
                n_concatamer += 1
                write_read(read, output, region, "concatamer", out_format)
                continue

            if read.reference_length < (region_length * min_overlap):
                n_short += 1
                write_read(read, output, region, "short", out_format)
                continue
            n_reads_region += 1
            write_read(read, output, region, "filtered", out_format)

    if tsv:
        stats_out_filename = os.path.join(
            output, "{}.tsv".format(stats_out_filename))
        write_tsv(n_total, n_unmapped, n_secondary, n_supplementary, n_ontarget,
                  n_concatamer, n_short, n_reads_region, incl_sec, stats_out_filename, region)


def write_tsv(n_total, n_unmapped, n_secondary, n_supplementary, n_ontarget, n_concatamer, n_short, n_reads_region, incl_sec, stats_out_filename, region):
    if n_total > 0:
        if incl_sec:
            filtered_perc = 100 * n_supplementary // n_total
        else:
            filtered_perc = 100 * (n_secondary + n_supplementary) // n_total

        unmapped_perc = 100 * n_unmapped // n_total
        secondary_perc = 100 * n_secondary // n_total
        supplementary_perc = 100 * n_supplementary // n_total
        ontarget_perc = 100 * n_ontarget // n_total

        if ontarget_perc > 0:
            concatermer_perc = 100 * n_concatamer // n_ontarget
            short_perc = 100 * n_short // n_ontarget

    with open(stats_out_filename, "a") as out_f:
        print(
            "format",
            "region",
            "reads_found",
            "reads_unmapped",
            "reads_secondary",
            "reads_supplementary",
            "reads_on_target",
            "reads_concatamer",
            "reads_short",
            "reads_filtered",
            "include_secondary",
            sep="\t",
            file=out_f
        )
        print(
            "count",
            region,
            n_total,
            n_unmapped,
            n_secondary,
            n_supplementary,
            n_ontarget,
            n_concatamer,
            n_short,
            n_reads_region,
            incl_sec,
            sep="\t",
            file=out_f
        )
        print(
            "%",
            region,
            "100%",
            unmapped_perc,
            secondary_perc,
            supplementary_perc,
            ontarget_perc,
            concatermer_perc,
            short_perc,
            filtered_perc,
            incl_sec,
            sep="\t",
            file=out_f
        )


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

    filter_reads(args)


if __name__ == "__main__":
    main()
