"""
This is a modified version of the code present in:
https://github.com/nanoporetech/pipeline-umi-amplicon/blob/master/lib/umi_amplicon_tools/extract_umis.py
"""

import argparse
import logging
import os

import edlib
import pysam
import sys


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
        "--max-error",
        dest="MAX_ERROR",
        type=int,
        default=2,
        help="Max edit distance for UMI",
    )
    parser.add_argument(
        "--adapter_length",
        dest="ADAPTER_LENGTH",
        type=int,
        default=200,
        help="Length of adapter",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="THREADS",
        type=int,
        default=1,
        help="Number of threads."
    )
    parser.add_argument(
        "--tsv",
        dest="TSV",
        action="store_true",
        help="write TSV output file"
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="OUT",
        type=str,
        required=False,
        help="Output directory"
    )
    parser.add_argument(
        "--output_format",
        dest="OUT_FORMAT",
        type=str,
        help="Choose fastq or fasta",
        default="fasta"
    )
    parser.add_argument(
        "--fwd-umi",
        dest="FWD_UMI",
        type=str,
        default="TTTVVVVTTVVVVTTVVVVTTVVVVTTT",
        help="Forward UMI sequence",
    )
    parser.add_argument(
        "--rev-umi",
        dest="REV_UMI",
        type=str,
        default="AAABBBBAABBBBAABBBBAABBBBAAA",
        help="Reverse UMI sequence",
    )
    parser.add_argument(
        "INPUT_FA",
        type=str,
        default="/dev/stdin",
        help="Filtered Reads"
    )

    args = parser.parse_args(argv)

    return args


def clip_entry(entry, umi_start_fwd, umi_start_rev, adapter_length, format):
    clip_fwd = umi_start_fwd
    if adapter_length == umi_start_rev:
        clip_rev = None
    else: 
        clip_rev = umi_start_rev - adapter_length
    
    entry.sequence = entry.sequence[clip_fwd:clip_rev]
    
    if format == "fastq":
        entry.quality = entry.quality[clip_fwd:clip_rev]
    
    return entry

def extract_umi(query_seq, query_qual, pattern, max_edit_dist, format, direction):
    umi_qual = None
    equalities = [("M", "A"), ("M", "C"), ("R", "A"), ("R", "G"), ("W", "A"), ("W", "T"), ("S", "C"), ("S", "G"), ("Y", "C"), ("Y", "T"), ("K", "G"), ("K", "T"), ("V", "A"), ("V", "C"),
                  ("V", "G"), ("H", "A"), ("H", "C"), ("H", "T"), ("D", "A"), ("D", "G"), ("D", "T"), ("B", "C"), ("B", "G"), ("B", "T"), ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T"),
                  ("m", "a"), ("m", "c"), ("r", "a"), ("r", "g"), ("w", "a"), ("w", "t"), ("s", "c"), ("s", "g"), ("y", "c"), ("y", "t"), ("k", "g"), ("k", "t"), ("v", "a"), ("v", "c"),
                  ("v", "g"), ("h", "a"), ("h", "c"), ("h", "t"), ("d", "a"), ("d", "g"), ("d", "t"), ("b", "c"), ("b", "g"), ("b", "t"), ("n", "a"), ("n", "c"), ("n", "g"), ("n", "t"), 
                  ("a", "A"), ("c", "C"), ("t", "T"), ("g", "G")]

    result = edlib.align(
        pattern,
        query_seq,
        task="path",
        mode="HW",
        k=max_edit_dist,
        additionalEqualities=equalities,
    )
    if result["editDistance"] == -1:
        return None, None, None, None

    edit_dist = result["editDistance"]
    locs = result["locations"][0]
    umi_start_pos = locs[0]
    umi_end_pos = locs[1] + 1
    umi = query_seq[umi_start_pos:umi_end_pos]
    
    if direction == "fwd":
        umi_start = umi_start_pos
    else:
        umi_start = umi_end_pos

    if format == "fastq":
        umi_qual = query_qual[locs[0]:locs[1]+1]

    return edit_dist, umi, umi_qual, umi_start


def extract_adapters(entry, adapter_length, format):
    read_5p_seq = None
    read_5p_qual = None
    read_3p_seq = None
    read_3p_qual = None

    if len(entry.sequence) > adapter_length:
        read_5p_seq = entry.sequence[:adapter_length]
        read_3p_seq = entry.sequence[-adapter_length:]
        
        if format == "fastq":
            read_5p_qual = entry.quality[:adapter_length]
            read_3p_qual = entry.quality[-adapter_length:]

    return read_5p_seq, read_3p_seq, read_5p_qual, read_3p_qual


def get_read_name(entry):
    return entry.name.split(";")[0]


def get_read_strand(entry):
    strand = entry.name.split("strand=")
    if len(strand) > 1:
        return strand[1]
    # second extraction only includes positive strand!
    else:
        return "+"


def combine_umis_fasta(seq_5p, seq_3p, strand):
    if strand == "+":
        return seq_5p + seq_3p
    else:
        return rev_comp(seq_3p) + rev_comp(seq_5p)


def combine_umis_fastq(seq_5p, seq_3p, qual_5p, qual_3p, strand):
    if strand == "+":
        return (seq_5p + seq_3p), (qual_5p + qual_3p)
    else:
        return (rev_comp(seq_3p) + rev_comp(seq_5p)), (rev_comp_qual(qual_3p) + rev_comp_qual(qual_5p))


def write_fasta(
    entry,
    strand,
    result_5p_fwd_umi_dist,
    result_3p_rev_umi_dist,
    result_5p_fwd_umi_seq,
    result_3p_rev_umi_seq,
    out,
):
    seq = combine_umis_fasta(result_5p_fwd_umi_seq,
                             result_3p_rev_umi_seq, strand)
    print(
        ">{};strand={};umi_fwd_dist={};umi_rev_dist={};umi_fwd_seq={};umi_rev_seq={};seq={}".format(
            entry.name,
            strand,
            result_5p_fwd_umi_dist,
            result_3p_rev_umi_dist,
            result_5p_fwd_umi_seq,
            result_3p_rev_umi_seq,
            entry.sequence,
        ),
        file=out,
    )
    print(seq, file=out)


def write_fastq(
    entry,
    strand,
    result_5p_fwd_umi_dist,
    result_3p_rev_umi_dist,
    result_5p_fwd_umi_seq,
    result_3p_rev_umi_seq,
    result_5p_fwd_umi_qual,
    result_3p_rev_umi_qual,
    out_f,
):
    seq, qual = combine_umis_fastq(result_5p_fwd_umi_seq, result_3p_rev_umi_seq,
                                   result_5p_fwd_umi_qual, result_3p_rev_umi_qual, strand)
    print(
        "@{};strand={};umi_fwd_dist={};umi_rev_dist={};umi_fwd_seq={};umi_rev_seq={};seq={};qual={}".format(
            entry.name,
            strand,
            result_5p_fwd_umi_dist,
            result_3p_rev_umi_dist,
            result_5p_fwd_umi_seq,
            result_3p_rev_umi_seq,
            entry.sequence,
            entry.quality
        ),
        file=out_f,
    )
    print(seq, file=out_f)
    print("+", file=out_f)
    print(qual, file=out_f)


def write_tsv(
        output_folder,
        output_file_name,
        max_pattern_dist,
        strand_stats,
        n_total,
        n_both_umi):

    fwd_rev_ratio = -1
    perc = -1.0
    if strand_stats["-"]:
        fwd_rev_ratio = strand_stats["+"] / strand_stats["-"]
    logging.info(
        "Found {} fwd and {} rev reads (ratio: {})".format(
            strand_stats["+"], strand_stats["-"], fwd_rev_ratio
        )
    )
    if n_total:
        perc = 100.0 * n_both_umi / n_total
        logging.info(
            "{}% of reads contained both UMIs with max {} mismatches".format(
                perc, max_pattern_dist
            )
        )
    tsv_file = os.path.join(output_folder, "{}.tsv".format(output_file_name))
    with open(tsv_file, "w") as tsv_f:
        print(
            "output_file",
            "max_pattern_distance",
            "detected_forward_strands",
            "detected_reverse_strands",
            "ratio_fwd_rvs",
            "total_reads",
            "included_reads",
            "percent_of_total_reads",
            sep="\t",
            file=tsv_f
        )
        print(
            output_file_name,
            max_pattern_dist,
            strand_stats["+"],
            strand_stats["-"],
            fwd_rev_ratio,
            n_total,
            n_both_umi,
            perc,
            file=tsv_f,
            sep="\t"
        )


def extract_umis(
    args
):
    adapter_length = args.ADAPTER_LENGTH
    max_pattern_dist = args.MAX_ERROR
    output_folder = args.OUT
    tsv = args.TSV
    input_file = args.INPUT_FA
    umi_fwd = args.FWD_UMI
    umi_rev = args.REV_UMI
    output_file_name = "detected_umis"
    format = args.OUT_FORMAT

    output_file = os.path.join(
        output_folder, "{}.{}".format(output_file_name, format))

    n_total = 0
    n_both_umi = 0
    strand_stats = {"+": 0, "-": 0}
    with pysam.FastxFile(input_file) as fh, open(output_file, "w") as out:
        for entry in fh:
            strand = get_read_strand(entry)
            entry.name = get_read_name(entry)
            n_total += 1

            read_5p_seq, read_3p_seq, read_5p_qual, read_3p_qual = extract_adapters(
                entry, adapter_length, format
            )

            if read_5p_seq is None or read_3p_seq is None:
                continue

            strand_stats[strand] += 1

            # Extract fwd UMI
            result_5p_fwd_umi_dist, result_5p_fwd_umi_seq, result_5p_fwd_umi_qual, umi_start_fwd = extract_umi(
                read_5p_seq, read_5p_qual, umi_fwd, max_pattern_dist, format, "fwd"
            )
            # Extract rev UMI
            result_3p_rev_umi_dist, result_3p_rev_umi_seq, result_3p_rev_umi_qual, umi_start_rev = extract_umi(
                read_3p_seq, read_3p_qual, umi_rev, max_pattern_dist, format, "rev"
            )

            if not result_5p_fwd_umi_seq or not result_3p_rev_umi_seq:
                continue

            clipped_entry = clip_entry(
                entry, umi_start_fwd, umi_start_rev, adapter_length, format
            )
            
            n_both_umi += 1

            if format == "fasta":
                write_fasta(
                    clipped_entry,
                    strand,
                    result_5p_fwd_umi_dist,
                    result_3p_rev_umi_dist,
                    result_5p_fwd_umi_seq,
                    result_3p_rev_umi_seq,
                    out
                )

            if format == "fastq":
                write_fastq(
                    clipped_entry,
                    strand,
                    result_5p_fwd_umi_dist,
                    result_3p_rev_umi_dist,
                    result_5p_fwd_umi_seq,
                    result_3p_rev_umi_seq,
                    result_5p_fwd_umi_qual,
                    result_3p_rev_umi_qual,
                    out
                )

        if tsv:
            write_tsv(
                output_folder,
                output_file_name,
                max_pattern_dist,
                strand_stats,
                n_total,
                n_both_umi
            )


def reverse_quality(sign):
    n_total_qual_signs = 32 + 126
    return chr(n_total_qual_signs - sign)


def rev_comp_qual(seq_qual):
    return "".join(reverse_quality(ord(element)) for element in reversed(seq_qual))


def rev_comp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement.get(base, base) for base in reversed(seq))


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

    extract_umis(args)


if __name__ == "__main__":
    main()
