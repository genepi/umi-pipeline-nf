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
        "--balance_strands",
        dest="BAL_STRANDS",
        action="store_true",
        help="Balance strands in clusters",
    )

    parser.add_argument(
        "--filter_strategy",
        dest="FILTER",
        type=str,
        default="random",
        help="[ random | quality ] Choose strategy to remove reads. Only if --balance_strands is set"
        # TODO Check vsearch: Reads might be sorted by length, so filtering would not be random!
    )

    parser.add_argument(
        "--min_reads_per_clusters",
        dest="MIN_CLUSTER_READS",
        type=int,
        default=20,
        help="Reads per cluster. Clusters with less reads will be discarded, clusters with more will be downsampled. 50% must be forward and 50% reverse reads",
    )

    parser.add_argument(
        "--max_reads_per_clusters",
        dest="MAX_CLUSTER_READS",
        type=int,
        default=60,
        help="Reads per cluster. Clusters with less reads will be discarded, clusters with more will be downsampled. 50% must be forward and 50% reverse reads",
    )

    parser.add_argument(
        "-o", "--output", dest="OUTPUT", required=True, help="Output folder"
    )

    parser.add_argument(
        "--tsv", dest="TSV", action="store_true", help="write TSV output file"
    )

    parser.add_argument(
        "--vsearch_consensus", dest="VSEARCH_CONSENSUS", required=True, type=str, help="VSearch consensus FASTX"
    )

    parser.add_argument(
        "--vsearch_folder", dest="VSEARCH_FOLDER", required=True, type=str, help="VSearch cluster folder"
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


def get_cluster_id(cluster):
    return cluster.name.split(";")[-1].split("=")[1]

def get_cluster_seq_n(cluster):
    return int(cluster.name.split(";")[-2].split("=")[1])

def get_read_seq(entry):
    return entry.name.split(";seq=")[1].split(";")[0]

def get_read_qual(entry):
    return entry.name.split(";qual=")[1]

def get_read_strand(entry):
    return entry.name.split("strand=")[1].split(";")[0]

def get_read_name(entry):
    return entry.name.split(";")[0]

def get_read_mean_qual(entry):
    qual = get_read_qual(entry)
    return get_mean_qual(qual)

def get_read_qual(entry):
    return entry.name.split("qual=")[1].split(";seqs=")[0]

def get_mean_qual(qual):
    return sum(map(lambda char: ord(char), qual)) / len(qual)

def get_split_reads(cluster_fasta_umi):
    reads_fwd = []
    reads_rev = []
    with pysam.FastxFile(cluster_fasta_umi) as fh:
        for entry in fh:
            strand = get_read_strand(entry)
            
            if strand == "+":
                reads_fwd.append(entry)
            else:
                reads_rev.append(entry)
    return reads_fwd, reads_rev

def get_sorted_reads(reads):
    return sorted(reads, key=get_read_mean_qual, reverse=True)

def get_filter_parameters(n_fwd, n_rev, min_reads, max_reads, balance_strands):
    if balance_strands:
        min_fwd = min_rev = min_reads // 2
        max_reads = min(n_fwd * 2, n_rev * 2, max_reads)
        max_fwd = max_rev = max_reads // 2
    else:
        min_fwd = 0
        min_rev = 0

        if n_fwd > n_rev:
            max_rev = min(n_rev, max_reads // 2)
            max_fwd = min(max_reads - max_rev, n_fwd)
        else:
            max_fwd = min(n_fwd, max_reads // 2)
            max_rev = min(max_reads - max_fwd, n_rev)

    return min_fwd, min_rev, max_fwd, max_rev


def get_filtered_reads(reads_fwd, reads_rev, reads_found, n_fwd, n_rev, min_reads, max_reads, balance_strands):
    skipped_fwd = 0
    skipped_rev = 0
    write_cluster = True

    min_fwd, min_rev, max_fwd, max_rev = get_filter_parameters(
        n_fwd, n_rev, min_reads, max_reads, balance_strands)

    if n_fwd < min_fwd or n_rev < min_rev or reads_found < min_reads:
        write_cluster = False
        skipped_fwd = n_fwd
        skipped_rev = n_rev
        return reads_fwd, reads_rev, write_cluster, skipped_fwd, skipped_rev
    if n_fwd > max_fwd:
        reads_fwd = reads_fwd[:max_fwd]
        skipped_fwd = n_fwd - max_fwd
    if n_rev > max_rev:
        reads_rev = reads_rev[:max_rev]
        skipped_rev = n_rev - max_rev

    return reads_fwd, reads_rev, write_cluster, skipped_fwd, skipped_rev

def polish_cluster(
    cluster_id,
    vsearch_folder,
    output_folder,
    min_reads,
    max_reads,
    filter,
    stats_out_filename,
    tsv,
    balance_strands,
    smolecule_out,
    format
):
    reads_found = 0
    reads_written = 0
    reads_written_fwd = 0
    reads_written_rev = 0

    reads_skipped = 0
    cluster_written = 0

    reads_fwd = []
    reads_rev = []

    cluster_fasta_umi = os.path.join(
        vsearch_folder, "cluster{}".format(cluster_id))
    parsed_cluster = os.path.join(
        output_folder, "cluster{}.{}".format(cluster_id, format))
    smolecule_file = os.path.join(
        output_folder, "smolecule{}.{}".format(cluster_id, format))

    parse_cluster(cluster_fasta_umi, parsed_cluster, format)

    reads_fwd, reads_rev = get_split_reads(cluster_fasta_umi)
    n_fwd = len(reads_fwd)
    n_rev = len(reads_rev)
    reads_found = n_fwd + n_rev

    # Fail fast if no stat file must be written
    if not tsv:
        if reads_found < min_reads:
            cluster_written = 0
            reads_written = 0
            reads_skipped = reads_found
            return cluster_written, reads_found, reads_skipped, reads_written
    else:
        if filter == "quality":
            reads_fwd = get_sorted_reads(reads_fwd)
            reads_rev = get_sorted_reads(reads_rev)

    #logging.info( isinstance(reads_fwd, list), isinstance(reads_rev, list), cluster_id)
    reads_fwd, reads_rev, write_cluster, reads_skipped_fwd, reads_skipped_rev = get_filtered_reads(
        reads_fwd, reads_rev, reads_found, n_fwd, n_rev, min_reads, max_reads, balance_strands)

    if write_cluster:
        cluster_written = 1
        reads_written_fwd = len(reads_fwd)
        reads_written_rev = len(reads_rev)
        
        reads = reads_fwd + reads_rev
        write_smolecule(cluster_id, reads, smolecule_file, format)
    else:
        cluster_written = 0

    if tsv:
        write_tsv_line(stats_out_filename, cluster_id, cluster_written, reads_found, n_fwd,
                       n_rev, reads_written_fwd, reads_written_rev, reads_skipped_fwd, reads_skipped_rev)

    return (
        cluster_written,
        reads_found,
        reads_skipped_fwd + reads_skipped_rev,
        reads_written_rev + reads_written_fwd
    )


def parse_cluster(cluster_fasta_umi, parsed_cluster, format):
    with open(parsed_cluster, "w") as out_f:
        with pysam.FastxFile(cluster_fasta_umi) as fh:
            for entry in fh:
                read_name = get_read_name(entry)
                seq = get_read_seq(entry)
                if format == "fastq":
                    qual = get_read_qual(entry)
                    write_fastq_entry(read_name, seq, qual, out_f)
                else:
                    write_fasta_entry(read_name, seq, out_f)


def write_smolecule(cluster_id, reads, smolecule_file, format):
    with open(smolecule_file, "w") as out_f:
        for n, entry in enumerate(reads):
            seq = get_read_seq(entry)
            read_name = "{}_{}".format(cluster_id, n)
            if format == "fastq":
                qual = get_read_qual(entry)
                write_fastq_entry(read_name, seq, qual, out_f)
            else:
                write_fasta_entry(read_name, seq, out_f)


def write_fastq_entry(read_name, read_seq, qual, out_f):
    print("@{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)
    print("+", file=out_f)
    print("{}".format(qual), file=out_f)


def write_fasta_entry(read_name, read_seq, out_f):
    print(">{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)


def write_tsv_line(stats_out_filename, cluster_id, cluster_written, reads_found, n_fwd, n_rev, reads_written_fwd, reads_written_rev, reads_skipped_fwd, reads_skipped_rev):
    with open(stats_out_filename, "a") as out_f:
        print(cluster_id,
              cluster_written,
              reads_found,
              n_fwd,
              n_rev,
              reads_written_fwd,
              reads_written_rev,
              reads_skipped_fwd,
              reads_skipped_rev,
              sep="\t",
              file=out_f,
              )

def parse_clusters(args):
    cluster_filename = args.VSEARCH_CONSENSUS
    min_read_per_cluster = args.MIN_CLUSTER_READS
    max_read_per_cluster = args.MAX_CLUSTER_READS
    filter = args.FILTER
    format = args.OUT_FORMAT
    vsearch_folder = args.VSEARCH_FOLDER
    output = args.OUTPUT
    balance_strands = args.BAL_STRANDS
    tsv = args.TSV

    smolecule_filename = "smolecule_clusters"
    stats_out_filename = "vsearch_cluster_stats"
    n_clusters = 0
    n_written = 0
    reads_found = 0
    reads_skipped = 0
    reads_written = 0

    if tsv:
        stats_out_filename = os.path.join(
            output, "{}.tsv".format(stats_out_filename))
        write_tsv_line(stats_out_filename, "cluster_id", "cluster_written", "reads_found", "reads_found_fwd",
                       "reads_found_rev", "reads_written_fwd", "reads_written_rev", "reads_skipped_fwd", "reads_skipped_rev")

    with pysam.FastxFile(cluster_filename) as fh:
        for cluster in fh:

            # cons_umi = None cluster.sequence.replace("T", "")
            cluster_id = get_cluster_id(cluster)

            a, b, c, d = polish_cluster(
                cluster_id,
                vsearch_folder,
                output,
                min_read_per_cluster,
                max_read_per_cluster,
                filter,
                stats_out_filename,
                tsv,
                balance_strands,
                smolecule_filename,
                format
            )
            n_clusters += 1
            n_written += a
            reads_found += b
            reads_skipped += c
            reads_written += d

    if n_written == 0 or reads_found == 0:
        raise RuntimeError(
            "ERROR - did not find a single cluster passing the min_read threashold!"
        )

    logging.info(
        "Clusters: {}% written ({})".format(
            n_written * 100.0 / n_clusters, n_written
        )
    )
    logging.info(
        "Reads: {} found".format(
            reads_found
        )
    )
    
    logging.info(
        "Reads: {} removed ({}%)".format(
            reads_skipped, reads_skipped * 100.0 // reads_found
        )
    )
    logging.info("Reads: {} written ({}%)".format(
        reads_written, reads_written * 100.0 // reads_found))


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

    try:
        parse_clusters(args)
    except RuntimeError as e:
        logging.error(e)


if __name__ == "__main__":
    main()
