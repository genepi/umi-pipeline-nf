import argparse
import logging
import os
import sys
import time

import threading as t

import pysam
import edlib


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
        "--max_dist_umi",
        dest="MAX_EDIT_DIST",
        type=int,
        default=2,
        help="Max distance of UMIs per cluster. Used to split cluster into subclusters",
    )
    
    parser.add_argument(
        "-o", "--output", dest="OUTPUT", required=True, help="Output folder"
    )

    parser.add_argument(
        "--tsv", dest="TSV", action="store_true", help="write TSV output file"
    )

    parser.add_argument(
        "--clusters", dest="CLUSTERS", nargs= "+", required=True, type=str, help="cluster fastxs"
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
    return cluster.split("cluster")[-1]

def get_read_seq(read):
    return read.name.split(";seq=")[1].split(";")[0]

def get_read_qual(read):
    return read.name.split(";qual=")[1]

def get_read_strand(read):
    return read.name.split(";strand=")[1].split(";")[0]

def get_read_name(read):
    return read.name.split(";")[0]

def get_read_mean_qual(read):
    qual = get_read_qual(read)
    return get_mean_qual(qual)

def get_mean_qual(qual):
    return sum(map(lambda char: ord(char), qual)) / len(qual)

def get_reads(cluster):
    residual_reads = []
    n_residual_reads = 0
    with pysam.FastxFile(cluster) as reads:
        for read in reads:
            residual_reads.append(read)
            n_residual_reads += 1
    return residual_reads, n_residual_reads

def get_split_cluster(reads, max_edit_dist):
    subcluster = []
    distant_reads = []
    parent = reads[0].sequence
    
    for read in reads:
        # calculate edit distance between parent and all other reads in the cluster
        result = edlib.align(
            parent,
            read.sequence,
            mode="HW",
            k=max_edit_dist
        )
        if result["editDistance"] == -1:
            distant_reads.append(read)
        else:
            subcluster.append(read)
    return subcluster, distant_reads                    


def get_split_reads(reads):
    reads_fwd = []
    reads_rev = []
    
    for read in reads:
        strand = get_read_strand(read)
        if strand == "+":
            reads_fwd.append(read)
        else:
            reads_rev.append(read)
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

def write_smolecule(cluster_id, reads, smolecule_file, format):
    with open(smolecule_file, "w") as out_f:
        for n, read in enumerate(reads):
            seq = get_read_seq(read)
            read_name = "{}_{}".format(cluster_id, n)
            if format == "fastq":
                qual = get_read_qual(read)
                write_fastq_read(read_name, seq, qual, out_f)
            else:
                write_fasta_read(read_name, seq, out_f)


def write_fastq_read(read_name, read_seq, qual, out_f):
    print("@{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)
    print("+", file=out_f)
    print("{}".format(qual), file=out_f)


def write_fasta_read(read_name, read_seq, out_f):
    print(">{}".format(read_name), file=out_f)
    print("{}".format(read_seq), file=out_f)

def write_subcluster(subcluster, subcluster_name):
    with open(subcluster_name, "w") as out_f:
        for read in subcluster:
            print(">{}".format(read.name), file=out_f)
            print("{}".format(read.sequence), file=out_f)
            
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

def parse_cluster(min_reads, max_reads, filter, format, cluster, output_folder, balance_strands, tsv, max_edit_dist, stats_out_filename):  
    reads_found = 0
    reads_found = 0
    reads_written_fwd = 0
    reads_written_rev = 0
    cluster_written = 0

    reads_fwd = []
    reads_rev = []
    cluster_id = get_cluster_id(cluster)
    
    n_subcluster = 0   
        
    smolecule_file = os.path.join(
        output_folder, "smolecule{}.{}".format(cluster_id, format))
    
    residual_reads, n_residual_reads = get_reads(cluster)
    
    while n_residual_reads >= min_reads:
            cluster_id_subcluster = "{}_{}".format(cluster_id, n_subcluster)
            
            subcluster, residual_reads = get_split_cluster(residual_reads, max_edit_dist)
            n_residual_reads = len(residual_reads)
            write_subcluster(
                subcluster,
                os.path.join(output_folder, "{}_{}".format(cluster, n_subcluster))
                )
            n_subcluster += 1
        
            reads_fwd, reads_rev = get_split_reads(subcluster)
            n_fwd = len(reads_fwd)
            n_rev = len(reads_rev)
            reads_found = n_fwd + n_rev

            if filter == "quality":
                reads_fwd = get_sorted_reads(reads_fwd)
                reads_rev = get_sorted_reads(reads_rev)

            reads_fwd, reads_rev, write_cluster, reads_skipped_fwd, reads_skipped_rev = get_filtered_reads(
                reads_fwd, reads_rev, reads_found, n_fwd, n_rev, min_reads, max_reads, balance_strands)

            if write_cluster:
                cluster_written = 1
                reads_written_fwd = len(reads_fwd)
                reads_written_rev = len(reads_rev)
                
                reads = reads_fwd + reads_rev
                write_smolecule(cluster_id_subcluster, reads, smolecule_file, format)
            else:
                cluster_written = 0

            if tsv:  
                write_tsv_line(stats_out_filename, cluster_id_subcluster, cluster_written, reads_found, n_fwd,
                n_rev, reads_written_fwd, reads_written_rev, reads_skipped_fwd, reads_skipped_rev)

def parse_cluster_wrapper(args):
    min_reads = args.MIN_CLUSTER_READS
    max_reads = args.MAX_CLUSTER_READS
    filter = args.FILTER
    format = args.OUT_FORMAT
    clusters = args.CLUSTERS
    output_folder = args.OUTPUT
    balance_strands = args.BAL_STRANDS
    tsv = args.TSV
    max_edit_dist = args.MAX_EDIT_DIST
    
    stats_out_filename = "split_cluster_stats"
    
    if tsv:
        stats_out_filename = os.path.join(
            output_folder, "{}.tsv".format(stats_out_filename))
        write_tsv_line(stats_out_filename, "cluster_id", "cluster_written", "reads_found", "reads_found_fwd",
                    "reads_found_rev", "reads_written_fwd", "reads_written_rev", "reads_skipped_fwd", "reads_skipped_rev")

    for cluster in clusters:
        parse_cluster_thread = t.Thread(target=parse_cluster, args=(
            min_reads, max_reads, filter, format, cluster, output_folder, balance_strands, tsv, max_edit_dist, stats_out_filename
        ))
        parse_cluster_thread.start()
        #parse_cluster(min_reads, max_reads, filter, format, cluster, output_folder, balance_strands, tsv, max_edit_dist, stats_out_filename)

    
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
        parse_cluster_wrapper(args)
    except RuntimeError as e:
        logging.error(e)


if __name__ == "__main__":
    main()
