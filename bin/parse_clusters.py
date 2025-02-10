#!/usr/bin/env python
"""
This is a modified version of the code present in:
https://github.com/nanoporetech/pipeline-umi-amplicon/blob/master/lib/umi_amplicon_tools/parse_clusters.py

In this version, we check that reads within each input cluster are “similar” by
clustering them based on their pairwise edit distances. Two reads are connected
if their edit distance is <= max_edit_dist. Connected components (using networkx)
are then used as subclusters.
"""

import argparse
import logging
import os
import sys
import time

import threading as t

import pysam
import edlib
import networkx as nx
from itertools import combinations


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
        help="Max edit distance allowed between reads in a subcluster",
    )
    
    parser.add_argument(
        "-o", "--output", dest="OUTPUT", required=True, help="Output folder"
    )

    parser.add_argument(
        "--tsv", dest="TSV", action="store_true", help="write TSV output file"
    )

    parser.add_argument(
        "--clusters", dest="CLUSTERS", nargs="+", required=True, type=str, help="cluster fastx files"
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
    """
    Reads the fastx file for a given cluster and returns a list of reads.
    """
    residual_reads = []
    n_residual_reads = 0
    with pysam.FastxFile(cluster) as reads:
        for read in reads:
            residual_reads.append(read)
            n_residual_reads += 1
    return residual_reads, n_residual_reads


def cluster_reads(reads, max_edit_dist):
    """
    Groups reads into subclusters using a graph approach.
    Each read is a node and an edge is added if the edit distance between two reads
    is <= max_edit_dist. Returns a list of subclusters (each subcluster is a list of reads).
    """
    G = nx.Graph()
    # add all reads as nodes
    for i, read in enumerate(reads):
        G.add_node(i)
    # add an edge for every pair that is within the edit distance threshold
    for i, j in combinations(range(len(reads)), 2):
        result = edlib.align(reads[i].sequence, reads[j].sequence, mode="HW", k=max_edit_dist)
        if result["editDistance"] != -1:
            G.add_edge(i, j)
    # Each connected component represents a subcluster
    subclusters = []
    for component in nx.connected_components(G):
        subcluster = [reads[i] for i in component]
        subclusters.append(subcluster)
    return subclusters


def get_split_reads(reads):
    """
    Splits reads into forward and reverse based on the strand.
    """
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
    """
    Sort reads based on the average quality (highest quality first).
    """
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
    """
    Depending on filtering parameters, returns (possibly trimmed) forward and reverse read lists,
    a flag whether to write the cluster, and counts of skipped reads.
    """
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
    """
    Write the combined (filtered) reads for a subcluster as either FASTA or FASTQ.
    """
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


def write_subcluster(subcluster, subcluster_filename):
    """
    Writes the reads in a subcluster to a file (in FASTA format).
    """
    with open(subcluster_filename, "w") as out_f:
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
              file=out_f)


def parse_cluster(min_reads, max_reads, filter, format, cluster, output_folder, balance_strands, tsv, max_edit_dist, stats_out_filename):
    """
    For each input cluster file, read the sequences, break them into subclusters such that
    sequences in each subcluster are similar (i.e. within max_edit_dist of at least one other read),
    and then perform filtering and write the results.
    """
    residual_reads, n_residual_reads = get_reads(cluster)
    cluster_id = get_cluster_id(cluster)
    
    # Cluster reads into subclusters based on pairwise edit distance
    subclusters = cluster_reads(residual_reads, max_edit_dist)
    
    for n_subcluster, subcluster in enumerate(subclusters):
        # Write the subcluster reads (for reference/debugging)
        subcluster_file = os.path.join(output_folder, "{}_subcluster_{}".format(cluster_id, n_subcluster))
        write_subcluster(subcluster, subcluster_file)
        
        reads_found = len(subcluster)
        reads_fwd, reads_rev = get_split_reads(subcluster)
        n_fwd = len(reads_fwd)
        n_rev = len(reads_rev)

        if filter == "quality":
            reads_fwd = get_sorted_reads(reads_fwd)
            reads_rev = get_sorted_reads(reads_rev)

        reads_fwd, reads_rev, write_cluster, reads_skipped_fwd, reads_skipped_rev = get_filtered_reads(
            reads_fwd, reads_rev, reads_found, n_fwd, n_rev, min_reads, max_reads, balance_strands
        )

        cluster_written = 0
        reads_written_fwd = 0
        reads_written_rev = 0

        cluster_id_subcluster = "{}_sub{}".format(cluster_id, n_subcluster)
        smolecule_file = os.path.join(
            output_folder, "smolecule{}.{}".format(cluster_id_subcluster, format)
        )

        if write_cluster:
            cluster_written = 1
            reads_written_fwd = len(reads_fwd)
            reads_written_rev = len(reads_rev)
            reads = reads_fwd + reads_rev
            write_smolecule(cluster_id_subcluster, reads, smolecule_file, format)

        if tsv:
            write_tsv_line(
                stats_out_filename,
                cluster_id_subcluster,
                cluster_written,
                reads_found,
                n_fwd,
                n_rev,
                reads_written_fwd,
                reads_written_rev,
                reads_skipped_fwd,
                reads_skipped_rev,
            )


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
            output_folder, "{}.tsv".format(stats_out_filename)
        )
        # write header
        write_tsv_line(
            stats_out_filename,
            "cluster_id",
            "cluster_written",
            "reads_found",
            "reads_found_fwd",
            "reads_found_rev",
            "reads_written_fwd",
            "reads_written_rev",
            "reads_skipped_fwd",
            "reads_skipped_rev",
        )

    for cluster in clusters:
        # You can run each cluster in its own thread if desired:
        # parse_cluster_thread = t.Thread(target=parse_cluster, args=(
        #     min_reads, max_reads, filter, format, cluster, output_folder, balance_strands, tsv, max_edit_dist, stats_out_filename
        # ))
        # parse_cluster_thread.start()
        parse_cluster(
            min_reads,
            max_reads,
            filter,
            format,
            cluster,
            output_folder,
            balance_strands,
            tsv,
            max_edit_dist,
            stats_out_filename,
        )


def main(argv=sys.argv[1:]):
    """
    Basic command line interface.
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
