import argparse
import logging
import os
import sys

import pyfastx
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
        "--min_reads_per_clusters",
        dest="MIN_CLUSTER_READS",
        type=int,
        default=20,
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
        "--cluster", dest="CLUSTER", required=True, type=str, help="cluster fastx"
    )

    args = parser.parse_args(argv)

    return args

def get_split_cluster(reads, max_edit_dist):
    subcluster = []
    distant_reads = []
    parent = reads[0].seq
    
    for read in reads:
        # calculate edit distance between parent and all other reads in the cluster
        result = edlib.align(
            parent,
            read.seq,
            mode="HW",
            k=max_edit_dist
        )
        if result["editDistance"] == -1:
            distant_reads.append(read)
        else:
            subcluster.append(read)
    return subcluster, distant_reads                    

def write_subcluster(subcluster, subcluster_name):
    with open(subcluster_name, "w") as out_f:
        for read in subcluster:
            print(">{}".format(read.name), file=out_f)
            print("{}".format(read.seq), file=out_f)

def parse_clusters(args):
    min_reads = args.MIN_CLUSTER_READS
    cluster = args.CLUSTER
    output_folder = args.OUTPUT
    max_edit_dist = args.MAX_EDIT_DIST
    
    n_subcluster = 0       
    residual_reads = pyfastx.Fasta(cluster)
    
    n_residual_reads = len(residual_reads)
    while n_residual_reads > min_reads:
        subcluster, residual_reads = get_split_cluster(residual_reads, max_edit_dist)
        n_residual_reads = len(residual_reads)
        if len(subcluster) > min_reads:
            write_subcluster(
                subcluster,
                os.path.join(output_folder, "{}_{}".format(cluster, n_subcluster))
                )
            n_subcluster += 1
    
    #if n_subcluster == 0:     
    #    # Emit empty file to satisfy nextflow downstrean process 
    #    open(os.path.join(output_folder, "{}_empty".format(cluster)), "w")

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
