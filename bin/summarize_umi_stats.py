#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Summarize UMI stats from multiple stats files.")
    parser.add_argument("stats_files", nargs="+", help="Input stats .tsv files (one or more)")
    parser.add_argument("--sample", required=True, help="Sample name to include in the summary")
    parser.add_argument("--target", required=True, help="Target name to include in the summary")
    parser.add_argument("-o", "--output", required=True, help="Output summary TSV file")
    return parser.parse_args()

def parse_stats_file(path):
    """
    Parse a stats file and return the relevant columns as a dictionary.
    """
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = [row for row in reader]

    if not rows:
        raise ValueError(f"Stats file is empty: {path}")

    # Extract only the relevant columns
    row = rows[0]
    return {
        "detected_forward_strands": int(row["detected_forward_strands"]),
        "detected_reverse_strands": int(row["detected_reverse_strands"]),
        "total_reads": int(row["total_reads"]),
        "included_reads": int(row["included_reads"]),
    }

def summarize_stats(stats_files):
    """
    Aggregate stats from multiple files and calculate percentages.
    """
    aggregated_stats = defaultdict(int)

    # Aggregate stats from all files
    for stats_file in stats_files:
        try:
            stats = parse_stats_file(stats_file)
            for key, value in stats.items():
                aggregated_stats[key] += value
        except ValueError as e:
            print(f"Warning: {e}")

    # Calculate percentages relative to total_reads
    total_reads = aggregated_stats["total_reads"]
    percentages = {
        key: (value / total_reads * 100) if total_reads > 0 else 0
        for key, value in aggregated_stats.items()
    }

    return aggregated_stats, percentages

def write_summary(output, sample, target, aggregated_stats, percentages):
    """
    Write the aggregated stats and percentages to the output file.
    """
    with open(output, "w", newline="") as outfh:
        writer = csv.writer(outfh, delimiter="\t")
        writer.writerow(["sample", "target", "metric", "count", "percent"])
        for key, value in aggregated_stats.items():
            writer.writerow([sample, target, key, value, f"{percentages[key]:.2f}"])

def main():
    args = parse_args()
    aggregated_stats, percentages = summarize_stats(args.stats_files)
    write_summary(args.output, args.sample, args.target, aggregated_stats, percentages)

if __name__ == "__main__":
    main()