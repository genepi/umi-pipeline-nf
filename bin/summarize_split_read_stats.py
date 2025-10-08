#!/usr/bin/env python3
import argparse
import csv
import os
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Summarize stats from multiple stats files into one summary file.")
    parser.add_argument("stats_files", nargs="+", help="Input stats .tsv files (one or more)")
    parser.add_argument("--sample", required=True, help="Sample name to include in the summary")
    parser.add_argument("--target", required=True, help="Target name to include in the summary")
    parser.add_argument("-o", "--output", required=True, help="Output summary TSV file")
    return parser.parse_args()

def parse_stats_file(path):
    """
    Parse a stats file and return the counts and percentages as dictionaries.
    """
    with open(path, newline='') as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = [row for row in reader if row]  # Skip empty rows

    if len(rows) < 3:
        raise ValueError(f"Invalid stats file format: {path}")

    header = rows[0]
    count_row = rows[1]
    percent_row = rows[2]

    # Convert count and percent rows to dictionaries
    counts = {header[i]: int(count_row[i]) if count_row[i].isdigit() else count_row[i] for i in range(2, len(header))}
    percents = {header[i]: float(percent_row[i]) if percent_row[i].replace('.', '', 1).isdigit() else percent_row[i] for i in range(2, len(header))}

    return counts, percents

def summarize_filter_read_stats(stats_files):
    """
    Aggregate counts from multiple stats files and recalculate percentages for all columns.
    """
    aggregated_counts = defaultdict(int)

    # Aggregate counts from all files
    for stats_file in stats_files:
        counts, _ = parse_stats_file(stats_file)
        for key, value in counts.items():
            if isinstance(value, int):  # Only aggregate numeric columns
                aggregated_counts[key] += value

    # Calculate total counts for recalculating percentages
    total_counts = aggregated_counts["reads_found"]

    # Recalculate percentages for all columns
    aggregated_percents = {}
    for key, value in aggregated_counts.items():
        aggregated_percents[key] = (value / total_counts * 100) if total_counts > 0 else 0

    return aggregated_counts, aggregated_percents

def write_summary(sample, target, output, aggregated_counts, aggregated_percents):
    """
    Write the aggregated counts and percentages to the output file.
    """
    with open(output, "w", newline="") as outfh:
        writer = csv.writer(outfh, delimiter="\t")

        # Write the header
        writer.writerow(["sample", "target", "metric", "count", "percent"])

        # Write the aggregated counts and percentages
        for key in aggregated_counts:
            writer.writerow([sample, target, key, aggregated_counts[key], f"{aggregated_percents.get(key, 0):.2f}"])

    print(f"Summary written to {output}")

def main():
    args = parse_args()
    aggregated_counts, aggregated_percents = summarize_filter_read_stats(args.stats_files)
    write_summary(args.sample, args.target, args.output, aggregated_counts, aggregated_percents)

if __name__ == "__main__":
    main()