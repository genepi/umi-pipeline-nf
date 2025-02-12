import argparse
import pandas as pd
import matplotlib.pyplot as plt

# Argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Generate per-sample cluster report")
    parser.add_argument("--input", required=True, help="Input TSV file")
    parser.add_argument("--output-pdf", required=True, help="Output PDF file")
    parser.add_argument("--output-tsv", required=True, help="Output tsv file")
    parser.add_argument("--min-reads", type=int, required=True, help="Minimum reads per cluster")
    parser.add_argument("--threshold", type=int, required=True, help="Threshold buffer for near-threshold clusters")
    return parser.parse_args()

# Load data
def load_data(file_path):
    return pd.read_csv(file_path, sep='\t')

# Generate report
def generate_report(data, min_reads, threshold, output_pdf, output_tsv):
    clusters_written = data[data['cluster_written'] == 1]
    clusters_skipped = data[(data['cluster_written'] == 0)]
    near_threshold = clusters_skipped[(clusters_skipped['reads_found'] >= (min_reads - threshold))]
    n_bins_written = max(clusters_written["reads_found"]) - min_reads
    n_bins_skipped = min_reads - threshold
    
    fig, axes = plt.subplots(1, 2, figsize=(8, 10), sharex=False)
    
    # Plot 1: Clusters Written vs. Reads Found
    axes[0].hist(clusters_written['reads_found'], bins = n_bins_written, color='blue', align = "left", rwidth = 0.9)
    axes[0].set_title("Clusters Written vs Reads Found")
    axes[0].set_xlabel("Cluster Size")
    axes[0].set_ylabel("Count")
    
    # Plot 2: Near-Threshold Clusters
    axes[1].hist(near_threshold['reads_found'], bins = n_bins_skipped, color='red', alpha=0.7, align = "left", rwidth = 0.9)
    axes[1].set_title("Clusters Near Threshold")
    axes[1].set_xlabel("Cluster Size [{}-{}]".format(threshold, min_reads))
    axes[1].set_ylabel("Count")
    
    plt.tight_layout()
    plt.savefig(output_pdf)
    clusters_written.to_csv(output_tsv, sep='\t', index=False)

# Main execution
if __name__ == "__main__":
    args = parse_args()
    data = load_data(args.input)
    generate_report(data, args.min_reads, args.threshold, args.output_pdf, args.output_tsv)
