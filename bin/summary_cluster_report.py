import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt

# Argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Update summary cluster stats dynamically")
    parser.add_argument("--cluster-stat", required=True, help="new cluster stat TSV files")
    parser.add_argument("--current-summary", required=True, help="current-summary summary TSV file")
    parser.add_argument("--output-tsv", required=True, help="Output summary TSV file")
    parser.add_argument("--output-pdf", required=True, help="Output summary pdf file")
    return parser.parse_args()

def parse_cluster_stat(cluster_stat_file):
    cluster_stats = pd.read_csv(cluster_stat_file, sep='\t')

    sample = os.path.basename(cluster_stat_file).split('_')[0]
    n_cluster = cluster_stats[cluster_stats['cluster_written'] == 1].shape[0]
    cluster_stats_parsed = pd.DataFrame({'sample': [sample], 'n_cluster': [n_cluster]})

    return cluster_stats_parsed


# Load and combine data without duplicating samples
def load_and_update_summary(cluster_stats_parsed, current_summary, output_tsv):
    # Load existing summary if it exists
    if os.path.exists(current_summary):
        summary = pd.read_csv(current_summary, sep='\t')
        # Remove old entries for this sample if present
        summary = summary[summary['sample'] != cluster_stats_parsed['sample']]
        summary = pd.concat([summary,cluster_stats_parsed], ignore_index = True)
    else:
        summary = cluster_stats_parsed
        
    print(summary)
    # Append new data and save
    summary.to_csv(output_tsv, sep='\t', index=False)
    
    
    return summary

# Generate horizontal histogram for the 'cluster_written' column
def plot_cluster_histogram(summary, output_pdf):
    plt.figure(figsize=(10, 6))
    
    summary.plot.barh(y = 'n_cluster', x = 'sample')
    plt.xlabel('Number of cluster')
    plt.ylabel('Sample')
    plt.title('Horizontal Histogram of Cluster Written for Each Sample')
    plt.legend(title='Sample')
    plt.tight_layout()
    plt.savefig(output_pdf)  # Save the plot to a PDF
    plt.close()


# Main execution
if __name__ == "__main__":
    args = parse_args()
    cluster_stats_parsed = parse_cluster_stat(args.cluster_stat)
    summary = load_and_update_summary(cluster_stats_parsed, args.current_summary, args.output_tsv)
    plot_cluster_histogram(summary, args.output_pdf)
