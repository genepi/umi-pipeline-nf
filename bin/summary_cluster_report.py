import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt

# Argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Update summary cluster stats dynamically")
    parser.add_argument("--cluster-stat", required=True, help="new cluster stat TSV files")
    parser.add_argument("--sample", required=True, help="sample information")
    parser.add_argument("--task-index", type=int, required=True, help="task index of the current process")
    parser.add_argument("--current-summary", required=True, help="current-summary summary TSV file")
    parser.add_argument("--output-tsv", required=True, help="Output summary TSV file")
    parser.add_argument("--output-pdf", required=True, help="Output summary pdf file")
    return parser.parse_args()

def parse_cluster_stat(cluster_stat_file, sample, task_index):
    cluster_stats = pd.read_csv(cluster_stat_file, sep='\t')
    n_cluster = cluster_stats[cluster_stats['cluster_written'] == 1].shape[0]
    cluster_stats_parsed = pd.DataFrame({'sample': [sample], 'n_cluster': [n_cluster], 'task_index': [task_index]})

    return cluster_stats_parsed, sample


# Load and combine data without duplicating samples
def load_and_update_summary(cluster_stats_parsed, sample, task_index, current_summary, output_tsv):
    # Load existing summary if it exists
    if os.path.exists(current_summary):
        summary = pd.read_csv(current_summary, sep='\t')
        # Remove old entries for this sample if present
        current_entry = summary[summary['sample'] == sample]

        if not current_entry.empty:
            if task_index > current_entry['task_index'].max():
                summary = summary[summary['sample'] != sample]

        summary = summary[summary['sample'] != sample]
        summary = pd.concat([summary,cluster_stats_parsed], ignore_index = True)
    else:
        summary = cluster_stats_parsed
        
    summary.to_csv(output_tsv, sep='\t', index=False)
    
    return summary

# Generate horizontal histogram for the 'cluster_written' column
def plot_cluster_histogram(summary, output_pdf):
    plt.figure(figsize=(10, 6))
    
    summary.plot.barh(y = 'n_cluster', x = 'sample')
    plt.xlabel('Number of cluster')
    plt.ylabel('Sample')
    plt.title('Number of Found Clusters for Each Sample')
    plt.legend(title='Sample')
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close()


# Main execution
if __name__ == "__main__":
    args = parse_args()
    
    sample = args.sample
    task_index = args.task_index
    cluster_stats_parsed, sample = parse_cluster_stat(args.cluster_stat, sample, task_index)
    summary = load_and_update_summary(cluster_stats_parsed, sample, task_index, args.current_summary, args.output_tsv)
    plot_cluster_histogram(summary, args.output_pdf)
