import argparse
import pandas as pd
import os
import matplotlib.pyplot as plt

# Argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="Update summary cluster stats dynamically")
    parser.add_argument("--cluster-stat", required=True, help="new cluster stat TSV files")
    parser.add_argument("--sample", required=True, help="sample information")
    parser.add_argument("--target", required=True, help="target information")
    parser.add_argument("--task-index", type=int, required=True, help="task index of the current process")
    parser.add_argument("--current-summary", required=True, help="current-summary summary TSV file")
    parser.add_argument("--output-tsv", required=True, help="Output summary TSV file")
    parser.add_argument("--output-pdf", required=True, help="Output summary pdf file")
    return parser.parse_args()

# Parse cluster statistics
def parse_cluster_stat(cluster_stat_file, sample, target, task_index):
    cluster_stats = pd.read_csv(cluster_stat_file, sep='\t')
    n_cluster = cluster_stats[cluster_stats['cluster_written'] == 1].shape[0]
    cluster_stats_parsed = pd.DataFrame({
        'sample': [sample],
        'target': [target],
        'n_cluster': [n_cluster],
        'task_index': [task_index]
    })
    return cluster_stats_parsed

# Load and update summary
def load_and_update_summary(cluster_stats_parsed, sample, target, task_index, current_summary, output_tsv):
    if os.path.exists(current_summary):
        summary = pd.read_csv(current_summary, sep='\t')
        current_entry = summary[(summary['sample'] == sample) & (summary['target'] == target)]
        if not current_entry.empty and task_index > current_entry['task_index'].max():
            summary = summary[~((summary['sample'] == sample) & (summary['target'] == target))]
        summary = pd.concat([summary, cluster_stats_parsed], ignore_index=True)
    else:
        summary = cluster_stats_parsed
    summary.to_csv(output_tsv, sep='\t', index=False)
    return summary

# Generate histograms split by target
def plot_cluster_histogram(summary, output_pdf):
    targets = summary['target'].unique()

    fig, axes = plt.subplots(len(targets), 1, figsize=(10, 6 * len(targets)), squeeze=False)

    for ax, target in zip(axes.flatten(), targets):
        target_data = summary[summary['target'] == target]
        target_data.plot.barh(y='n_cluster', x='sample', ax=ax, legend=False)
        ax.set_xlabel('Number of Clusters')
        ax.set_ylabel('Sample')
        ax.set_title(f'Number of Found Clusters for Target: {target}')

    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close()

# Main execution
if __name__ == "__main__":
    args = parse_args()

    cluster_stats_parsed = parse_cluster_stat(args.cluster_stat, args.sample, args.target, args.task_index)
    summary = load_and_update_summary(cluster_stats_parsed, args.sample, args.target, args.task_index, args.current_summary, args.output_tsv)
    plot_cluster_histogram(summary, args.output_pdf)
