process SUMMARY_CLUSTER_STATS {    
    publishDir "${params.output}/cluster_stats/", mode: 'copy'
    
    input:
    path cluster_stats_files
    
    output:
    path "summary_cluster_stats.tsv"
    path "summary_cluster_report.pdf"
    
    script:
    """
    python generate_summary_report.py \
        --input ${cluster_stats_files} \
        --output-tsv summary_cluster_stats.tsv \
        --output-pdf summary_cluster_report.pdf 
    """
}