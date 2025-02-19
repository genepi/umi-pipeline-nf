process SUMMARY_CLUSTER_STATS {   
    maxForks 1
    publishDir "${params.output}/cluster_stats/", mode: 'copy'  
    
    input:
    tuple val(sample), val(type), path(smolecule_cluster_stats)
    path umi_cluster_stats_summary

    output:
    path "summary_cluster_stats.tsv"
    path "summary_cluster_report.pdf"
    
    script:
    def current_summary = "${launchDir}/${params.output}/cluster_stats/summary_cluster_stats.tsv"
    """
    python ${umi_cluster_stats_summary} \
        --cluster-stat ${smolecule_cluster_stats} \
        --sample ${sample} \
        --task-index ${task.index} \
        --current-summary $current_summary \
        --output-tsv summary_cluster_stats.tsv \
        --output-pdf summary_cluster_report.pdf 
    """
}