process SUMMARY_CLUSTER_STATS {   
    maxForks 1
    publishDir "${params.output}/cluster_stats/", mode: 'copy'  

    input:
    tuple val( sample ), val( target ), path( smolecule_cluster_stats )
    path cache_dir
    path umi_cluster_stats_summary

    output:
    path "summary_cluster_stats.tsv"
    path "summary_cluster_report.pdf"
    
    script:
    def current_summary = "${cache_dir}/summary_cluster_stats.tsv"
    """
        if [[ "${task.index}" == "1" ]]; then
            mkdir -p ${cache_dir}
        fi
        
        python ${umi_cluster_stats_summary} \
            --cluster-stat ${smolecule_cluster_stats} \
            --sample ${sample} \
            --task-index ${task.index} \
            --current-summary $current_summary \
            --output-tsv summary_cluster_stats.tsv \
            --output-pdf summary_cluster_report.pdf 

        cp summary_cluster_stats.tsv ${cache_dir}
    """
}