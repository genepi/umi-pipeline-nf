cluster_stats_file="cluster_stats.tsv"
process MERGE_CLUSTER_STATS {
    publishDir "${params.output}/${sample}/stats/${type}", pattern: "*tsv", mode: 'copy'

  input:
    tuple val( sample ), val( target ), path( cluster_stats )
    val( type )
  output:
    tuple val( "${sample}" ), val( "${target}" ), path ( "*tsv" )
        
  script:
  """
    echo -e "cluster_id\tcluster_written\treads_found\treads_found_fwd\treads_found_rev\treads_written_fwd\treads_written_rev\treads_skipped_fwd\treads_skipped_rev" \
    > ${cluster_stats_file}
    cat ${cluster_stats} >> ${cluster_stats_file}
  """
}