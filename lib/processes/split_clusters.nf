vsearch_dir="vsearch_clusters"

process SPLIT_CLUSTER {
    publishDir "${params.output}/${sample}/clustering/${type}/${vsearch_dir}", mode: 'copy'

  input:
    tuple val( sample ), val( target ), path( cluster )
    val (type )
    path umi_split_cluster_python
  output:
    tuple val( "${sample}" ), val( "${target}" ), path( "cluster*" ), emit:cluster_fastas
 
  script:
  """
        python ${umi_split_cluster_python} \
         --min_reads_per_cluster ${params.min_reads_per_cluster} \
         --max_dist_umi ${params.max_dist_umi} \
         --cluster ${cluster} \
         -o .
  """
}