process SPLIT_CLUSTER {
    publishDir "${params.output}/${sample}/clustering/${type}/split_clusters", mode: 'copy'

  input:
    tuple val( sample ), val( target ), path( cluster )
    val (type )
    path umi_split_cluster_python
  output:
    tuple val( "${sample}" ), val( "${target}" ), path( "${cluster}_*" ), optional: true, emit:split_cluster_fastas
 
  script:
    def min_reads_per_cluster = type == "raw" ? "${params.min_reads_per_cluster}" : 1
  """
        python ${umi_split_cluster_python} \
         --min_reads_per_cluster $min_reads_per_cluster \
         --max_dist_umi ${params.max_dist_umi} \
         --cluster ${cluster} \
         -o . 
  """
}