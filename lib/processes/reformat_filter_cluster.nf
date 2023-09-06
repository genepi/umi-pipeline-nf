process REFORMAT_FILTER_CLUSTER {
    publishDir "${params.output}/${sample}/stats/${type}", pattern: "*.tsv", mode: 'copy'
    publishDir "${params.output}/${sample}/clustering/${type}", pattern: "smolecule*", mode: 'copy'
      
    input:
        tuple val( sample ), val( target ), path( cluster_fastq )
        val( type )
        path umi_parse_clusters_python

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "smolecule*"), emit: smolecule_cluster_fastq
        path( "*.tsv" )

    script:
        def balance_strands = params.balance_strands ? "--balance_strands" : ""
        def write_report = params.write_reports ? "--tsv" : ""
    """
        python ${umi_parse_clusters_python} \
          --filter_strategy ${params.filter_strategy_clusters} \
          --min_reads_per_clusters ${params.min_reads_per_cluster} \
          --max_reads_per_clusters ${params.max_reads_per_cluster} \
          --cluster ${cluster_fastq} \
          --output_format ${params.output_format} \
          $balance_strands \
          $write_report \
          -o .
    """
}