process POLISH_CLUSTER {
    tag "${sample}"

    input:
        tuple val( sample ), val( target ), path( smolecule_clusters_fastq )
        val ( type )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_clusters_fastq.baseName}_consensus.fastq" ), emit: consensus_fastq

    script:
    """
        medaka smolecule \
          --threads ${params.threads} \
          --length 50 \
          --depth 2 \
          --model ${params.medaka_model} \
          --method spoa . \
          --qualities \
          --save_features \
          ${smolecule_clusters_fastq} 2> smolecule.log
        
        mv consensus.fastq ${smolecule_clusters_fastq.baseName}_consensus.fastq
    """
}
