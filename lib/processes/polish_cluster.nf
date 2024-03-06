cpus=2
process POLISH_CLUSTER {
    cpus "${cpus}"
    tag "${sample}"

    input:
        tuple val( sample ), val( target ), path( smolecule_clusters_fastq )
        val ( type )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_clusters_fastq.baseName}_consensus.fastq" ), emit: consensus_fastq

    script:
    """
        medaka smolecule \
          --threads $cpus \
          --length 50 \
          --depth 2 \
          --model ${params.medaka_model} \
          --method spoa . \
          --qualities \
          ${smolecule_clusters_fastq} 2> smolecule.log
        
        mv consensus.fastq ${smolecule_clusters_fastq.baseName}_consensus.fastq
    """
}
