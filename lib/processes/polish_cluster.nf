process POLISH_CLUSTER {
    tag "${sample}"
    conda "/root/miniconda3/envs/medaka"

    input:
        tuple val( sample ), val( target ), path( smolecule_clusters_fasta )
        val ( type )

    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${smolecule_clusters_fasta.baseName}_consensus.fasta" ), emit: consensus_fasta

    script:
    """
        medaka smolecule \
          --threads ${params.threads} \
          --length 50 \
          --depth 2 \
          --model ${params.medaka_model} \
          --method spoa . \
          ${smolecule_clusters_fasta} 2> smolecule.log
        
        mv consensus.fasta ${smolecule_clusters_fasta.baseName}_consensus.fasta
    """
}
