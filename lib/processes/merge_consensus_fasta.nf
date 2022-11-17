merged_fasta="merged_consensus.fasta"
process MERGE_CONSENSUS_FASTA {
    input:
        tuple val( sample ), val( target ), path( sample )
    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${merged_fasta}"), emit: merged_consensus_fasta
    """
        cat ${sample} > ${merged_fasta}
    """
}