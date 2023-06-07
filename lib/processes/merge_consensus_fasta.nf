merged_fasta="merged_consensus.fasta"
process MERGE_CONSENSUS_FASTA {
    input:
        tuple val( sample ), val( target ), path( consensus_fastas )
    output:
        tuple val( "${sample}" ), val( "${target}" ), path( "${merged_fasta}"), emit: merged_consensus_fasta
    """
        cat ${consensus_fastas} > ${merged_fasta}
    """
}