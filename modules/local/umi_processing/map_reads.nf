process MAP_READS {
    publishDir "${params.output}/${sample}/${type}/align/", mode: 'copy'

    input:
        tuple val( sample ), path( consensus_fastq )
        val( type )
        path reference
    output:
        tuple val( "${sample}"), path ( "*.bam" ), path ( "*.bam.bai" ), emit: bam_consensus

    script:
    """
        minimap2 \
          ${params.minimap2_param} \
          -t ${params.threads} \
          ${reference} \
          ${consensus_fastq} | 
        samtools sort \
          -@ ${params.threads} \
          -o ${consensus_fastq.baseName}.bam - && \
        samtools index \
          -@ ${params.threads} \
          ${consensus_fastq.baseName}.bam
    """
}