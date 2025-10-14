process MAP_READS {
    publishDir "${params.output}/${sample}/${type}/align/", mode: 'copy', enabled: "${params.verbose}"

    input:
        tuple val( sample ), path( fastq )
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
          ${fastq} |
        samtools sort \
          -@ ${params.threads} \
          -o ${fastq.baseName}.bam - && \
        samtools index \
          -@ ${params.threads} \
          ${fastq.baseName}.bam
    """
}