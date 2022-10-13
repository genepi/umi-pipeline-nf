process MAP_1D {
    publishDir "${params.output}/${fastq.baseName}/align", mode: 'copy'

    input:
        path fastq
        path reference 
    output:
        tuple val( "${fastq.baseName}" ), path ( "1d.bam" ) , emit: bam_1d
        path "1d.bam.bai", emit: bai_1d

    """
        catfishq --max_n 0 ${fastq} | minimap2 ${params.minimap2_param} -t ${params.threads} ${reference} - | samtools sort -@ 5 -o 1d.bam - && samtools index -@ ${params.threads} 1d.bam
    """
}

// Note: --max_n was a config paramter in the snakemakepipeline 
// which was not used in the end. default was 0, which is used now
// Step could be used for subsample but make extra process for that