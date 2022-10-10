process MAP_1D {
    publishDir "${params.output}", mode: 'copy'

    input:
        path fastq
        path reference 
    output:
        path "align/1d.bam", emit: bam_1d
        path "align/1d.bam.bai", emit: bai_1d

    """
        mkdir -p align
        catfishq --max_n 0 ${fastq} | \\
        minimap2 ${params.minimap2_param} -t ${params.threads} ${reference} - | \\
        samtools sort -@ 5 -o align/1d.bam - && samtools index -@ ${params.threads} align/1d.bam
    """
}

// Note: --max_n was a config paramter in the snakemakepipeline 
// which was not used in the end. default was 0, which is used now
// Step could be used for subsample but make extra process for that