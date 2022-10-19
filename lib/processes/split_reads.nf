stats_filename = "umi_filter_reads_stats.txt"
process SPLIT_READS {

    publishDir "${params.output}/${sample}/stats/${type}", mode: 'copy', pattern: "${stats_filename}"
    publishDir "${params.output}/${sample}/fasta_filtered/${type}", mode: 'copy', pattern: "*.fastq"

    input:
        tuple val( sample ), val( target ), path ( bam ) , path ( bam_bai )
        path bed
        val( type )
        path python_filter_reads
    output:
        path "${stats_filename}"
        tuple val ( "${sample}" ), val( "target" ), path ( "*.fastq" ), emit: split_reads_fastq

    """
        python ${python_filter_reads} --min_overlap ${params.min_overlap} -o . ${bed} ${bam} 2>&1 | \
        tee ${stats_filename}
    """
}