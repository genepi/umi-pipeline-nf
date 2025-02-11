process CONTINUE_PIPELINE {
  publishDir "${params.output}/continue/${params.output_format}_umi/raw/", mode: 'copy'
  publishDir "${params.input}/barcode_continue/", mode: 'copy'

    input:
    val continue_condition

    output:
    path "continue.fastq"

    script:
    """
      touch continue.fastq
    """

}