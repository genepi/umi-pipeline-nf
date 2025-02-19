process CONTINUE_PIPELINE {
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