process FREEBAYES {
    publishDir "${params.output}/${params.variant_caller}/${sample}/", mode: 'copy'

    input:
      tuple val( sample ), val( target ), path( bam ), path( bai )
      val( type )
      path reference
      path reference_fai

    output:
      path "${type}.vcf", emit: variants

    script:
    """
      freebayes --version
      freebayes -f ${reference} ${bam} | vcfallelicprimitives -kg >${type}.vcf
    """
}
