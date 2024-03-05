process FREEBAYES {
    publishDir "${params.output}/${params.variant_caller}/${sample}/${type}", mode: 'copy'

    input:
      tuple val( sample ), val( target ), path( bam ), path( bai )
      val( type )
      path reference
      path reference_fai

    output:
      path "${type}.vcf", emit: variants

    script:
    """
      freebayes -f ${reference} -F 0.009 -p 80 ${bam} | vcfallelicprimitives -kg >${type}.vcf
    """
}
