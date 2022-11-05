process FREEBAYES {
    publishDir "${params.output}/${params.variant_caller}/${sample}/", mode: 'copy'
  input:
    tuple val( sample ), val( target ), path( bam ), path( bai )
    val( type )
    path reference
  output:
    path "${type}.vcf", emit: variants
  script:
  """
    freebayes -f ${reference} ${bam} | vcfallelicprimitives -kg > ${type}.vcf
  
  """
}

//freebayes -f ${reference} ${bam} | vcfallelicprimitives -kg > ${type}.vcf

/*freebayes -f ${reference} --haplotype-length 0 --min-alternate-count 1 \
      --min-alternate-fraction 0 --pooled-continuous --report-monomorphic ${bam} > ${type}.vcf
*/    