process MUTSERVE {
    publishDir "${params.output}/${params.variant_caller}/${sample}/", mode: 'copy'
  input:
    tuple val( sample ), val( target ), path( bam ), path( bai )
    val( type )
    path reference
  output:
    path "${type}.vcf", emit: variants
  script:
  """
    mutserve call --output ${type}.vcf --write-raw --reference ${reference} --deletions --contig-name ${sample}
  """
}