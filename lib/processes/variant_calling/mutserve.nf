process MUTSERVE {
    publishDir "${params.output}/${params.variant_caller}/${sample}/", mode: 'copy'
  input:
    tuple val( sample ), val( target ), path( bam ), path( bai )
    val( type )
    path bed
    path reference
  output:
    path "${type}.txt", emit: variants
    path "${type}_raw.txt", emit: variants_raw
  script:
  """
    mutserve call --output ${type}.vcf --write-raw --reference ${reference} --deletions --contig-name \$(awk '{ print \$1 }' ${bed}) ${bam}
  """
}



