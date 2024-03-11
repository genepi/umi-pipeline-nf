process MUTSERVE {
    publishDir "${params.output}/${params.variant_caller}/${sample}/${type}", mode: 'copy'
  
    input:
      tuple val( sample ), val( target ), path( bam ), path( bai )
      val( type )
      path bed
      path reference
      path reference_fai
    
    output:
      path "${type}.txt", emit: variants
      path "${type}_raw.txt", emit: variants_raw
    
    script:
    """
      java -jar /opt/mutserve_LPA_adapted.jar call  \
      --output ${type}.vcf \
      --write-raw \
      --reference ${reference} \
      --insertions \
      --deletions \
      --contig-name \$(awk '{ print \$1 }' ${bed}) \
      ${bam}
    """
}



