process MUTSERVE {
    publishDir "${params.output}/${params.variant_caller}/${sample}/${target}/${type}", mode: 'copy'
  
    input:
      tuple val( target ), path( bed ), val( sample ), path( bam ), path( bai )
      val( type )
      path reference
      path reference_fai
    
    output:
      path "${type}.txt", emit: variants
      path "${type}_raw.txt", emit: variants_raw
      path "${type}_parsed.txt", emit: variants_parsed
    
    script:
    def min_variant_level = 0.0085
    """
      java -jar /opt/mutserve_LPA_adapted.jar call  \
        --output ${type}.vcf \
        --write-raw \
        --reference ${reference} \
        --insertions \
        --deletions \
        --contig-name \$(awk '{ print \$1 }' ${bed}) \
      ${bam} && \
      echo -e "Sample\tPosition\tReference\tVariant\tVariant-level" > ${type}_parsed.txt && \
      tail -n +2 ${type}_raw.txt | \
      awk -F'\t' '(\$3 == \$4 && \$27 >= $min_variant_level) || \$3 != \$4' | \
      awk -F'\t' 'BEGIN {OFS="\t"} {print \$1, \$2, \$3, (\$3 == \$4 ? \$5 : \$4), (\$3 != \$4 ? \$25 : \$27)}' >> ${type}_parsed.txt

    """
}
