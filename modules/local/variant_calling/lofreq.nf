process LOFREQ {
  publishDir "${params.output}/${params.variant_caller}/${sample}/${target}/${type}", mode: 'copy'

  input:
  tuple val(sample), val(target), path(bam), path(bai)
  val type
  path reference
  path reference_fai

  output:
  tuple val("${sample}"), path("${type}.vcf"), emit: variants

  script:
  """
      lofreq call \
      --ref ${reference} \
      --out ${type}.vcf \
      --call-indels \
      --min-cov 5 \
      --no-default-filter \
      ${bam}
    """
}
