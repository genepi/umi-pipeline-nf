process  {
  input:
    tuple val( sample ), val( target ), path( merged_fastq )
  output:
    tuple tuple val( "${sample.baseName}" ), val( "target" ), path( "${subsampled_fastq}"), emit: subsampled_fastq
  script:
  """
    
  """
}