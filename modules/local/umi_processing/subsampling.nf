process SUBSAMPLING {
  publishDir "${params.output}/${sample}/${type}/subsampling", pattern: '*.fastq', mode: 'copy'
  publishDir "${params.output}/${sample}/${type}/subsampling/stats", pattern: '*.tsv', mode: 'copy'
  
  input:
    tuple val( sample ), path( merged_fastq )
    val( type )
  
  output:
    tuple val( "${sample}" ), path( "${merged_fastq.baseName}_subsampled.fastq"), emit: subsampled_fastq
    path( "${merged_fastq.baseName}_subsampled.tsv" )
  
  script:
  """
    echo -e "seed\tsubsampling_readnumber\n${params.subsampling_seed}\t${params.subsampling_readnumber}" > ${merged_fastq.baseName}_subsampled.tsv
    seqtk sample \
      -s ${params.subsampling_seed} \
      ${merged_fastq} \
      ${params.subsampling_readnumber} > ${merged_fastq.baseName}_subsampled.fastq 
  """
}