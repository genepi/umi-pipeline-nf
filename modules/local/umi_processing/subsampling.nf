process SUBSAMPLING {
  publishDir "${params.output}/${sample}/subsampling", pattern: '*.fastq', mode: 'copy'
  publishDir "${params.output}/${sample}/stats", pattern: '*.tsv', mode: 'copy'
  
  input:
      tuple val( sample ), val( target ), path( merged_fastq )
  
  output:
      tuple val( "${sample}" ), val( "${target}" ), path( "${merged_fastq.baseName}_subsampled.fastq"), emit: subsampled_fastq
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