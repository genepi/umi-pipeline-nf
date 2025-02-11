process GLUE_CLUSTERS {
    tag "${sample}"
    publishDir "${params.output}/${sample}/clustering/combined", pattern: "${sample}_combined_clusters_*.${params.output_format}", mode: 'copy'

    input:
    tuple val(sample), val (type), path(cluster_files)

    output:
    // Emit all files created in the "combined_files" folder.
    tuple val(sample), val (type), path("*${params.output_format}"), emit: glued_clusters

    script:
    """
    #!/bin/bash
    set -euo pipefail

    clusters_per_file=${params.clusters_per_sample}

    # Initialize file and cluster counters.
    file_index=1
    cluster_count=0
    output_file="${sample}_combined_clusters_\${file_index}.${params.output_format}"
    
    # Loop over each provided cluster file.
    for file in ${cluster_files.join(' ')}; do
      # If the current output file already holds the allowed number of clusters,
      # then start a new file.
      if [ "\$cluster_count" -ge "\$clusters_per_file" ]; then
         file_index=\$((file_index+1))
         cluster_count=0
         output_file="${sample}_combined_clusters_\${file_index}.${params.output_format}"
      fi

      # Increase the count for the current file.
      cluster_count=\$((cluster_count+1))
      
      # Process the current cluster file, appending the cluster number to each header.
      if [ "${params.output_format}" = "fasta" ]; then
        # For FASTA, header lines start with ">"
        awk -v cnum=\$cluster_count '{ if(/^>/) { print \$0"_"cnum } else { print \$0 } }' "\$file" >> "\$output_file"
      elif [ "${params.output_format}" = "fastq" ]; then
        # For FASTQ, only modify the header line in each 4-line record.
        awk -v cnum=\$cluster_count '{
          if ((NR-1) % 4 == 0) {
            print \$0"_"cnum
          } else {
            print \$0
          }
        }' "\$file" >> "\$output_file"

      fi
    done
    """
}
