consensus_fasta="consensus.fasta"
vsearch_dir="vsearch_clusters"
detected_umis_file_name="dected_umis.fastq"
process CLUSTER {
    publishDir "${params.output}/${sample.baseName}/clustering/${type}", pattern: "cluster*", mode: 'copy'
        
    input:
        path sample
    output:
        tuple val( "${sample.baseName}" ), path( "cluster*" ), optional: true, emit:cluster_fastas
        
    script:
    """
        vsearch \
        --clusterout_id \
        --clusters cluster \
        --consout ${consensus_fasta} \
        --minseqlength ${params.min_length} \
        --maxseqlength ${params.max_length} \
        --threads ${params.threads} \
        --cluster_fast ${sample}/${detected_umis_file_name} \
        --clusterout_sort \
        --gapopen 0E/5I \
        --gapext 0E/2I \
        --mismatch -8 \
        --match 6 \
        --iddef 0 \
        --minwordmatches 0 \
        --qmask none \
        --id ${params.vsearch_sequence_identity}
    """
}
