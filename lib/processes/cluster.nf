consensus_fasta="consensus.fasta"
detected_umis_file_name="detected_umis.fastq"

process CLUSTER {
    // publishDir "${params.output}/${sample}/clustering/${type}", pattern: "${consensus_fasta}", mode: 'copy'
    
    input:
        path sample
        val ( type )
    output:
        tuple val( "${sample.baseName}" ), val( "lpa5104" ), path( "${consensus_fasta}" ), emit:consensus_fasta
        tuple val( "${sample.baseName}" ), val( "lpa5104" ), path( "cluster*" ), emit:cluster_fastas
        
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
