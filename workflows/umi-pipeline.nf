nextflow.enable.dsl = 2

include { LIVE_UMI_PROCESSING       } from './live-umi-processing.nf'
include { OFFLINE_UMI_PROCESSING    } from './offline-umi-processing.nf'
include { UMI_POLISHING             } from './umi-polishing.nf'
include { VARIANT_CALLING           } from './variant_calling.nf'
include { COPY_BED                  } from '../modules/local/copy_bed.nf'

workflow UMI_PIPELINE {

    main:

        //validate input parameters
        WorkflowMain.validate(params)
        
        // file paths
        bed                         = file("${params.bed}", checkIfExists: true)
        reference                   = file("${params.reference}", checkIfExists: true)
        reference_fai               = file( "${params.reference_fai}", checkIfExists: true)

        // python scripts
        umi_filter_reads            = file( "${projectDir}/bin/filter_reads.py", checkIfExists: true)
        umi_extract                 = file( "${projectDir}/bin/extract_umis.py", checkIfExists: true)
        umi_parse_clusters          = file( "${projectDir}/bin/parse_clusters_old.py", checkIfExists: true)
        umi_reformat_consensus      = file( "${projectDir}/bin/reformat_consensus.py", checkIfExists: true )
        umi_cluster_report          = file( "${projectDir}/bin/cluster_report.py", checkIfExists: true )
        umi_cluster_stats_summary   = file( "${projectDir}/bin/summary_cluster_report.py", checkIfExists: true )

        // subdirectory and file prefixes
        raw                         = "raw"
        consensus                   = "consensus"
        final_consensus             = "final"
        n_parsed_cluster            = [:]
        // cluster_summary_output_path = "${params.output}/cluster_stats/summary_cluster_stats.tsv"
        def extracted_fastq_cache_dir = new File (".nextflow/cache/${workflow.sessionId}/extracted_fastq_cache_dir")
        extracted_fastq_cache_dir.mkdir()
        extracted_fastq_cache_dir_nf = file( extracted_fastq_cache_dir )


        COPY_BED( bed )

        COPY_BED.out.bed
            .set{ bed_ch }

        if ( params.live ){        
            LIVE_UMI_PROCESSING(
                raw,
                reference,
                umi_filter_reads,
                extracted_fastq_cache_dir_nf,
                umi_extract,
                umi_parse_clusters,
                umi_cluster_report,
                umi_cluster_stats_summary,
                bed_ch
                )
            
            LIVE_UMI_PROCESSING.out.processed_umis
                .set{processed_umis}

        } else {            
            OFFLINE_UMI_PROCESSING(
                raw,
                reference,
                umi_filter_reads,
                umi_extract,
                umi_parse_clusters,
                bed_ch
            )

            OFFLINE_UMI_PROCESSING.out.processed_umis
                .set{ processed_umis }
        }

        UMI_POLISHING(
            processed_umis,
            n_parsed_cluster,
            consensus,
            final_consensus,
            reference,
            umi_extract,
            umi_reformat_consensus
        )

        VARIANT_CALLING(
            UMI_POLISHING.out.consensus_bam,
            UMI_POLISHING.out.final_consensus_bam,
            consensus,
            final_consensus,
            reference,
            reference_fai,
            bed_ch
        )
}