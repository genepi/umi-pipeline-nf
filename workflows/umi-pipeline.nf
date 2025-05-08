nextflow.enable.dsl = 2

include { LIVE_UMI_PROCESSING       } from './live-umi-processing.nf'
include { OFFLINE_UMI_PROCESSING    } from './offline-umi-processing.nf'
include { UMI_POLISHING             } from './umi-polishing.nf'
include { VARIANT_CALLING           } from './variant_calling.nf'
include { PARSE_BED                 } from '../modules/local/parse_bed.nf'

workflow UMI_PIPELINE {

    main:

        //validate input parameters
        WorkflowMain.validate(params)
        
        // file paths
        bed                         = file("${params.bed}", checkIfExists: true)
        reference                   = file("${params.reference}", checkIfExists: true)
        reference_fai               = file("${params.reference_fai}", checkIfExists: true)

        // python scripts
        umi_filter_reads            = file( "${projectDir}/bin/filter_reads.py", checkIfExists: true)
        umi_extract                 = file( "${projectDir}/bin/extract_umis.py", checkIfExists: true)
        umi_parse_clusters          = file( "${projectDir}/bin/parse_clusters.py", checkIfExists: true)
        umi_reformat_consensus      = file( "${projectDir}/bin/reformat_consensus.py", checkIfExists: true )
        umi_cluster_report          = file( "${projectDir}/bin/cluster_report.py", checkIfExists: true )
        umi_cluster_stats_summary   = file( "${projectDir}/bin/summary_cluster_report.py", checkIfExists: true )
        umi_parse_bam               = file( "${projectDir}/bin/parse_cluster_alignment.py", checkIfExists: true)

        // subdirectory and file prefixes
        raw                         = "raw"
        consensus                   = "consensus"
        final_consensus             = "final"
        n_parsed_cluster            = [:]
        // cluster_summary_output_path = "${params.output}/cluster_stats/summary_cluster_stats.tsv"
        def extracted_fastq_cache_dir = new File (".nextflow/cache/${workflow.sessionId}/extracted_fastq_cache_dir")
        extracted_fastq_cache_dir.mkdir()
        extracted_fastq_cache_dir_nf = file( extracted_fastq_cache_dir )
        def cluster_summary_cache_dir = new File (".nextflow/cache/${workflow.sessionId}/cluster_summary_cache_dir")
        cluster_summary_cache_dir.mkdir()
        cluster_summary_cache_dir_nf = file( cluster_summary_cache_dir )


        // Use Nextflow to split input bed file
        Channel.fromPath( bed )
        .splitText()
        .map { line ->
            def fields = line.tokenize('\t')
            def target = fields[3].trim()  // Trim to remove trailing \n or spaces
            return tuple(target, line)
        }
        .set{ bed_channel }
        PARSE_BED( bed_channel )

        PARSE_BED.out.bed_channel
            .set{ bed_ch }

        bed_ch.view()

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
                cluster_summary_cache_dir_nf,
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
                umi_cluster_report,
                umi_cluster_stats_summary,
                cluster_summary_cache_dir_nf,
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
            umi_parse_bam,
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