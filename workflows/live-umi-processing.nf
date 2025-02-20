nextflow.enable.dsl = 2

include {CONTINUE_PIPELINE} from '../modules/local/umi_processing/continue_pipeline.nf'
include {MERGE_FASTQ} from '../modules/local/umi_processing/merge_input.nf'
include {MAP_READS} from '../modules/local/map_reads.nf'
include {SPLIT_READS} from  '../modules/local/umi_processing/split_reads.nf'
include {DETECT_UMI_FASTQ} from '../modules/local/umi_processing/detect_umi_fastq.nf'
include {CLUSTER} from '../modules/local/umi_processing/cluster.nf'
include {REFORMAT_FILTER_CLUSTER} from '../modules/local/umi_processing/reformat_filter_cluster.nf'
include {CLUSTER_STATS_LIVE} from '../modules/local/umi_processing/cluster_stats_live.nf'
include {SUMMARY_CLUSTER_STATS} from '../modules/local/umi_processing/summary_cluster_stats.nf'

workflow LIVE_UMI_PROCESSING {
    take:
        raw
        reference
        umi_filter_reads
        extracted_fastq_cache_dir_nf
        umi_extract
        umi_parse_clusters
        umi_cluster_report
        umi_cluster_stats_summary
        cluster_summary_cache_dir_nf
        bed


    main:
       println("We are live!")

        Channel
            .watchPath("${params.output}/CONTINUE")
            .take(1)
            .set{ continue_ch }

        CONTINUE_PIPELINE( continue_ch )

        Channel
            .fromPath("${params.input}/barcode*/*.fastq")
            .set{ existing_fastqs }
        
        Channel
            .watchPath("${params.input}/barcode*/*.fastq", 'create, modify')
            .until { it.getFileName().toString().toLowerCase().contains("continue") } 
            .set { watched_fastqs }
            
        existing_fastqs
            .concat( watched_fastqs )
            .map{ 
                fastq -> 
                def barcode = fastq.parent.name
                tuple(barcode, fastq)
                }
            .splitFastq( by: params.chunk_size , file: true)
            .set{ chunked_input_fastqs }

                
        MERGE_FASTQ( chunked_input_fastqs )
        MAP_READS( MERGE_FASTQ.out.merged_fastq, raw, reference )
        SPLIT_READS( MAP_READS.out.bam_consensus, bed, raw, umi_filter_reads )

        SPLIT_READS.out.split_reads_fastx
        .filter{ _sample, _target, fastq -> fastq.countFastq() > params.min_reads_per_barcode }
        .set{ split_reads_filtered }

        DETECT_UMI_FASTQ( split_reads_filtered, extracted_fastq_cache_dir_nf, raw, umi_extract )

        CLUSTER( DETECT_UMI_FASTQ.out.umi_extract_fastq, raw )

        // Filters the clusters to only keep cluser with more or equal than min_reads_per_cluster, but keeps the grouping per sample
        CLUSTER.out.cluster_fastas
            .map { barcode, target, clusters -> 
                def filtered_clusters = clusters.findAll { fasta -> fasta.countFasta() >= params.min_reads_per_cluster }
                filtered_clusters ? [barcode, target, filtered_clusters] : null
            }
            .filter { it != null }
            .set{ cluster_fastas }
        
        REFORMAT_FILTER_CLUSTER( cluster_fastas, raw, umi_parse_clusters )

        CLUSTER_STATS_LIVE( REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_stats, umi_cluster_report )
        SUMMARY_CLUSTER_STATS( REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_stats, cluster_summary_cache_dir_nf, umi_cluster_stats_summary)

        REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs
            .combine( continue_ch )
            .map { sample, type, fastqs, task_index, _continue_file ->
                tuple(sample, type, fastqs, [task_index])
            }
            .filter { _sample, _type, fastqs, _task_index -> fastqs instanceof List }
            .groupTuple(by: [0, 1], sort: { it[3] })
            .map { sample, type, fastqs, _task_index ->
                tuple(sample, type, fastqs[0])
            }
            .set { processed_umis }


    emit:
        processed_umis

}