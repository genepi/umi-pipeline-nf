include {MERGE_FASTQ} from '../modules/local/umi_processing/merge_input.nf'
include {SUBSAMPLING} from '../modules/local/umi_processing/subsampling.nf'
include {MAP_READS} from '../modules/local/umi_processing/map_reads.nf'
include {SPLIT_READS} from  '../modules/local/umi_processing/split_reads.nf'
include {DETECT_UMI_CONSENSUS_FASTQ as DETECT_UMI_FASTQ} from '../modules/local/umi_polishing/detect_umi_consensus_fastq.nf'
include {CLUSTER} from '../modules/local/umi_processing/cluster.nf'
include {REFORMAT_FILTER_CLUSTER} from '../modules/local/umi_processing/reformat_filter_cluster.nf'
include {CLUSTER_STATS} from '../modules/local/umi_processing/cluster_stats.nf'
include {SUMMARY_CLUSTER_STATS} from '../modules/local/umi_processing/summary_cluster_stats.nf'


workflow OFFLINE_UMI_PROCESSING {

    take:
        raw
        reference
        umi_filter_reads
        umi_extract
        umi_parse_clusters
        umi_cluster_report
        umi_cluster_stats_summary
        cluster_summary_cache_dir_nf

        bed

    main:        

        def fastq_channel = params.single_sample ?
            Channel.fromPath("${params.input}/*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: true) :
            Channel.fromPath("${params.input}/barcode*/*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: true)
        
        fastq_channel
            .map{ 
                fastqs -> 
                def barcode = fastqs.parent.name
                tuple(barcode, fastqs)
            }
            .set{ existing_fastqs }

        if( params.subsampling ){
            existing_fastqs
            .groupTuple( by: 0 ) 
            .set{ existing_fastqs_grouped }
            
            MERGE_FASTQ( existing_fastqs_grouped )
            SUBSAMPLING( MERGE_FASTQ.out.merged_fastq, raw )
            SUBSAMPLING.out.subsampled_fastq
            .set { input_fastqs }
        } else {
            MERGE_FASTQ( existing_fastqs )
            .set { input_fastqs }
        }


        input_fastqs
            .splitFastq( by: params.chunk_size , file: true)
            .set{ chunked_input_fastqs }

        MAP_READS( chunked_input_fastqs, raw, reference )
        
        MAP_READS.out.bam_consensus
        .combine(bed)
        .set{ bam_consensus_bed_sets }

        SPLIT_READS( bam_consensus_bed_sets, raw, umi_filter_reads )

        SPLIT_READS.out.split_reads_fastx
        .filter{ _sample, _target, fastq -> fastq.countFastq() > params.min_reads_per_barcode }
        .set{ split_reads_filtered }

        DETECT_UMI_FASTQ( split_reads_filtered, raw, umi_extract )
        
        DETECT_UMI_FASTQ.out.umi_extract_fastq
        .groupTuple( by: [0, 1])
        .set{ extracted_umis }

        CLUSTER( extracted_umis, raw )

        // Filters the clusters to only keep cluser with more or equal than min_reads_per_cluster, but keeps the grouping per sample
        CLUSTER.out.cluster_fastas
            .map { barcode, target, clusters -> 
                def filtered_clusters = clusters.findAll { fasta -> fasta.countFasta() >= params.min_reads_per_cluster }
                filtered_clusters ? [barcode, target, filtered_clusters] : null
            }
            .filter { it != null }
            .set{ cluster_fastas }

        REFORMAT_FILTER_CLUSTER( cluster_fastas, raw, umi_parse_clusters )

        CLUSTER_STATS( REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_stats, raw, umi_cluster_report )
        SUMMARY_CLUSTER_STATS( REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_stats, cluster_summary_cache_dir_nf, umi_cluster_stats_summary)

        REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs
        .filter{ _sample, _type, fastqs, _task_index -> fastqs instanceof List}            
        .map{ sample, type, fastqs, _task_index ->
                tuple(sample, type, fastqs)
        }
        .set{ processed_umis }

        emit:
            processed_umis

}
