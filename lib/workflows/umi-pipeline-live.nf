    nextflow.enable.dsl = 2

    requiredParams = [
        'input', 'reference', 'reference_fai', 'bed', 'output'
    ]

    for (param in requiredParams) {
        if (params[param] == null) {
        exit 1, "Parameter ${param} is required."
        }
    }

    // file paths
    bed                         = file("${params.bed}", checkIfExists: true)
    reference                   = file("${params.reference}", checkIfExists: true)
    reference_fai               = file( "${params.reference_fai}", checkIfExists: true)

    // python scripts
    umi_filter_reads            = file( "${projectDir}/bin/filter_reads.py", checkIfExists: true)
    umi_extract                 = file( "${projectDir}/bin/extract_umis.py", checkIfExists: true)
    umi_parse_clusters          = file( "${projectDir}/bin/parse_clusters.py", checkIfExists: true)
    umi_reformat_consensus      = file( "${projectDir}/bin/reformat_consensus.py", checkIfExists: true )
    umi_cluster_report          = file( "${projectDir}/bin/cluster_report.py", checkIfExists: true )

    // subdirectory and file prefixes
    raw                         = "raw"
    consensus                   = "consensus"
    final_consensus             = "final"
    n_parsed_cluster            = [:]


    ////////////////////
    // BEGIN PIPELINE //
    ////////////////////

    include {COPY_BED} from '../processes/copy_bed.nf'
    include {CONTINUE_PIPELINE} from '../processes/continue_pipeline.nf'
    include {MERGE_FASTQ} from '../processes/merge_input.nf'
    include {MERGE_CONSENSUS_FASTQ} from '../processes/merge_consensus_fastq.nf'
    include {MAP_READS; MAP_READS as MAP_CONSENSUS; MAP_READS as MAP_FINAL_CONSENSUS} from '../processes/map_reads.nf'
    include {SPLIT_READS} from  '../processes/split_reads.nf'
    include {DETECT_UMI_FASTQ; DETECT_UMI_FASTQ as DETECT_UMI_CONSENSUS_FASTQ} from '../processes/detect_umi_fastq.nf'
    include {CLUSTER as CLUSTER_CONSENSUS} from '../processes/cluster.nf'
    include {CLUSTER_LIVE} from '../processes/cluster_live.nf'
    include {REFORMAT_FILTER_CLUSTER} from '../processes/reformat_filter_cluster.nf'
    include {CLUSTER_STATS_LIVE} from '../processes/cluster_stats_live.nf'
    include {GLUE_CLUSTERS} from '../processes/glue_clusters.nf'
    include {POLISH_CLUSTER} from '../processes/polish_cluster.nf'
    include {FILTER_CONSENSUS_FASTQ} from '../processes/filter_consensus_fastq.nf'
    include {REFORMAT_CONSENSUS_CLUSTER} from '../processes/reformat_consensus_cluster.nf'
    include {LOFREQ as LOFREQ_CONSENSUS; LOFREQ as LOFREQ_FINAL_CONSENSUS} from '../processes/lofreq.nf'
    include {MUTSERVE as MUTSERVE_CONSENSUS; MUTSERVE as MUTSERVE_FINAL_CONSENSUS} from '../processes/mutserve.nf'
    include {FREEBAYES as FREEBAYES_CONSENSUS; FREEBAYES as FREEBAYES_FINAL_CONSENSUS} from '../processes/freebayes.nf'


    // SUB-WORKFLOWS
    workflow UMI_PIPELINE_LIVE {

        main:

            println("We are live!")

            COPY_BED( bed )

            Channel
                .watchPath("${params.output}/CONTINUE")
                .take(1)
                .set{ continue_ch }

            CONTINUE_PIPELINE( continue_ch )
            

            Channel
            .fromPath("${params.input}/barcode*/*.fastq")
            .set{ existing_files_ch }

            Channel
            .watchPath("${params.input}/barcode*/*.fastq", 'create, modify')
            .until { it.getFileName().toString().toLowerCase().contains("continue") } 
            .set { watched_files_ch }

            existing_files_ch
            .concat( watched_files_ch )
            .map{ 
                fastq -> 
                def barcode = fastq.parent.name
                tuple(barcode, fastq)
                }
            .splitFastq( by: params.chunk_size , file: true)
            .set{ ch_input_files }
            
            MERGE_FASTQ( ch_input_files )
            .set { merged_fastq }

            MAP_READS( merged_fastq, raw, reference )
            SPLIT_READS( MAP_READS.out.bam_consensus, COPY_BED.out.bed, raw, umi_filter_reads )

            SPLIT_READS.out.split_reads_fastx
            .filter{ _sample, _target, fastq -> fastq.countFastq() > 0 }
            .set{ splt_reads_filtered }

            DETECT_UMI_FASTQ( splt_reads_filtered, raw, umi_extract )

            Channel
                .watchPath("${params.output}/*/${params.output_format}_umi/raw/*.fastq", 'create, modify')
                .until { it.getFileName().toString().toLowerCase().contains("continue") }
                .map { path -> 
                    def barcode = path.parent.parent.parent.name
                    def files = path.parent.listFiles().findAll { it.name.endsWith('.fastq') }
                    tuple(barcode, "target", files)
                }
                .set { cluster_ch }
            
            CLUSTER_LIVE( cluster_ch, raw )

            CLUSTER_LIVE.out.cluster_fastas
                .map { barcode, target, clusters -> 
                    def filtered_clusters = clusters.findAll { fasta -> fasta.countFasta() > params.min_reads_per_cluster }
                    filtered_clusters ? [barcode, target, filtered_clusters] : null
                }
                .filter { it != null }
                .set{ cluster_fastas }
            
            REFORMAT_FILTER_CLUSTER( cluster_fastas, raw, umi_parse_clusters )

            // Assuming that REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs produces tuples of (sample, type, fastqFiles)
            REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs
                .filter { _sample, _type, fastqs -> fastqs instanceof List && fastqs.size() > 0 }  // only pass non-empty lists
                .set { smolecule_cluster_fastqs_list }

            // Launch the reporting process for each sample
            CLUSTER_STATS_LIVE( REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_stats, umi_cluster_report )
            
            //SUMMARY_CLUSTER_STATS( CLUSTER_STATS_LIVE.out.cluster_stats)

/*

            REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastqs
            .combine( continue_ch )
            .map { sample, type, fastqs, _continue_file -> tuple(sample, type, fastqs) }
            .filter{ _sample, _type, fastqs -> fastqs.class == ArrayList}
            .set{ smolecule_cluster_fastqs_list }     

            Channel
                .watchPath("${params.output}/CONTINUE")
                .flatMap { 
                        // The CONTINUE file is present now – list all matching smolecule files
                        Channel.fromPath("${params.output}//clustering/raw/smolecule*")
                }
                .view()
                .map{ 
                    smolecule -> 
                    def barcode = smolecule.parent.parent.parent.name
                    tuple(barcode, smolecule)
                }
                .groupTuple()
                .set{ smolecule_cluster_fastqs_list }

            GLUE_CLUSTERS(smolecule_cluster_fastqs_list)
                .map { sample, target, clusters -> tuple(sample, target, clusters instanceof List ? clusters : [clusters]) }
                .transpose(by: 2)
                .set { glued_clusters }

            POLISH_CLUSTER( glued_clusters, consensus )
            
            POLISH_CLUSTER.out.consensus_fastq
            //.map{ sample, type, fastq -> tuple( groupKey(sample, n_parsed_cluster.get("$sample")), type, fastq) }
            .groupTuple( )
            .set{ merge_consensus }

            
            if ( params.output_format == "fastq"){
                MERGE_CONSENSUS_FASTQ(merge_consensus, consensus)
                FILTER_CONSENSUS_FASTQ(MERGE_CONSENSUS_FASTQ.out.merged_consensus_fastq, consensus)
                FILTER_CONSENSUS_FASTQ.out.filtered_consensus_fastq
                .set{ consensus_fastq }
            } else {
                MERGE_CONSENSUS_FASTQ(merge_consensus, consensus)
                .set{ consensus_fastq }
            }

            MAP_CONSENSUS( consensus_fastq, consensus, reference )
            DETECT_UMI_CONSENSUS_FASTQ( consensus_fastq, consensus, umi_extract )
            CLUSTER_CONSENSUS( DETECT_UMI_CONSENSUS_FASTQ.out.umi_extract_fastq , consensus )
            REFORMAT_CONSENSUS_CLUSTER( CLUSTER_CONSENSUS.out.consensus_fasta, final_consensus, umi_reformat_consensus )
            MAP_FINAL_CONSENSUS( REFORMAT_CONSENSUS_CLUSTER.out.consensus_fastq, final_consensus, reference )
            
            if( params.call_variants ){
                if( params.variant_caller == "lofreq" ){
                    LOFREQ_CONSENSUS( MAP_CONSENSUS.out.bam_consensus, consensus, reference, reference_fai )
                    LOFREQ_FINAL_CONSENSUS( MAP_FINAL_CONSENSUS.out.bam_consensus, final_consensus, reference, reference_fai )
                }else if( params.variant_caller == "mutserve"){
                    MUTSERVE_CONSENSUS( MAP_CONSENSUS.out.bam_consensus, consensus, COPY_BED.out.bed, reference, reference_fai )
                    MUTSERVE_FINAL_CONSENSUS( MAP_FINAL_CONSENSUS.out.bam_consensus, final_consensus, COPY_BED.out.bed, reference, reference_fai )
                }else if( params.variant_caller == "freebayes"){
                    FREEBAYES_CONSENSUS( MAP_CONSENSUS.out.bam_consensus, consensus, reference, reference_fai )
                    FREEBAYES_FINAL_CONSENSUS( MAP_FINAL_CONSENSUS.out.bam_consensus, final_consensus, reference, reference_fai )
                }else{
                    exit 1, "${params.variant_caller} is not a valid option. \nPossible variant caller are <lofreq/mutserve/freebayes>"
                
                }
            }  
            
*/
    }


    //////////////////
    // END PIPELINE //
    //////////////////

    // WORKFLOW TRACING # what to display when the pipeline finishes
    // eg. with errors
    workflow.onError {
        log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
    }

    // eg. in general
    workflow.onComplete {

        log.info ""
        log.info "         Pipeline execution summary"
        log.info "         ---------------------------"
        log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
        log.info "         Profile      : ${workflow.profile}"
        log.info "         Launch dir   : ${workflow.launchDir}"    
        log.info "         Work dir     : ${workflow.workDir} ${!params.debug && workflow.success ? "(cleared)" : "" }"
        log.info "         Status       : ${workflow.success ? "success" : "failed"}"
        log.info "         Error report : ${workflow.errorReport ?: "-"}"
        log.info ""

        
        // run a small clean-up script to remove "work" directory after successful completion 
        if (!params.debug && workflow.success) {
            ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
        
    }