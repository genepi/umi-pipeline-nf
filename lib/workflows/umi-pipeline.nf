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
umi_split_cluster_python    = file( "${projectDir}/bin/split_cluster.py", checkIfExists: true )

// subdirectory and file prefixes
raw                         = "raw"
consensus                   = "consensus"
final_consensus             = "final"
barcode_sizes               = [:]


// STAGE CHANNELS
// to remove also barcode01 use :~/.*barcode((0[2-9])|([1-9][0-9]))/ 
// Remove barcode01 and uncalssified from the input fastq folder
Channel.fromPath("${params.input}/*", type: 'dir')
    .filter( ~/.*barcode(([0-9][0-9]))/ )
    .set { fastq_files_ch }

////////////////////
// BEGIN PIPELINE //
////////////////////

include {COPY_BED} from '../processes/copy_bed.nf'
include {MERGE_FASTQ} from '../processes/merge_input.nf'
include {MERGE_CONSENSUS_FASTQ} from '../processes/merge_consensus_fastq.nf'
include {SUBSAMPLING} from '../processes/subsampling.nf'
include {MAP_READS; MAP_READS as MAP_CONSENSUS; MAP_READS as MAP_FINAL_CONSENSUS} from '../processes/map_reads.nf'
include {SPLIT_READS} from  '../processes/split_reads.nf'
include {DETECT_UMI_FASTQ} from '../processes/detect_umi_fastq.nf'
include {CLUSTER; CLUSTER as CLUSTER_CONSENSUS} from '../processes/cluster.nf'
include {SPLIT_CLUSTER; SPLIT_CLUSTER as SPLIT_CLUSTER_CONSENSUS} from '../processes/cluster.nf'
include {REFORMAT_FILTER_CLUSTER} from '../processes/reformat_filter_cluster.nf'
include {POLISH_CLUSTER} from '../processes/polish_cluster.nf'
include {DETECT_UMI_CONSENSUS_FASTQ} from '../processes/detect_umi_consensus_fastq.nf'
include {REFORMAT_CONSENSUS_CLUSTER} from '../processes/reformat_consensus_cluster.nf'
include {LOFREQ} from '../processes/variant_calling/lofreq.nf'
include {MUTSERVE} from '../processes/variant_calling/mutserve.nf'
include {FREEBAYES} from '../processes/variant_calling/freebayes.nf'

// SUB-WORKFLOWS
workflow UMI_PIPELINE {

    main:
        COPY_BED( bed )

        if( params.subsampling ){
            MERGE_FASTQ( fastq_files_ch )
            SUBSAMPLING( MERGE_FASTQ.out.merged_fastq )
            SUBSAMPLING.out.subsampled_fastq
            .set { merged_fastq }
        } else {
            MERGE_FASTQ( fastq_files_ch )
            MERGE_FASTQ.out.merged_fastq
            .set { merged_fastq }
        }

        merged_fastq
        .filter { sample, target, fastq_file -> fastq_file.countFastq() > params.min_reads_per_barcode }
        .set { merged_filtered_fastq }

        MAP_READS( merged_filtered_fastq, raw, reference )
        SPLIT_READS( MAP_READS.out.bam_consensus, COPY_BED.out.bed, raw, umi_filter_reads )
        DETECT_UMI_FASTQ( SPLIT_READS.out.split_reads_fastx, raw, umi_extract )
        CLUSTER( DETECT_UMI_FASTQ.out.umi_extract_fastq, raw )

        // filter clusters that are below the minimal reads per cluster threshold
        // usually up to 80% of the clusters are below the threshold
        CLUSTER.out.cluster_fastqs
        .filter { sample, target, cluster_fastq -> cluster_fastq.countFastq() > params.min_reads_per_cluster}
        .set { cluster_fastqs }

        SPLIT_CLUSTER( cluster_fastqs, raw, umi_split_cluster_python )

        SPLIT_CLUSTER.out.split_cluster_fastqs
        .filter { sample, target, split_cluster_fastq -> split_cluster_fastq.countFastq() > params.min_reads_per_cluster}
        .set { split_cluster_fastqs }

        REFORMAT_FILTER_CLUSTER( split_cluster_fastqs, raw, umi_parse_clusters )

        // count number of final clusters per sample to end cluster polishing sooner
        // TODO test for small number of clusters if pipeline gets stuck
        REFORMAT_FILTER_CLUSTER.out.smolecule_cluster_fastq
        .groupTuple(by: 1)
        .map{ sample, type, fastqs -> barcode_sizes.put("$sample", fastqs.size)}

        POLISH_CLUSTER( flatten_smolecule_fastqs, consensus )

        POLISH_CLUSTER.out.consensus_fastq
        .map{ sample, type, fastq -> tuple( groupKey(sample, barcode_sizes.get("$sample")), type, fastq) }
        .groupTuple()
        .set { merge_consensus }

        MERGE_CONSENSUS_FASTQ(merge_consensus)
        
        MAP_CONSENSUS( MERGE_CONSENSUS_FASTQ.out.merged_consensus_fastq, consensus, reference )
        DETECT_UMI_CONSENSUS_FASTQ( MERGE_CONSENSUS_FASTQ.out.merged_consensus_fastq, consensus, umi_extract )
        CLUSTER_CONSENSUS( DETECT_UMI_CONSENSUS_FASTQ.out.umi_extract_fastq , consensus )
        REFORMAT_CONSENSUS_CLUSTER( CLUSTER_CONSENSUS.out.consensus_fastq, final_consensus, umi_reformat_consensus )
        MAP_FINAL_CONSENSUS( REFORMAT_CONSENSUS_CLUSTER.out.consensus_fastq, final_consensus, reference )
        
        if( params.call_variants ){
            if( params.variant_caller == "lofreq" ){
                LOFREQ( MAP_FINAL_CONSENSUS.out.bam_consensus, final_consensus, reference, reference_fai )
            }else if( params.variant_caller == "mutserve"){
                MUTSERVE( MAP_FINAL_CONSENSUS.out.bam_consensus, final_consensus, COPY_BED.out.bed, reference, reference_fai )
            }else if( params.variant_caller == "freebayes"){
                FREEBAYES( MAP_FINAL_CONSENSUS.out.bam_consensus, final_consensus, reference, reference_fai )
            }else{
                exit 1, "${params.variant_caller} is not a valid option. \nPossible variant caller are <lofreq/mutserve/freebayes>"
            
            }
        }
        
        
        
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