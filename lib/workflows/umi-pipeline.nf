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

// subdirectory and file prefixes
raw                         = "raw"
consensus                   = "consensus"
final_consensus             = "final"
n_parsed_cluster            = [:]


// STAGE CHANNELS
// to remove also barcode01 use :~/.*barcode((0[2-9])|([1-9][0-9]))/ 
// Remove barcode01 and uncalssified from the input fastq folder
Channel.fromPath("${params.input}/*", type: 'dir')
    .filter( ~/.*barcode(([0-9][0-9]))/ )
    .set { detected_umis_ch }

////////////////////
// BEGIN PIPELINE //
////////////////////

include {COPY_BED} from '../processes/copy_bed.nf'
include {MERGE_FASTQ} from '../processes/merge_input.nf'
include {MERGE_CONSENSUS_FASTQ} from '../processes/merge_consensus_fastq.nf'
include {SUBSAMPLING} from '../processes/subsampling.nf'
include {MAP_READS; MAP_READS as MAP_CONSENSUS; MAP_READS as MAP_FINAL_CONSENSUS} from '../processes/map_reads.nf'
include {SPLIT_READS} from  '../processes/split_reads.nf'
include {DETECT_UMI_FASTQ; DETECT_UMI_FASTQ as DETECT_UMI_CONSENSUS_FASTQ} from '../processes/detect_umi_fastq.nf'
include {CLUSTER; CLUSTER as CLUSTER_CONSENSUS} from '../processes/cluster.nf'
include {REFORMAT_FILTER_CLUSTER} from '../processes/reformat_filter_cluster.nf'


// SUB-WORKFLOWS
workflow UMI_PIPELINE {

    main:
        CLUSTER( detected_umis_ch, raw )

        CLUSTER.out.cluster_fastas
        .transpose( by: 1 )
        .filter{ sample, fasta -> fasta.countFasta() >= params.min_reads_per_cluster }
        .groupTuple()
        .set{ filtered_clusters }

        REFORMAT_FILTER_CLUSTER( filtered_clusters, raw, umi_parse_clusters )
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