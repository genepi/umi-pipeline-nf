nextflow.enable.dsl = 2

requiredParams = [
    'input', 'reference', 'bed'
]

for (param in requiredParams) {
    if (params[param] == null) {
      exit 1, "Parameter ${param} is required."
    }
}

// DEFINE PATHS
bed = file("${params.bed}", checkIfExists: true)
reference = file("${params.reference}", checkIfExists: true)

// python scripts
umi_filter_reads = file( "${projectDir}/bin/filter_reads.py", checkIfExists: true)
umi_extract = file( "${projectDir}/bin/extract_umis.py", checkIfExists: true)
umi_parse_clusters = file( "${projectDir}/bin/parse_clusters.py", checkIfExists: true)

// STAGE CHANNELS
fastq_files_ch = Channel.fromPath("${params.input}/*", type: 'dir')

// file_prefixes
raw = ""
consensus = "consensus"
final_consensus = "final"


////////////////////
// BEGIN PIPELINE //
////////////////////

// INCLUDES # here you must give the relevant process files from the lib directory 
include {COPY_BED} from '../processes/copy_bed.nf'
include {MAP_1D} from '../processes/map_1d.nf'
include {SPLIT_READS} from  '../processes/split_reads.nf'
include {DETECT_UMI_FASTA; DETECT_UMI_FASTA as DETECT_UMI_CONSENSUS_FASTA} from '../processes/detect_umi_fasta.nf'
include {CLUSTER, CLUSTER as CLUSTER_CONSENSUS} from '../processes/cluster.nf'
include {REFORMAT_FILTER_CLUSTER, REFORMAT_FILTER_CLUSTER as REFORMAT_FILTER_CLUSTER_CONSENSUS} from '../processes/reformat_filter_cluster.nf'
include {POLISH_CLUSTER} from '../processes/polish_cluster.nf'
include {MAP_CONSENSUS; MAP_CONSENSUS as MAP_FINAL_CONSENSUS} from '../processes/map_consensus.nf'
//include {DETECT_UMI_CONSENSUS_FASTA} from '../processes/detect_umi_consensus_fasta.nf'

// SUB-WORKFLOWS
workflow UMI_PIPELINE {

    main:

        COPY_BED( bed )
        MAP_1D( fastq_files_ch, reference )
        SPLIT_READS( MAP_1D.out.bam_1d, COPY_BED.out.bed, umi_filter_reads )
        DETECT_UMI_FASTA( SPLIT_READS.out.split_reads_fastq, raw, umi_extract )
        CLUSTER( DETECT_UMI_FASTA.out.umi_extract_fasta, raw )
        REFORMAT_FILTER_CLUSTER( CLUSTER.out.consensus_fasta, CLUSTER.out.vsearch_dir, umi_parse_clusters)
        POLISH_CLUSTER( REFORMAT_FILTER_CLUSTER.out.smolecule_clusters_fasta )
        MAP_CONSENSUS( POLISH_CLUSTER.out.consensus_fasta, reference )
        CLUSTER_CONSENSUS( DETECT_UMI_CONSENSUS_FASTA.out.umi_extract_fasta )

        DETECT_UMI_CONSENSUS_FASTA( POLISH_CLUSTER.out.consensus_fasta, final_consensus, umi_extract )
        MAP_FINAL_CONSENSUS( DETECT_UMI_CONSENSUS_FASTA.out.umi_extract_fasta, reference )

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
