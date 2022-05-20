#!/usr/bin/env nextflow

// ENABLE DSL2
nextflow.enable.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         ===============================================
          E C S E Q - t e m p l a t e   P I P E L I N E
         ===============================================
         ~ version ${workflow.manifest.version}

         Usage: 
              nextflow run ecseq/dnaseq [OPTIONS]...

         Options: GENERAL
              --input [path/to/input/dir]     [REQUIRED] Provide the directory containing fastq file(s) in "*{1,2}.fastq.gz" format

              --reference [path/to/ref.fa]    [REQUIRED] Provide the path to the reference genome in fasta format

              --output [STR]                  A string that can be given to name the output directory. [default: "."]


         Options: MODIFIERS
              --SE                            Indicate to the pipeline whether fastq files are SE reads in "*.fastq.gz" format. [default: off]

              --FastQC                        Generate FastQC report of trimmed reads. [default: off]

              --bamQC                         Generate bamQC report of alignments. [default: off]

              --keepReads                     Keep trimmed fastq reads. [default: off]


         Options: TRIMMING
              --forward                       Forward adapter sequence. [default: "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"]

              --reverse                       Reverse adapter sequence. [default: "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"]

              --minQual                       Minimum base quality threshold. [default: 20]

              --minLeng                       Minimum read length threshold. [default: 25]

              --minOver                       Minimum overlap threshold. [default: 3]


         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    


         Example: 
              nextflow run ecseq/dnaseq \
              --input /path/to/input/dir --reference /path/to/genome.fa \
              --FastQC --bamQC --keepReads

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT VERSION AND EXIT
if(params.version){
    println """\
         ===============================================
          E C S E Q - t e m p l a t e   P I P E L I N E
         ===============================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}


// DEFINE PATHS # these are strings which are used to define input Channels,
// but they are specified here as they may be referenced in LOGGING
file1 = file("/path/to/file1.file", checkIfExists: true, glob: false)
file2 = file("${params.some_parameter}", checkIfExists: true, glob: false)



// PRINT STANDARD LOGGING INFO
log.info ""
log.info "         ==============================================="
log.info "          E C S E Q - t e m p l a t e   P I P E L I N E"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ===============================================" }
else {
log.info "         ===============================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         input dir    : ${workflow.profile.tokenize(",").contains("test") ? "-" : "${params.input}"}"
log.info "         reference    : ${params.reference}"
log.info "         output dir   : ${params.output}"
log.info "         mode         : ${params.SE ? "single-end" : "paired-end"}"
log.info ""
log.info "         ==============================================="
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



////////////////////
// STAGE CHANNELS //
////////////////////

/*
 *   Channels are where you define the input for the different
 *    processes which make up a pipeline. Channels indicate
 *    the flow of data, i.e. the "route" that a file will take.
 */

// STAGE BAM FILES FROM TEST PROFILE # this establishes the test data to use with -profile test
if ( workflow.profile.tokenize(",").contains("test") ){

        include { check_test_data } from './lib/functions.nf' params(readPaths: params.readPaths, singleEnd: params.SE)
        INPUT = check_test_data(params.readPaths, params.SE)

} else {

    // STAGE INPUT CHANNELS # this defines the normal input when test profile is not in use
    INPUT = Channel...

}



////////////////////
// BEGIN PIPELINE //
////////////////////

/*
 *   Workflows are where you define how different processes link together. They
 *    may be modularised into "sub-workflows" which must be named eg. 'DNAseq'
 *    and there must always be one MAIN workflow to link them together, which
 *    is always unnamed.
 */

// INCLUDES # here you must give the relevant process files from the lib directory 
include {process_1;process_2;process_3;process_4;process_5} from './lib/process.nf' params(params)

// SUB-WORKFLOWS
workflow 'template' {

    // take the initial Channels and paths
    take:
        INPUT
        file1
        file2

    // here we define the structure of our workflow i.e. how the different processes lead into each other
    // eg. process(input1, input2, etc.)
    // eg. process.out[0], process.out[1], etc.
    // index numbers [0],[1],etc. refer to different outputs defined for processes in process.nf
    // You can forgo the index number if there is only 1 output.
    // ALWAYS PAY ATTENTION TO CARDINALITY!!

    main:
        // process_1 perhaps begins with the INPUT Channel (defined above)
        process_1(INPUT)
        // then process_2 perhaps works on the first output of process_1
        process_2(process_1.out[0])

        // process_3 perhaps works on the individual files only
        process_3(file1,file2)

        // perhaps now we need to combine different process outputs with some kind of Channel operator
        combined_outputs = process_1.out[0].combine(process_3.out)

        // maybe process_4 then uses the combined_channels as well as input files
        process_4(combined,file1,file2)

        // then finally process_5 works on the output of process_4
        process_5(process_4.out[0])

}

// MAIN WORKFLOW 
workflow {

    // call sub-workflows eg. WORKFLOW(Channel1, Channel2, Channel3, etc.)
    main:
        DNAseq(INPUT, file1, file2)

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
