#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { validateParameters; paramsHelp; paramsSummaryLog  } from 'plugin/nf-validation'
include { UMI_PIPELINE                                      } from './workflows/umi-pipeline.nf'

workflow {

     if (params.help) {
          def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
          def String command = "nextflow run genepi/${workflow.manifest.name} -r v0.2.1 -profile test,docker"
          log.info paramsHelp(command) + citation
          exit 0
     }

     // Validate input parameters
     if (params.validate_params) {
          validateParameters()
     }

     // Print summary of supplied parameters
     log.info paramsSummaryLog(workflow)
     
     // Run the workflow
     UMI_PIPELINE()
}
