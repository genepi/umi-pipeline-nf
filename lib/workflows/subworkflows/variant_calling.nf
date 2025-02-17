include {LOFREQ as LOFREQ_CONSENSUS; LOFREQ as LOFREQ_FINAL_CONSENSUS} from '../../processes/lofreq.nf'
include {MUTSERVE as MUTSERVE_CONSENSUS; MUTSERVE as MUTSERVE_FINAL_CONSENSUS} from '../../processes/mutserve.nf'
include {FREEBAYES as FREEBAYES_CONSENSUS; FREEBAYES as FREEBAYES_FINAL_CONSENSUS} from '../../processes/freebayes.nf'


workflow VARIANT_CALLING {
    take:
        consensus_bam
        final_consensus_bam
        consensus
        final_consensus
        reference 
        reference_fai
        bed 

    main:
        if( params.call_variants ){
            if( params.variant_caller == "lofreq" ){
                LOFREQ_CONSENSUS( consensus_bam, consensus, reference, reference_fai )
                LOFREQ_FINAL_CONSENSUS( final_consensus_bam, final_consensus, reference, reference_fai )
            }else if( params.variant_caller == "mutserve"){
                MUTSERVE_CONSENSUS( consensus_bam, consensus, bed, reference, reference_fai )
                MUTSERVE_FINAL_CONSENSUS( final_consensus_bam, final_consensus, bed, reference, reference_fai )
            }else if( params.variant_caller == "freebayes"){
                FREEBAYES_CONSENSUS( consensus_bam, consensus, reference, reference_fai )
                FREEBAYES_FINAL_CONSENSUS( final_consensus_bam, final_consensus, reference, reference_fai )
            }else{
                exit 1, "${params.variant_caller} is not a valid option. \nPossible variant caller are <lofreq/mutserve/freebayes>"
            
            }
        }      
}