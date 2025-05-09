    include {LOFREQ} from '../modules/local/variant_calling/lofreq.nf'
    include {MUTSERVE} from '../modules/local/variant_calling/mutserve.nf'
    include {FREEBAYES} from '../modules/local/variant_calling/freebayes.nf'

    workflow VARIANT_CALLING {
        take:
            bam
            type
            reference 
            reference_fai
            bed 

        main:
            if( params.variant_caller == "lofreq" ){
                LOFREQ( bam, type, reference, reference_fai )
            }else if( params.variant_caller == "mutserve"){
                MUTSERVE( bam, type, bed, reference, reference_fai )
            }else if( params.variant_caller == "freebayes"){
                FREEBAYES( bam, type, reference, reference_fai )
            }else{
                exit 1, "${params.variant_caller} is not a valid option. \nPossible variant caller are <lofreq/mutserve/freebayes>"            
            }

    }