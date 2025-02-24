import nextflow.Nextflow

class WorkflowMain {

    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis, please cite:\n\n" +
            "Nanopore sequencing with unique molecular identifiers enables accurate mutation analysis " +
            "and haplotyping in the complex lipoprotein(a) KIV-2 VNTR. Genome Med 16, 117 (2024). " +
            "https://doi.org/10.1186/s13073-024-01391-8";
    }

    public static void validate(params) {
        def requiredParams = [
            'input', 'reference', 'reference_fai', 'bed', 'output'
        ]

        requiredParams.each {
            param ->
            if (params[param] == null) {
                Nextflow.error("Parameter ${param} is required to run the pipeline.")
            }
        }

    }

    public static String version(workflow) {
        return "" +
            "          ===================                          \n" +
            "            ${workflow.manifest.name.toUpperCase()}    \n" +
            "          ===================                          \n" +
            "            version ${workflow.manifest.version}       \n";
    }

    public static String onComplete(workflow, baseDir, params) {
        // run a small clean-up script to remove "work" directory after successful completion 
        if (workflow.success) {
            if (params.live){
                def stop = new File("${params.input}/barcode_continue/")
                stop.deleteDir()
            }
            ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() 
        }
    }
}