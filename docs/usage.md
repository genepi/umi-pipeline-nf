# ecSeq-DNAseq Usage
This document describes the parameter options used by the pipeline.

* [Running the pipeline](#running-the-pipeline)
* [Inputs and outputs](#inputs-and-outputs)
    * [`--input`](#--input-arg-required)
    * [`--reference`](#--reference-arg-required)
    * [`--reference_fai`](#--reference_fai-arg-required)
    * [`--bed`](#--bed-arg-required)
    * [`--output`](#--output-arg-required)
* [Workflow modifying parameters](#workflow-modifying-parameters)
    * [`--subsampling`](#--subsampling)
    * [`--call_variants`](#--call_variants)
* [Read filtering parameters](#read-filtering-parameters)
    * [`--min_read_length`](#--min_read_length-arg)
    * [`--min_qscore`](#--min_qscore-arg)
* [Subsampling parameters](#trimming-parameters)
    * [`--subsampling_seed`](#--subsampling_seed-arg)
    * [`--subsampling_readnumber`](#--subsampling_readnumber-arg)
* [Variant calling parameters](#variant-calling-parameters)
    * [`--variant_caller`](#--variant_caller-arg-required-lofreq--freebayes--mutserve)
* [Advanced parameters](#advanced-parameters)
    * [`--min_reads_per_barcode`](#--min_reads_per_barcode-arg)
    * [`--umi_errors`](#--umi_errors-arg)
    * [`--min_reads_per_cluster`](#--min_reads_per_cluster-arg)
    * [`--max_reads_per_cluster`](#--max_reads_per_cluster-arg)
    * [`--filter_strategy_clusters`](#--filter_strategy_clusters-arg-random--quality)
    * [`--output_format`](#--output_format-arg-fasta--fastq)
    * [`--min_overlap`](#--min_overlap-arg)
    * [`--include_secondary_reads`](#--include_secdondary_reads-arg)
    * [`--balance_strands`](#--balance_strands-arg)
    * [`--medaka_model`](#--medaka_model-arg)
    * [`--fwd_umi`](#--fwd_umi-arg)
    * [`--rev_umi`](#--rev_umi-arg)
    * [`--min_length`](#--min_length-arg)
    * [`--max_length`](#--max_length-arg)
    * [`--minimap_param`](#--minimap_param-arg)
    * [`--write_reports`](#--write_reports-arg)
    * [`--threads`](#--threads-arg)
* [Additional parameters](#additional-parameters)
    * [`--debug`](#--debug)
    * [`--version`](#--version)
    * [`--help`](#--help)
* [Software dependencies](#software-dependencies)
    * [`-profile`](#-profile-arg)
    * [`-with-conda`](#-with-conda-arg)
* [Other command line parameters](#other-command-line-parameters)
    * [`-work-dir`](#-work-dir-arg)
    * [`-params-file`](#-params-file-arg)
    * [`-config`](#-config-arg)
    * [`-resume`](#-resume-arg)
    * [`-name`](#-name-arg)

## Workflow

## Running the pipeline
The main command for running the pipeline is as follows:

```bash
nextflow run genepi/umi-pipeline-nf [OPTIONS]
```

Note that the pipeline will create files in your working directory:

```bash
work/           # Directory containing the nextflow working files
.nextflow.log   # Log file from Nextflow
.nextflow/      # Nextflow cache and history information
```

## Inputs and Outputs

### `--input <ARG>` [REQUIRED]
Specify the path to the directory containing the demultiplexed reads. [default: "null"]

### `--output <ARG>` [REQUIRED]
Specify the name of the output directory to save the results. [default: "null"]

### `--reference <ARG>` [REQUIRED]
Specify the path to the reference genome in fasta format. [default: "null"]

### `--reference_fai <ARG>` [REQUIRED]
Specify the path to the reference index file. [default: "null"]

### `--bed <ARG>` [REQUIRED]
Specify the path to the bed file containing the target name, start and end position. [default: "null"]

## Workflow Modifying Parameters 

### `--subsampling`
Specify if the raw reads per barcode should be subsampled. [default: false]

### `--call_variants`
Specify if the variants in the final consensus sequences should be called. [default: false]
Note: If set to true, the (`variant caller`)[ must be specified

## Read Filtering Parameters

### `--min_read_length <ARG>`
Specify the minimal read length of the raw reads. [default: 0]

### `--min_qcore <ARG>`
Specify the minimal Q-score of the raw reads. [default: 0]

## Subsampling Parameters

### `--subsampling_seed <ARG>`
Specify the seed for pseudo-random subsampling. [default: 11]

### `--subsampling_readnumber <ARG>`
Specify the number of reads that should be preserved after subsampling. [default: 100,000]

## Variant Calling Parameters

### `--variant_caller <ARG>` [REQUIRED] [lofreq | freebayes | mutserve]
Specify the variant caller that should be used for variant calling.  [default: null]

## Advanced Parameters

### `--min_reads_per_barcode <ARG>`
Specify the minimal number of raw reads per barcode. [default: 1000]  
Note: Barcodes with fewer reads will be filtered before the analysis.

### `--umi_errors <ARG>`
Specify the number of deviating positions per read of each UMI cluster. [default: 3]

### `--min_reads_per_cluster <ARG>`
Specify the minimal number of reads per cluster. [default: 20]
Note: Clusters with fewer reads will be excluded from the analysis. 

### `--max_reads_per_cluster <ARG>`
Specify the maximal number of reads per cluster. [default: 60]
Note: Clusters with more reads will be downsampled  
according to the [`cluster filtering strategy`](#--filter_strategy_clusters-arg-random-quality).

### `--output_format <ARG>` [fasta | fastq]
Specify the output format until the cluster filtering step. [default: "fastq"]  
Note: Only the "fastq" option can be used to filter reads by the quality and requires FASTQ input reads.  

### `--filter_strategy_clusters <ARG>` [random | quality]
Specify filter strategy to downsample cluster above the maximal reads per cluster parameter. [default: "quality"]  
Note: "quality" filter strategy only goes in hand with "fastq" output format. 

### `--write-reports <ARG>`
Specify if reports of the pipeline should be saved. [default: true]

### `--threads <ARG>`
Specify the number of threads that are allocated for the pipeline. [default: all available processors - 1]

### `--min_overlap <ARG>` [0-1]
Specify the minimal overlap of the mapped raw reads to the reference sequence. [default: 0.90]

### `--include_secondary_reads <ARG>`
Specify if secondary mappings should be included in the analysis. [default: false]

### `--balance_strands <ARG>`
Specify if the number of forward and reverse reads per cluster should be equalized. [default: true]

### `--medaka_model <ARG>`
Specify the medaka model that is used for cluster polishing. [default: "r1041_e82_400bps_hac_g615"]  
Note: The models are specific for Chemistry, basecalling algorithm and sequencing speed. 

### `--fwd_umi <ARG>`
Specify the pattern of the forward UMI primer. [default: "TTTVVVVTTVVVVTTVVVVTTVVVVTTT"]

### `--rev_umi <ARG>`
Specify the pattern of the reverse UMI primer. [default: "AAABBBBAABBBBAABBBBAABBBBAAA"]

### `--min_length <ARG>`
Specify the minimal length of the extracted UMI sequences. [default: 40]

### `--min_length <ARG>`
Specify the maximal length of the extracted UMI sequences. [default: 60]

### `--minimap2_param <ARG>`
Specify the minimap2 parameters. [default: "-ax map-ont -k 13]

## Additional Parameters

### `--debug`
Specify in order to prevent Nextflow from clearing the work dir cache following a successful pipeline completion. [default: off]

### `--version`
When called with `nextflow run genepi/umi-pipeline-nf --version` this will display the pipeline version and quit.

### `--help`
When called with `nextflow run genepi/umi-pipeline-nf --help` this will display the parameter options and quit.

## Software Dependencies

### `-profile <ARG>`
Use this parameter to choose a preset configuration profile. Profiles available with the pipeline are:

* `standard`
    * The default profile, is used if `-profile` is not specified.
    * Uses sensible resource allocation for, runs using the `local` executor (native system calls) and expects all software to be installed and available on the `$PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles below.
* `conda`
    * Builds a conda environment from the environment.yml file provided by the pipeline
    * Requires conda to be installed on your system.
* `docker`
    * Launches a docker image pulled from genepi/umi-pipeline-nf (ADAPT)
    * Requires docker to be installed on your system. 
* `singularity`
    * Launches a singularity image pulled from genepi/umi-pipeline-nf (ADAPT)
    * Requires singularity to be installed on your system.
* `custom`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config for process resource allocation.

If you wish to provide your own package containers it is possible to do so by setting the `standard` or `custom` profile and then providing your custom package with the command line flags below. These are not required with the other profiles.

### `-with-conda <ARG>`
Flag to enable conda. You can provide either a pre-built environment or a *.yml file.

## Other command line parameters

### `-work-dir <ARG>`
Specify the path to a custom work directory for the pipeline to run with (eg. on a scratch directory)

### `-params-file <ARG>`
Provide a file with specified parameters to avoid typing them out on the command line. This is useful for carrying out repeated analyses. A template params file [`assets/params.config`](../assets/params.config) has been made available in the pipeline repository.

### `-config <ARG>`
Provide a custom config file for adapting the pipeline to run on your own computing infrastructure. A template config file [`assets/custom.config`](../assets/custom.config) has been made available in the pipeline repository. This file can be used as a boilerplate for building your own custom config.

### `-resume [<ARG>]`
Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where they got previously. Give a specific pipeline name as an argument to resume it, otherwise, Nextflow will resume the most recent. NOTE: This will not work if the specified run finished successfully and the cache was automatically cleared. (see: [`--debug`](#--debug))

### `-name <ARG>`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
