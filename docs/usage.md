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
    * [`--variant_caller`](#--variant_caller-arg-required)
* [Advanced parameters](#advanced-parameters)
    * [`--min_reads_per_barcode`](#--min_reads_per_barcode-arg)
    * [`--umi_errors`](#--umi_errors-arg)
    * [`--min_reads_per_cluster`](#--min_reads_per_cluster-arg)
    * [`--max_reads_per_cluster`](#--max_reads_per_cluster-arg)
    * [`--filter_strategy_clusters`](#--filter_strategy_clusters-arg)
    * [`--output_format`](#--output_format-arg)
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
    * [`-with-docker`](#-with-docker-arg)
    * [`-with-singularity`](#-with-singularity-arg)
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
nextflow run ecSeq/DNAseq [OPTIONS]
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
Name the output directory to save the results. [default: "null"]

### `--reference <ARG>` [REQUIRED]
Specify the path to the reference genome in fasta format. [default: "null"]

### `--reference_fai <ARG>` [REQUIRED]
Specify the path to the reference index file. [default: "null"]

### `--bed <ARG>` [REQUIRED]
Specify the path to the bed file containing the target name, start and end position. [default: "null"]

## Advanced Parameters

### `--output_format <ARG>` [fasta | fastq]
Specify the output format until the cluster filtering step.  
Note: Only the "fastq" option can be used to filter reads by quality and requires FASTQ input reads. [default: "fasta"]  

### `--filter_strategy_clusters <ARG>` [random | quality]
Specify filter strategy to downsample cluster above the maximal reads per cluster parameter.  
Note: "quality" filter strategy only goes in hand with "fastq" output format. 


## Additional Parameters

### `--debug`
Specify in order to prevent Nextflow from clearing the work dir cache following a successful pipeline completion. [default: off]

### `--version`
When called with `nextflow run ecseq/dnaseq --version` this will display the pipeline version and quit.

### `--help`
When called with `nextflow run ecseq/dnaseq --help` this will display the parameter options and quit.

## Software Dependencies

There are different ways to provide the required software dependencies for the pipeline. The recommended method is to use the Conda, Docker or Singularity profiles as provided by the pipeline. 

### `-profile <ARG>`
Use this parameter to choose a preset configuration profile. Profiles available with the pipeline are:

* `standard`
    * The default profile, used if `-profile` is not specified.
    * Uses sensible resource allocation for , runs using the `local` executor (native system calls) and expects all software to be installed and available on the `$PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles below.
* `conda`
    * Builds a conda environment from the environment.yml file provided by the pipeline
    * Requires conda to be installed on your system.
* `docker`
    * Launches a docker image pulled from ecseq/dnaseq
    * Requires docker to be installed on your system. 
* `singularity`
    * Launches a singularity image pulled from ecseq/dnaseq
    * Requires singularity to be installed on your system.
* `custom`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config for process resource allocation.

If you wish to provide your own package containers it is possible to do so by setting the `standard` or `custom` profile, and then providing your custom package with the command line flags below. These are not required with the the other profiles.

### `-with-conda <ARG>`
Flag to enable conda. You can provide either a pre-built environment or a *.yaml file.

### `-with-docker <ARG>`
Flag to enable docker. The image will automatically be pulled from Dockerhub.

### `-with-singularity <ARG>`
Flag to enable use of singularity. The image will automatically be pulled from the internet. If running offline, follow the option with the full path to the image file.

## Other command line parameters

### `-work-dir <ARG>`
Specify the path to a custom work directory for the pipeline to run with (eg. on a scratch directory)

### `-params-file <ARG>`
Provide a file with specified parameters to avoid typing them out on the command line. This is useful for carrying out repeated analyses. A template params file [`assets/params.config`](../assets/params.config) has been made available in the pipeline repository.

### `-config <ARG>`
Provide a custom config file for adapting the pipeline to run on your own computing infrastructure. A template config file [`assets/custom.config`](../assets/custom.config) has been made available in the pipeline repository. This file can be used as a boilerplate for building your own custom config.

### `-resume [<ARG>]`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. Give a specific pipeline name as an argument to resume it, otherwise Nextflow will resume the most recent. NOTE: This will not work if the specified run finished successfully and the cache was automatically cleared. (see: [`--debug`](#--debug))

### `-name <ARG>`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
