[<img width="200" align="right" src="docs/images/ecseq.jpg">](https://www.ecseq.com)
[![Nextflow](https://img.shields.io/badge/nextflow-20.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/ecseq/dnaseq.svg)](https://hub.docker.com/r/ecseq/dnaseq)

umi-pipeline-nf Pipeline
======================

**umi-pipeline-nf** is based on a [snakemake pipeline](https://github.com/nanoporetech/pipeline-umi-amplicon) provided by [Oxford Nanopore Technologies (ONT)](https://nanoporetech.com/). To increase efficiency and usability the pipeline was transferred to [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Overview
`umi-pipeline-nf` creates highly accurate single-molecule consensus sequences based on amplicon data tagged by unique molecular identifiers (UMIs). The pipeline can be run for the whole fastq_pass folder of your nanopore run and per default, the output are the aligned consensus sequences in bam file format. 
Additional flags can be set to perform a variant calling ( [freebayes](https://github.com/freebayes/freebayes), [lofreq](http://csb5.github.io/lofreq/) or [mutserve](https://mitoverse.readthedocs.io/mutserve/mutserve/) )

> See the [output documentation](docs/output.md) for more details of the results.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/)

2. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run genepi/umi-pipeline-nf -profile test,docker
```

3. Start running your own analysis!
3.1 Download and adapt the config/custom.config with paths to your data (relative and absolute paths possible)

```bash
nextflow run genepi/umi-pipeline-nf -r default_clustering_strategy -latest -c <custom.config> -profile docker 
```

> See the [usage documentation](docs/usage.md) for all of the available options when running the pipeline.


### Credits

These scripts were originally written for use by [GENEPI](https://genepi.i-med.ac.at/), by ([@StephanAmstler](https://github.com/AmstlerStephan)).
The version corresponds to the latest update of the master thesis. 
