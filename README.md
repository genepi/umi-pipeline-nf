[![Nextflow](https://img.shields.io/badge/nextflow-20.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)

Umi-pipeline-nf
======================

**Umi-pipeline-nf** creates highly accurate single-molecule consensus sequences for unique molecular identifier (UMI)-tagged amplicon data.  
The pipeline can be run for the whole fastq_pass folder of your nanopore run and, per default, outputs the aligned consensus sequences of each UMI cluster in bam file. The optional variant calling creates a vcf file for all variants that are found in the consensus sequences.
umi-pipeline-nf is inspired by a snakemake-based analysis pipeline ([ONT UMI analysis pipeline](https://github.com/nanoporetech/pipeline-umi-amplicon); originally developed by [Karst et al, Nat Biotechnol 18:165–169, 2021](https://www.nature.com/articles/s41592-020-01041-y)). We transferred the pipeline to [Nextflow](https://www.nextflow.io) and included [additional functionalities](#main-adaptations).  

## Workflow

1. Input reads are aligned against a reference genome.
2. The flanking UMI sequences of all reads are extracted.
3. The extracted UMIs are used to cluster the reads.
4. Per cluster, highly accurate consensus sequences are created.
5. The consensus sequences are aligned against the reference sequenced.
6. An optional variant calling step can be performed.

> See the [output documentation](docs/output.md) for a detailed overview of the pipeline and its output files.

## Main Adaptations

* It comes with docker containers making **installation simple, portable** and **results highly reproducible**.
* The pipeline is **optimized for parallelization**.
* **Additional UMI cluster splitting** step to remove admixed UMI clusters.
* Read filtering strategy per UMI cluster was adapted to **preserve the highest quality reads**.
* **Three commonly used variant callers** ([freebayes](https://github.com/freebayes/freebayes), [lofreq](http://csb5.github.io/lofreq/) or [mutserve](https://mitoverse.readthedocs.io/mutserve/mutserve/)) are supported by the pipeline.
* The raw reads can be optionally **subsampled**.
* The raw reads can be **filtered by read length and quality**.
 
> See the [usage documentation](docs/usage.md) for all of the available parameters of the pipeline.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/).

2. Download the pipeline and test it on a [minimal dataset](data/info.txt) with a single command.

```bash
nextflow run genepi/umi-pipeline-nf -r v0.1.0 -profile test,docker
```

3. Start running your own analysis!  
3.1 Download and adapt the config/custom.config with paths to your data (relative and absolute paths possible).

```bash
nextflow run genepi/umi-pipeline-nf -r v0.1.0 -c <custom.config> -profile docker 
```


### Credits

The pipeline was written by ([@StephanAmstler](https://github.com/AmstlerStephan)).  
Nextflow template pipeline: [EcSeq](https://github.com/ecSeq).  
Original snakemake-based pipeline: [nanoporetech/pipeline-umi-amplicon](https://github.com/nanoporetech/pipeline-umi-amplicon).  
Original workflow: [SorenKarst/longread_umi](https://github.com/SorenKarst/longread_umi).
