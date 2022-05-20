# ecSeq-DNAseq Output
This document describes the output produced by the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Read Trimming](#read-trimming) - read trimming with cutadapt
* [Read Alignment](#read-alignment) - mapping trimmed reads with bowtie2
* [Quality Control](#quality-control) - generating FastQC and bamQC reports
* [Pipeline Info](#pipeline-info) - reports from nextflow about the pipeline run

## Read Trimming
Input reads will be trimmed with cutadapt, and output to a new directory if specified by the user.

**Output directory: `./trimming`**


## Read Alignment
Depending on which options are specified to the pipeline, trimmed reads will be aligned with bowtie2.

**Output directory: `./mapping`**


## Quality Control
Following trimming, the pipeline will generate FastQC reports for each new set of reads. Following alignment, the pipeline will generate bamQC reports for each BAM file.

**Output directory: `./trimming`**
**Output directory: `./mapping/<sample>`**


## Pipeline Info
Nextflow has several built-in reporting tools that give information about the pipeline run.

**Output directory: `./`**

* `dag.svg`
  * DAG graph giving a diagrammatic view of the pipeline run.
  * NB: If [Graphviz](http://www.graphviz.org/) was not installed when running the pipeline, this file will be in [DOT format](http://www.graphviz.org/content/dot-language) instead of SVG.
* `report.html`
  * Nextflow report describing parameters, computational resource usage and task bash commands used.
* `timeline.html`
  * A waterfall timeline plot showing the running times of the workflow tasks.
* `trace.txt`
  * A text file with machine-readable statistics about every task executed in the pipeline.