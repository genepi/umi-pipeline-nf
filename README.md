[![Nextflow](https://img.shields.io/badge/nextflow-20.07.1-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)

Umi-pipeline-nf
======================

**Umi-pipeline-nf** creates highly accurate single-molecule consensus sequences for unique molecular identifier (UMI)-tagged amplicons from nanopore sequencing data.  
The pipeline processes FastQ files (typically from the `fastq_pass` folder of your nanopore run) and outputs high-quality aligned consensus sequences in BAM format for each UMI cluster. The optional variant calling creates a vcf file for all variants that are found in the consensus sequences.  
The newest version of the pipeline supports live analysis of the clusters during sequencing and seemless polishing of the clusters as soon as enough clusters are found.

Umi-pipeline-nf originated from a Snakemake-based analysis pipeline ([pipeline-umi-amplicon](https://github.com/nanoporetech/pipeline-umi-amplicon); originally developed by [Karst et al, Nat Biotechnol 18:165–169, 2021](https://www.nature.com/articles/s41592-020-01041-y)). We have migrated the pipeline to [Nextflow](https://www.nextflow.io) and incorporated several optimizations and additional functionalities.


![Workflow](docs/images/umi-pipeline-nf_metro-map.svg)

## Workflow

The pipeline is organized into four main subworkflows, each with its own processing steps and outputs:

1. **LIVE UMI PROCESSING**  
   - **Purpose:** Real-time processing of raw FastQ files.
   - **Steps:**  
     - Merge and filter raw FastQ files.
     - Align reads to the reference genome.
     - Extract UMI sequences.
     - Cluster UMI-tagged reads.
   - **Outputs:**  
     - Processed UMI clusters are passed on to later stages.
     - Raw alignment files (e.g., in `<output>/<barcodeXX>/raw/align/` or `<output>/<barcodeXX>/<target>/fastq_filtered/raw/`).
     - Filtered FastQ files and clustering statistics.

    **To stop the pipeline when it's in live mode, create a CONTINUE file in the output directory:**  
    `touch <output>/CONTINUE`

2. **OFFLINE UMI PROCESSING**  
   - **Purpose:** Batch processing with an optional subsampling step.
   - **Steps:**  
     - Merge and filter FastQ files.
     - Optionally subsample the merged reads.
     - Perform alignment, UMI extraction, and clustering similar to LIVE processing.
   - **Outputs:**  
     - Processed UMI clusters.
     - Alignment and subsampling reports (e.g., in `<output>/<barcodeXX>/raw/subsampling/` and `<output>/<barcodeXX>/<target>/fastq_filtered/raw/`).

3. **UMI POLISHING**  
   - **Purpose:** Refine UMI clusters to generate high-quality consensus sequences.
   - **Steps:**  
     - Polish clusters using medaka.
     - Realign consensus sequences to the reference genome.
     - Re-extract and re-cluster UMIs from consensus reads.
     - Parse final consensus clusters.
   - **Outputs:**  
     - Consensus BAM and FastQ files (e.g., in `<output>/<barcodeXX>/<target>/align/consensus/` and `<output>/<barcodeXX>/<target>/fastq/consensus/`).
     - Polishing logs and detailed cluster statistics.

4. **VARIANT CALLING**  
   - **Purpose:** Identify genetic variants from the consensus data.
   - **Steps:**  
     - Perform variant calling using one of the supported callers: [freebayes](https://github.com/freebayes/freebayes), [lofreq](http://csb5.github.io/lofreq/), or [mutserve](https://mitoverse.readthedocs.io/mutserve/mutserve/).
   - **Outputs:**  
     - VCF files with variant calls (e.g., in `<output>/<barcodeXX>/<target>/<freebayes/mutserve/lofreq>/`).

> See the [output documentation](docs/output.md) for a detailed overview of the pipeline outputs and directory structure.

## Main Adaptations

* It comes with a docker/singularity container making **installation simple, easy to use on clusters** and **results highly reproducible**.
* The pipeline is **optimized for parallelization**.
* **Additional UMI cluster splitting** step to remove admixed UMI clusters.
* Read filtering strategy per UMI cluster was adapted to **preserve the highest quality reads**.
* **Three commonly used variant callers** ([freebayes](https://github.com/freebayes/freebayes), [lofreq](http://csb5.github.io/lofreq/) or [mutserve](https://mitoverse.readthedocs.io/mutserve/mutserve/)) are supported by the pipeline.
* The raw reads can be optionally **subsampled**.
* The raw reads can be **filtered by read length and quality**.
* **GPU acceleration for cluster polishing by Medaka** is available when using the `docker` profile. The GPU driver, [nvidia-toolkit](https://developer.nvidia.com/cuda-toolkit), and [nvidia-container-toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) must be installed. Tested with an RTX 4080 SUPER GPU (16 GB). Note: GPU acceleration is not compatible with cluster profiles.
* Allows multi line bed files to run the pipeline for several targets at once.
 
> See the [usage documentation](docs/usage.md) for all of the available parameters of the pipeline.

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/).

2. Download the pipeline and test it on a [minimal dataset](data/info.txt) with a single command.

```bash
nextflow run genepi/umi-pipeline-nf -r v1.0.1 -profile test,docker
```

3. Start running your own analysis!  
3.1 Download and adapt the config/custom.config with paths to your data (relative and absolute paths possible).

```bash
nextflow run genepi/umi-pipeline-nf -r v1.0.1 -c <custom.config> -profile custom,<docker,singularity> 
```

## Citation 

If you use the pipeline please cite [our Paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01391-8):

Amstler S, Streiter G, Pfurtscheller C, Forer L, Di Maio S, Weissensteiner H, Paulweber B, Schoenherr S, Kronenberg F, Coassin S. Nanopore sequencing with unique molecular identifiers enables accurate mutation analysis and haplotyping in the complex lipoprotein(a) KIV-2 VNTR. Genome Med 16, 117 (2024). https://doi.org/10.1186/s13073-024-01391-8


### Credits

The pipeline was written by [@StephanAmstler](https://github.com/AmstlerStephan).  
Nextflow template pipeline: [EcSeq](https://github.com/ecSeq).  
Snakemake-based ONT pipeline for UMI nanopore sequencing analysis: [nanoporetech/pipeline-umi-amplicon](https://github.com/nanoporetech/pipeline-umi-amplicon).  
UMI-corrected nanopore sequencing analysis first shown by: [SorenKarst/longread_umi](https://github.com/SorenKarst/longread_umi).
