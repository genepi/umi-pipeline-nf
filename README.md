[![Nextflow](https://img.shields.io/badge/nextflow-25.04.8-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![CI Tests](https://github.com/genepi/umi-pipeline-nf/actions/workflows/ci-tests.yml/badge.svg)](https://github.com/genepi/umi-pipeline-nf/actions/workflows/ci-tests.yml)

Umi-pipeline-nf
======================

**Umi-pipeline-nf** creates highly accurate single-molecule consensus sequences for unique molecular identifier (UMI)-tagged amplicons from nanopore sequencing data.
The pipeline processes FastQ files (typically from the `fastq_pass` folder of your nanopore run) and outputs high-quality aligned consensus sequences in BAM format for each UMI cluster. The optional variant calling creates a vcf file for all variants that are found in the consensus sequences.
The newest version of the pipeline supports live analysis of the clusters during sequencing and seamless polishing of the clusters as soon as enough clusters are found.

Umi-pipeline-nf is conceptually based on the Snakemake-based UMI analysis pipeline by ONT ([pipeline-umi-amplicon](https://github.com/nanoporetech/pipeline-umi-amplicon), which relies on the workflow developed by [Karst et al, Nat Biotechnol 18:165–169, 2021](https://www.nature.com/articles/s41592-020-01041-y)). We newly implemented umi-pipeline-nf in [Nextflow](https://www.nextflow.io) and performed a complete redesign and comprehensive functional extension. Umi-pipeline-nf is a fundamentally restructured, extensively optimized and extended pipeline to analyse UMI-tagged nanopore data. It incorporates new modules, GPU-accelerated polishing, two different polishing strategies, real-time sequencing integration, enhanced scalability, simple portability and improved usability.

![Workflow](docs/images/umi-pipeline-nf_metro-map.jpg)

## Workflow

The pipeline is organized into four main subworkflows, each with its own processing steps and outputs:

1. **LIVE UMI PROCESSING**
   - **Purpose:** Real-time processing of raw FastQ files.
   - **Steps:**
     - Merge and filter raw FastQ files.
     - Align reads to the reference genome.
     - Extract UMI sequences.
     - Cluster UMI-tagged reads.
   - **Outputs (in verbose):**
     - Processed UMI clusters are passed on to later stages.
     - Raw alignment files (e.g., in `<output>/<barcodeXX>/raw/align/` or `<output>/<barcodeXX>/<target>/fastq_filtered/raw/`).
     - Filtered FastQ files and clustering statistics.

    **To stop the pipeline when it's in live mode, create a CONTINUE file in the output directory:**
    ```bash
    touch <output>/CONTINUE
    ```
    **Note**
    Nextflow needs write permission in the output directory of MinKNOW -> Add nextflow to the minknow user group (Root rights required to change permissions).
    ```bash
    sudo usermod -aG minknow nextflow
    ```

2. **OFFLINE UMI PROCESSING**
   - **Purpose:** Batch processing with an optional subsampling step.
   - **Steps:**
     - Merge and filter FastQ files.
     - Optionally subsample the merged reads.
     - Perform alignment, UMI extraction, and clustering similar to LIVE processing.
   - **Outputs (in verbose):**
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

   ### Polishing Strategies

   Two complementary polishing strategies are implemented in the pipeline:

   **1. POA-based polishing (graph-based consensus):**
   - Implemented via Medaka’s *smolecule* workflow.
   - All reads within a UMI cluster are aligned to each other using a Partial Order Alignment (POA) algorithm, which represents bases as nodes in a directed acyclic graph (DAG).
   - This captures substitutions, insertions, and deletions efficiently while preserving positional context.
   - The POA alignment is used to generate an initial consensus, which is then refined with Medaka’s neural network.
   - Terminal UMI sequences are retained, enabling optional second-round polishing.
   - **Strengths:** Highly accurate, effective for targets lacking a close reference or with structural variation.
   - **Trade-offs:** Computationally more intensive in terms of time and memory.

   **2. Reference-based polishing (fast reference-aligned consensus):**
   - Reads are directly aligned to the reference sequence using minimap2.
   - Consensus is generated with Medaka’s consensus and stitch modules.
   - By bypassing the POA graph construction, this approach reduces RAM usage (up to 5×) and runtime (up to 30×).
   - **Strengths:** Highly efficient and scalable, suitable when a sufficiently similar reference is available.
   - **Trade-offs:** Less robust for targets with large indels, rearrangements, or high sequence heterogeneity.

4. **VARIANT CALLING**
   - **Purpose:** Identify genetic variants from the consensus data.
   - **Steps:**
     - Perform variant calling using one of the supported callers: [freebayes](https://github.com/freebayes/freebayes), [lofreq](http://csb5.github.io/lofreq/), or [mutserve](https://mitoverse.readthedocs.io/mutserve/mutserve/).
   - **Outputs:**
     - VCF files with variant calls (e.g., in `<output>/<barcodeXX>/<target>/<freebayes/mutserve/lofreq>/`).

> See the [output documentation](docs/output.md) for a detailed overview of the pipeline outputs and directory structure.

## Live Demo

The **live_demo** is a self-contained example that simulates a nanopore sequencing run with UMI-tagged reads. It allows users to **test the live mode of the pipeline without actual sequencing data** and observe the real-time processing, clustering, and polishing of reads.

### Purpose

- Demonstrates the **LIVE UMI PROCESSING** workflow in action.
- Helps users understand how the pipeline handles **incoming reads**, performs **UMI clustering**, and polishes consensus sequences on-the-fly.
- Useful for training, testing, or debugging the pipeline without connecting to a sequencing device.

### How It Works

1. **Simulated Read Generation**
   - A small input FASTQ file is repeatedly copied into barcode-specific folders (`barcode01`, `barcode02`, …) at regular intervals.
   - This mimics the arrival of reads from a live nanopore run.

2. **Live Pipeline Execution**
   - The pipeline is started in background mode using the `--live` flag.
   - As new reads appear in the input directories, the pipeline automatically processes them, extracts UMIs, clusters reads, and updates consensus sequences.

3. **Stopping the Demo**
   - To gracefully stop the pipeline, create a `CONTINUE` file in the output directory:
     ```bash
     touch <output>/CONTINUE
     ```
   - The pipeline will finish processing the currently queued reads and terminate.

4. **Cleanup**
   - After the pipeline finishes, temporary input and output files can be safely removed using the included cleanup script.

### Running the Live Demo (Make sure Nextflow is installed on your system)

```bash
git clone https://github.com/genepi/umi-pipeline-nf.git
cd umi-pipeline-nf
bash live_demo/run_demo.sh --iterations 4 --interval 20
```

- `--iterations`: Number of times the input FASTQ is copied to simulate new reads (default: 5).
- `--interval`: Time in seconds between iterations (default: 30).

> The demo runs a complete simulated live sequencing session, showing how **UMI extraction, clustering, and consensus polishing** happen in real time.

### Output

- Simulated input reads are placed in `live/data/fastq_pass/barcode*/`.
- Processed consensus sequences and clustering statistics appear in `live/output/barcode*/`.
- Logs and intermediate files are visible if `--verbose` is enabled.

## Main features and differences to other UMI pipelines

* It comes with a docker/singularity container making **installation simple, easy to use on clusters** and **results highly reproducible**.
* The pipeline is **optimized for parallelization**.
* **Additional UMI cluster splitting** step to remove admixed UMI clusters.
* Read filtering strategy per UMI cluster was adapted to **preserve the highest quality reads**.
* **Three commonly used variant callers** ([freebayes](https://github.com/freebayes/freebayes), [lofreq](http://csb5.github.io/lofreq/) or [mutserve](https://mitoverse.readthedocs.io/mutserve/mutserve/)) are supported by the pipeline.
* The raw reads can be optionally **subsampled**.
* The raw reads can be **filtered by read length and quality**.
* **GPU acceleration for cluster polishing by Medaka** is available when using the `docker` profile. Tested with an NVIDIA RTX 4080 SUPER GPU (16 GB).
* Two polishing strategies are supported:
  * POA-based polishing for maximum accuracy, especially in regions with high variation or no close reference.
  * Reference-based polishing for faster, memory-efficient consensus building when a reliable reference is available.
* Allows **multi line bed files** to run the pipeline for several targets at once.
* Supports **live analysis of the clusters during sequencing and seamless polishing** of the clusters as soon as enough clusters are found.
* Use **verbose parameter** to output all intermediate files.
* **Optional context primer** help improving the UMI extraction for non-unique UMI sequences.

To see all available parameters run
```bash
nextflow run genepi/umi-pipeline-nf -r v1.0.0 --help
```

## Quick Start

1. Install [`nextflow`](https://www.nextflow.io/).

2. Download the pipeline and test it on a [minimal dataset](data/info.txt) with a single command.

```bash
nextflow run genepi/umi-pipeline-nf -r v1.0.0 -profile test,docker
```

3. Start running your own analysis!
3.1 Download and adapt the config/custom.config with paths to your data (relative and absolute paths possible).

```bash
nextflow run genepi/umi-pipeline-nf -r v1.0.0 -c <custom.config> -profile custom,<docker,singularity>
```

## Citation

If you use the pipeline please cite [our Paper](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-024-01391-8):

Amstler S, Streiter G, Pfurtscheller C, Forer L, Di Maio S, Weissensteiner H, Paulweber B, Schoenherr S, Kronenberg F, Coassin S. Nanopore sequencing with unique molecular identifiers enables accurate mutation analysis and haplotyping in the complex lipoprotein(a) KIV-2 VNTR. Genome Med 16, 117 (2024). https://doi.org/10.1186/s13073-024-01391-8

### Credits

The pipeline was written by [@StephanAmstler](https://github.com/AmstlerStephan).
Nextflow template pipeline: [EcSeq](https://github.com/ecSeq).
Snakemake-based ONT pipeline for UMI nanopore sequencing analysis: [nanoporetech/pipeline-umi-amplicon](https://github.com/nanoporetech/pipeline-umi-amplicon).
UMI-corrected nanopore sequencing analysis first shown by: [SorenKarst/longread_umi](https://github.com/SorenKarst/longread_umi).

