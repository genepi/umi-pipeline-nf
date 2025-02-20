# UMI Pipeline Subworkflow Output Overview

This document details the outputs produced by each subworkflow of the UMI analysis pipeline. The pipeline is split into distinct subworkflows, each handling a specific part of the processing—from raw input to variant calling. This guide explains the outputs, the processing steps/modules used, and how these outputs are organized.

> The overall workflow definitions are provided in the merged workflows file. :contentReference[oaicite:0]{index=0}  
> The module structure (with locations and script names) is detailed in the module structure file. :contentReference[oaicite:1]{index=1}  
> The output directory layout is shown in the output structure file. :contentReference[oaicite:2]{index=2}

## Table of Contents

1. [LIVE UMI PROCESSING](#live-umi-processing)
2. [OFFLINE UMI PROCESSING](#offline-umi-processing)
3. [UMI POLISHING](#umi-polishing)
4. [VARIANT CALLING](#variant-calling)
5. [Overall Directory Structure](#overall-directory-structure)
6. [Conclusion](#conclusion)

---

## 1. LIVE UMI PROCESSING

**Purpose:**  
This subworkflow is designed to process raw FastQ files in real time. It watches for new input files and a "continue" signal, processes the files in chunks, and produces early clustering outputs.

**Key Steps & Modules:**

- **CONTINUE_PIPELINE:**  
  Monitors a specific file (e.g., a file with "continue" in its name) to trigger processing.

- **MERGE_FASTQ:**  
  Reads FastQ files from barcode directories, splits them into chunks (using a parameterized chunk size), and merges the chunks.

- **MAP_READS:**  
  Maps the merged (or chunked) FastQ files to a reference genome, creating raw alignment BAM files.

- **SPLIT_READS:**  
  Splits the mapped BAM files into subsets (often to isolate regions of interest).

- **DETECT_UMI_FASTQ:**  
  Extracts UMI (Unique Molecular Identifier) information from the split FastQ files.

- **CLUSTER:**  
  Clusters reads based on UMI tags. A filtering step ensures only clusters with a minimum number of reads (set by a parameter) are kept.

- **REFORMAT_FILTER_CLUSTER:**  
  Reformats and further filters the clusters to prepare them for downstream analysis.

- **CLUSTER_STATS_LIVE & SUMMARY_CLUSTER_STATS:**  
  Generate real-time and summary statistics for the clustering process.

**Outputs:**

- **Processed UMIs:**  
  The final output is a channel called `processed_umis` which aggregates UMI clusters for each sample. This data is forwarded to later stages (polishing and variant calling).

- **Raw Alignment Files:**  
  Under each barcode (e.g., `barcode02/align/raw`), you’ll find raw BAM files (e.g., `filtered_no_umis.1.bam`) and their corresponding index files.

- **Filtered FastQ Files:**  
  Directories like `barcode02/fastq_filtered/raw` contain the FastQ files after initial filtering.

- **Statistics:**  
  Subdirectories (e.g., `barcode02/stats/raw`) contain files such as TSV reports summarizing the performance of read splitting and clustering.

> For additional details on the processing logic, see the merged workflows file. :contentReference[oaicite:3]{index=3}

---

## 2. OFFLINE UMI PROCESSING

**Purpose:**  
The offline subworkflow processes the input FastQ files in a batch mode. Unlike the live mode, it does not wait for a continue signal and may optionally perform subsampling.

**Key Steps & Modules:**

- **MERGE_FASTQ:**  
  Similar to the live workflow, FastQ files from each barcode are merged.

- **SUBSAMPLING (optional):**  
  If enabled (via `params.subsampling`), the merged FastQ is subsampled before proceeding.

- **MAP_READS:**  
  Maps the (merged or subsampled) FastQ files to the reference genome, generating BAM files.

- **SPLIT_READS:**  
  Splits the mapped reads to isolate target regions or consensus areas.

- **DETECT_UMI_FASTQ:**  
  Extracts UMI information from the split reads.

- **CLUSTER:**  
  Clusters reads based on their UMI information and filters out clusters that do not meet the minimum read threshold.

- **REFORMAT_FILTER_CLUSTER:**  
  Reformats the cluster outputs to produce a standardized `processed_umis` set.

**Outputs:**

- **Processed UMIs:**  
  As with the live mode, the final output is a set of `processed_umis` which is passed to the polishing stage.

- **Alignment Files:**  
  BAM files (from mapping) are stored in directories like `barcodeXX/align/raw`.

- **Filtered FastQ and Clustering Results:**  
  The outputs include filtered FastQ files (in `fastq_filtered/raw`) and clustering outputs. Clustering results include FASTA files with consensus sequences in `clustering/consensus` and raw clustering outputs in `clustering/raw`.

- **Statistical Reports:**  
  TSV files (e.g., `split_cluster_stats.tsv`) in directories such as `stats/raw` summarize key metrics from this processing stage.

> Module-level details and script locations are provided in the module structure file. :contentReference[oaicite:4]{index=4}

---

## 3. UMI POLISHING

**Purpose:**  
This subworkflow refines the UMI clusters from the previous steps to generate high-quality consensus sequences. It glues together clusters, performs polishing, and produces consensus outputs that will serve as input for variant calling.

**Key Steps & Modules:**

- **GLUE_CLUSTERS:**  
  Combines clusters when needed, ensuring that clusters belonging to the same sample are properly aggregated.

- **POLISH_CLUSTER:**  
  Polishes the consensus sequences from the clusters, improving sequence accuracy.

- **MERGE_CONSENSUS_FASTQ:**  
  Merges polished clusters into a single consensus FastQ file per sample.

- **FILTER_CONSENSUS_FASTQ (conditional):**  
  Filters the consensus FastQ file if the output format is set to FastQ.

- **MAP_CONSENSUS / MAP_FINAL_CONSENSUS:**  
  Maps both the consensus and the final polished FastQ files to the reference genome, generating consensus BAM files.

- **DETECT_UMI_CONSENSUS_FASTQ & CLUSTER_CONSENSUS:**  
  Further process consensus sequences to extract UMI data and perform a secondary clustering if needed.

- **REFORMAT_CONSENSUS_CLUSTER:**  
  Reformats consensus cluster outputs into final, standardized FASTA/FASTQ files.

**Outputs:**

- **Consensus BAM Files:**  
  - `consensus_bam`: Generated from mapping the consensus FastQ file.  
  - `final_consensus_bam`: Produced after the final mapping of the polished consensus.

- **Consensus FastQ Files:**  
  Located in directories such as `barcodeXX/fastq/consensus` (e.g., `masked_consensus.fastq` and `merged_consensus.fastq`) and the final FastQ files in `barcodeXX/fastq/final`.

- **Polishing Logs and Reports:**  
  Log files (e.g., `barcodeXX_combined_clusters_1_smolecule.log`) and other reports from the polishing process are stored in `polishing/consensus`.

> For further insight into how consensus sequences are derived, refer to the merged workflows and module structure files. :contentReference[oaicite:5]{index=5} :contentReference[oaicite:6]{index=6}

---

## 4. VARIANT CALLING

**Purpose:**  
The final subworkflow is dedicated to variant calling. It takes the consensus BAM files from the polishing stage and identifies genetic variants using one of several variant callers, based on user parameters.

**Key Steps & Modules:**

- **Variant Caller Selection:**  
  The workflow supports multiple variant callers:
  - **LOFREQ:** Processes both consensus and final BAM files.
  - **MUTSERVE:** Uses a bed file and reference to call variants.
  - **FREEBAYES:** Calls variants based on consensus data.
  
- **Input BAM Files:**  
  Both `consensus_bam` and `final_consensus_bam` are used as inputs to the variant calling process.

- **Output Generation:**  
  The selected variant caller produces VCF files that list the genetic variants.

**Outputs:**

- **VCF Files:**  
  For each sample, VCF files are generated. For example, when using FREEBAYES, you will see:
  - Under `freebayes/barcode03/consensus`: A file such as `consensus.vcf`
  - Under `freebayes/barcode03/final`: A file such as `final.vcf`
  
- **Caller-Specific Logs/Reports:**  
  Depending on the caller, additional reports or logs may be created detailing variant calling metrics.

> The logic for variant caller selection and the corresponding outputs is explained in the merged workflows file. :contentReference[oaicite:7]{index=7}

---

## Overall Directory Structure

Below is an overview of the directory structure created after a complete pipeline run. This structure illustrates where outputs from each subworkflow are stored.

output_directory/
├── barcodeXX/
│   ├── align/
│   │   ├── raw/           # Raw alignments (e.g., filtered_umi_KIV2A.1.bam)
│   │   ├── consensus/     # Consensus alignments (e.g., masked_consensus.bam)
│   │   └── final/         # Final alignments (e.g., final.bam)
│   ├── clustering/
│   │   ├── raw/           # Raw clustering outputs
│   │   └── consensus/     # Consensus clustering outputs (e.g., consensus.fasta)
│   ├── fastq/
│   │   ├── consensus/     # Consensus FastQ files (e.g., masked_consensus.fastq)
│   │   └── final/         # Final FastQ files (e.g., final.fastq)
│   ├── fastq_filtered/
│   │   └── raw/           # Filtered FastQ files
│   ├── fastq_umi/
│   │   └── consensus/     # UMI-specific consensus FastQ files
│   ├── polishing/
│   │   └── consensus/     # Polishing outputs and logs
│   └── stats/
│       ├── consensus/     # Consensus statistics (e.g., masked_consensus_umis.tsv)
│       ├── raw/           # Raw statistics
│       └── target/        # Target-specific stats (e.g., cluster reports)
└── freebayes/           # Variant calling outputs (if using freebayes)
    ├── barcodeXX/
    │   ├── consensus/     # VCF files from consensus BAM (e.g., consensus.vcf)
    │   └── final/         # VCF files from final BAM (e.g., final.vcf)


> The above structure is derived from the output structure file. :contentReference[oaicite:8]{index=8}

---

## Conclusion

This guide splits the UMI pipeline outputs by subworkflow:

- **LIVE UMI PROCESSING** and **OFFLINE UMI PROCESSING** handle initial input processing, mapping, and clustering to produce `processed_umis`.  
- **UMI POLISHING** refines these clusters to create high-quality consensus sequences and associated BAM/FastQ files.  
- **VARIANT CALLING** uses the consensus data to call genetic variants, outputting VCF files.

Using this document, users can navigate the complex output directory and understand which results are produced by which subworkflow, helping in troubleshooting, quality control, and downstream analysis.

