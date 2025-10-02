# UMI Pipeline Subworkflow Output Overview

This document details the outputs produced by each subworkflow of the UMI analysis pipeline. The pipeline is split into distinct subworkflows, each handling a specific part of the processing—from raw input to variant calling. This guide explains the outputs, the processing steps/modules used, and how these outputs are organized.

![Workflow](images/umi-pipeline-nf_metro-map.svg)

## Pipeline Overview

The pipeline performs the following major steps:

- **Merge and Filter Reads:**  
  Input FastQ files are merged and filtered (using tools similar to *catfishq*) to remove low-quality reads.

- **Subsample Reads (Offline Only):**  
  Optionally, merged and filtered reads are subsampled using tools like *seqtk sample*.  
  *Output:* Subsampled data is saved with statistics in `<output>/<barcodeXX>/raw/subsampling/`.

- **Align Reads:**  
  Reads are aligned to a reference genome using [minimap2](https://github.com/lh3/minimap2).  
  *Output:* Alignment files are stored in `<output>/<barcodeXX>/raw/align`.

- **Extract UMI Sequences:**  
  UMI sequences are extracted from the aligned reads using a Python script, embedding the raw read information in the header.  
  *Output (only in verbose):* Extracted UMIs are saved in `<output>/<barcodeXX>/<target><fasta/fastq>_umi/<type>/`.

- **Cluster UMI Sequences:**  
  Reads are clustered by their UMI sequences (using vsearch) and then filtered based on a defined Hamming distance (default: 3). Clusters with too few reads are discarded, and overly abundant clusters may be downsampled.  
  *Output (only in verbose):* Cluster files are stored in `<output>/<barcodeXX>/<target>clustering/<type>/` and detailed cluster stats in `<output>/<barcodeXX>/<target>stats/<type>/`.

- **Polish Cluster:**  
  Clusters are polished (using medaka) to generate consensus sequences per cluster.  
  *Output (only in verbose):* Consensus sequences are produced in `<output>/<barcodeXX>/<target>clustering/<type>/clusters`.

- **Align Consensus Reads:**  
  Consensus sequences are realigned to the reference genome to verify quality.  
  *Output:* Aligned consensus files are stored in `<output>/<barcodeXX>/<target>align/<type>/`.

- **Extract & Cluster Consensus UMIs:**  
  UMIs are re-extracted from consensus reads and then re-clustered and parsed using a Python script.  
  *Output (only in verbose):* Final consensus clusters are parsed and saved in `<output>/<barcodeXX>/<target>align/<type>/`.

- **Call Variants:**  
  Variant calling is performed on the final consensus sequences using one of several supported callers (LOFREQ, MUTSERVE, or FREEBAYES).  
  *Output:* VCF files are produced in `<output>/<freebayes/mutserve/lowfreq>/<barcodeXX>/<target>/<type>/`.

- **Pipeline Info:**  
  Nextflow’s built-in reporting generates files such as a DAG visualization (`dag.svg`), a detailed report (`report.html`), a timeline (`timeline.html`), and a trace file (`trace.txt`).  
  *Output:* These files are stored in `<output>/nextflow_stats/`.

## Table of Contents

1. [LIVE UMI PROCESSING](#live-umi-processing)
2. [OFFLINE UMI PROCESSING](#offline-umi-processing)
3. [UMI POLISHING](#umi-polishing)
4. [VARIANT CALLING](#variant-calling)
5. [Overall Directory Structure](#overall-directory-structure)
6. [Conclusion](#conclusion)

---

## LIVE UMI PROCESSING

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


---

## OFFLINE UMI PROCESSING

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


---

## UMI POLISHING

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


---

## VARIANT CALLING

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


---

## Overall Directory Structure

Below is an overview of the directory structure created after a complete pipeline run. This structure illustrates where outputs from each subworkflow are stored.

``` plaintext
output_directory/  
├── barcodeXX/                  # results for BarcodeXX
│   ├── raw/                    # target unspecific outputs (only in verbose)
│   │   ├── align/              # Raw alignments (e.g., filtered_barcodeXX.1.bam) (only in verbose)    
│   │   └── subsampling/        # subsampled reads
│   ├── targetX/                # results for targetX
│   │   ├── align
│   │   │   ├── consensus/      # Consensus alignments (e.g., masked_consensus.bam)  
│   │   │   └── final/          # Final alignments (e.g., final.bam )  
│   │   ├── clustering/  
│   │   │   ├── raw/            # Raw clustering outputs (only in verbose)  
│   │   │   └── consensus/      # Consensus clustering outputs (e.g., consensus.fasta) (only in verbose)  
│   │   ├── fastq/  
│   │   │   ├── consensus/      # Consensus FastQ files (e.g., masked_consensus.fastq)  
│   │   │   └── final/          # Final FastQ files (e.g., final.fastq)  
│   │   ├── fastq_filtered/  
│   │   │   └── raw/            # Filtered FastQ files (only in verbose) 
│   │   ├── fastq_umi/  
│   │   │   └── consensus/      # UMI-specific consensus FastQ files (only in verbose)  
│   │   ├── polishing/  
│   │   │   └── consensus/      # Polishing outputs and logs (only in verbose)  
│   │   └── stats/  
│   │       ├── consensus/      # Consensus statistics (e.g., masked_consensus_umis.tsv)  
│   │       └── raw/            # Raw statistics  
│   └── targetY/                # results for targetY
└── freebayes/                  # Variant calling outputs (if using freebayes)  
    └── barcodeXX/  
        ├── targetX/            # results for targetX 
        │   ├── consensus/      # VCF files from consensus BAM (e.g., consensus.vcf)  
        │   └── final/          # VCF files from final BAM (e.g., final.vcf)  
        └── targetY/

```
---

## Conclusion

This guide splits the UMI pipeline outputs by subworkflow:

- **LIVE UMI PROCESSING** and **OFFLINE UMI PROCESSING** handle initial input processing, mapping, and clustering to produce `processed_umis`.  
- **UMI POLISHING** refines these clusters to create high-quality consensus sequences and associated BAM/FastQ files.  
- **VARIANT CALLING** uses the consensus data to call genetic variants, outputting VCF files.

Using this document, users can navigate the complex output directory and understand which results are produced by which subworkflow, helping in troubleshooting, quality control, and downstream analysis.

