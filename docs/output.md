# umi-pipeline-nf Output
This document describes the output produced by the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

[![Workflow_flow_chart](https://mermaid.ink/img/pako:eNqlWGFvmzAQ_SvI-9JKWUXSZM2itRJtyRYpabqQTNrGhFwwDRoBhk23rLS_fQesIRgTJxv5gGP7vTvf2c-5PCI7dAgaINcPf9pLHDNlfmkGZqDAQ5O7-xhHSyWhJLa8IEqY8tVEiom-FROypxyz2JGJ3rlhwBTq_SbnHfVCWcCoMspGTXRcgihb-6QCVVzP9weUxeF3MnilqmqraL_-6TlsOVCjX2KTZW_R4WLKflgRpvToqOrOmwulHLQcL4aVHB8f8wR3xBEgobcZEhOXxCSwiQC4GePhJHCUosUFGwfYX1NiUbyKIEiVgAMHsZkXBnmWXryoIvg0tCENWjFDMfIZglTwFIekg8c-Pz9v5WJ7jVyaGjMUE-xQPtZlK4pDm0AO7TBaZ_n6yrNUv11Nbz9bl_p1Zdu-cKxIfE8Kd2o0E332XreGmjH_KMRCyvI1e8F9DWssLg1tcjse3bwX28VRscy6Ve3WmunatSG2CeZYA9IAc_MdWIcw2D1WsvLyBeMawbU-16_m1mIyyletCVlsP6GMxDXw1XhhzPWZEAOnIIxXGPLu-YBt5Jjpw-lsos2t4WgMXNYuzij0PbpspLqdjkfGh50MWQ7sMKAkoIk4D1fTG0O_MRbSeG5o5JHdcMpjvMO7v-uSeLiJe-mfNPSlf7LYuR6cfEkEh6MbbSzx8gHHHg7ASexnZ8laJQwU_oHUCRdzQ5990vdicWNC7vCa1P0aznT9Uvus7-eNHwJTXRzGU6Cp6kLZgpg4XqbT2xrxyFNsjZnoSYTmfKkxcOMVFsHFppyfX6TwTmsS2iyudcyWBDVbyxU1xxaSDgef0FQkuzs1OWMQR1MS7BwoEOqdKp77m1t3CvNpXa9l7v4bxctQ1W1RoAWDFRAv8vJroAL_Kw-NklSZ3CDre-t_hayq5zK9r0ArQn4QslnCD5F7UQR3OVSbIo5qTbQPEfhagDi13kvSuePHy41csXKC9EXQU5ni7824EfdUKv97cxZCn-6-B5oVb_NLP1eR8nd_vkNS2V49mKMxnQcxyRPxX3SCJPwXH5-Dg6O2JaZQgBWN7XISahe-nOEn1aYUE_jyM99TtTqIuxwVnin7oBaCC2SFPQfq8scMZyK2JCuo2wbQdIiLEx_qaTN4gqk4YaGxDmw0cLFPSQslkYMZufYwlJOrTW-EAzR4RL_QoN3unPRUtdfrd952T-HVbqE1dJ_0e321q7a73bPeab_de2qh32EIDOrJ267abZ_1umr_tHP2ptNrIQIHKYwnxX8H-V8IuYkvOYDFCXn6A0LNnOg?type=png)](https://mermaid-js.github.io/mermaid-live-editor/edit#pako:eNqlWGFvmzAQ_SvI-9JKWUXSZM2itRJtyRYpabqQTNrGhFwwDRoBhk23rLS_fQesIRgTJxv5gGP7vTvf2c-5PCI7dAgaINcPf9pLHDNlfmkGZqDAQ5O7-xhHSyWhJLa8IEqY8tVEiom-FROypxyz2JGJ3rlhwBTq_SbnHfVCWcCoMspGTXRcgihb-6QCVVzP9weUxeF3MnilqmqraL_-6TlsOVCjX2KTZW_R4WLKflgRpvToqOrOmwulHLQcL4aVHB8f8wR3xBEgobcZEhOXxCSwiQC4GePhJHCUosUFGwfYX1NiUbyKIEiVgAMHsZkXBnmWXryoIvg0tCENWjFDMfIZglTwFIekg8c-Pz9v5WJ7jVyaGjMUE-xQPtZlK4pDm0AO7TBaZ_n6yrNUv11Nbz9bl_p1Zdu-cKxIfE8Kd2o0E332XreGmjH_KMRCyvI1e8F9DWssLg1tcjse3bwX28VRscy6Ve3WmunatSG2CeZYA9IAc_MdWIcw2D1WsvLyBeMawbU-16_m1mIyyletCVlsP6GMxDXw1XhhzPWZEAOnIIxXGPLu-YBt5Jjpw-lsos2t4WgMXNYuzij0PbpspLqdjkfGh50MWQ7sMKAkoIk4D1fTG0O_MRbSeG5o5JHdcMpjvMO7v-uSeLiJe-mfNPSlf7LYuR6cfEkEh6MbbSzx8gHHHg7ASexnZ8laJQwU_oHUCRdzQ5990vdicWNC7vCa1P0aznT9Uvus7-eNHwJTXRzGU6Cp6kLZgpg4XqbT2xrxyFNsjZnoSYTmfKkxcOMVFsHFppyfX6TwTmsS2iyudcyWBDVbyxU1xxaSDgef0FQkuzs1OWMQR1MS7BwoEOqdKp77m1t3CvNpXa9l7v4bxctQ1W1RoAWDFRAv8vJroAL_Kw-NklSZ3CDre-t_hayq5zK9r0ArQn4QslnCD5F7UQR3OVSbIo5qTbQPEfhagDi13kvSuePHy41csXKC9EXQU5ni7824EfdUKv97cxZCn-6-B5oVb_NLP1eR8nd_vkNS2V49mKMxnQcxyRPxX3SCJPwXH5-Dg6O2JaZQgBWN7XISahe-nOEn1aYUE_jyM99TtTqIuxwVnin7oBaCC2SFPQfq8scMZyK2JCuo2wbQdIiLEx_qaTN4gqk4YaGxDmw0cLFPSQslkYMZufYwlJOrTW-EAzR4RL_QoN3unPRUtdfrd952T-HVbqE1dJ_0e321q7a73bPeab_de2qh32EIDOrJ267abZ_1umr_tHP2ptNrIQIHKYwnxX8H-V8IuYkvOYDFCXn6A0LNnOg)

* [Merge and Filter Reads](#merge-and-filter-reads) - Merge and filter input fastq files
* [Subsample Reads](#subsampling) - Subsample merged and filtered reads
* [Align Reads](#align-reads) - Align all reads to the reference genome
* [Separate Amplicons](#separate-amplicons) - Separate all reads into amplicons
* [Extract UMI](#extract-umi-sequences) - Extract UMI sequences of all reads
* [Cluster UMI](#cluster-umi-sequences) - Cluster UMI sequences per amplicon  
* [Extract Consensus Reads](#extract-conensus-sequences) - Extract consensus reads of UMI clusters  
* [Align Consensus Reads](#align-consensus-reads) - Align all consensus reads to the reference genome 
* [Call Variants](#call-variants) - Perform a variant calling of consensus reads (optional)
* [Haplotyping](#haplotyping) - Perform haplotyping of consensus reads (optional)
* [Pipeline Info](#pipeline-info) - reports from nextflow about the pipeline run

## Merge and Filter Reads
Input fastq files will be merged and filtered - if specified by user - using catfishq.

## Subsampling
If specified by user the merged and filtered reads will be subsampled using seqtk sample.
**Output directory: `<output>/<barcodeXX>/subsampling/`**

Information regarding subsampling will be saved in a tsv file 
**Output directory: `<output>/<barcodeXX>/stats/`**

## Align Reads
Merged - filtered and subsampled - reads will be aligned to the provided reference genome using minimap2.
**Output directory: `<output>/<barcodeXX>/align/raw/`**

## Separate Amplicons

## Extract UMI Sequences

## Cluster UMI Sequences

## Extract Consensus Sequences

## Align Consensus Reads

## Call Variants

## Haplotyping

## Pipeline Info
Nextflow has several built-in reporting tools that give information about the pipeline run.

**Output directory: `<output>/nextflow_stats/`**

* `dag.svg`
  * MMD file, which can be visualized using [Mermaid](https://mermaid-js.github.io/mermaid/#/)
* `report.html`
  * Nextflow report describing parameters, computational resource usage and task bash commands used.
* `timeline.html`
  * A waterfall timeline plot showing the running times of the workflow tasks.
* `trace.txt`
  * A text file with machine-readable statistics about every task executed in the pipeline.