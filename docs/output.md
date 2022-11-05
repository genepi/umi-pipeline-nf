# umi-pipeline-nf Output
This document describes the output produced by the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

[![Workflow](https://mermaid.ink/img/pako:eNqlWPFv2jgU_lci3y-txFWBwcbQtVLahh0SlB6BSdtyitzEKdGFJBcn21jT_u334qyE2A6GW_ghxvb3vef37M88npAbewSNkB_G39w1TjNteW1HdqTBQ_OHxxQnay2nJHWCKMkz7YuNNBv9XU0on3rMyc5s9IcfR5lGgx_ksqdfaSsY1SblqI3OaxDNtiFpQDU_CMMRzdL4HzL6Tdf1TtX-_VvgZeuRnnyXm6x7qw4f0-xfJ8GUnp013Xl7pdWDjheksJLz83Oe4IF4EiT0tkNS4pOURC6RAHdjPJxEnla1uGDjCIdbShyKNwkEqRFw4CBuFsQRy9KrF00En4YupMGoZmgWmyFJBU9xSjp47MvLy14u9tfIpak1QynBHuVjXbeSNHYJ5NCNk22Zry88S_Pbzfz-k3Nt3ja27SvHhqSPpHJHoJmZiw-mMzas5V9SLKSMrTmIHgWstbq2jNn9dHL3QW4XJ9UyRavGvbMwjVtLbhPMZS1IC8wtD2A9ksHucfJNwBaMBYJbc2neLJ3VbMJWbUhZ3DCnGUkF8M10ZS3NhRQDpyBONxjyHoSAbeVYmOP5YmYsnfFkClzOIc4kDgO6bqW6n08n1p8HGcocuHFESURzeR5u5neWeWetlPHc0agju-NUx_iAdz_XpfBwF_faP2Xoa_9UsfMDOPmKCI4nd8ZU4eVXnAY4AidxWJ4lZ5NnoPBfiUi4Wlrm4qN5FIufEvKAt0T0a7wwzWvjk3mcN2EMTKI4TOdA09SFugUx8YJSp_c14omn2Buz0bMMzfkiMHDjDRbJxaZdXl4V8C4ECW0XVxGzJ0Ht1piiMmwl6XDwCS1ksntQk0sGeTQVwWZAiVAfVHHmL7PuVeYLUa9V7v4_itehptuyQEsGGyBe5NXXQAP-Ux5aJakxuUXWj9b_BllTz1V634A2hPwkZLuEnyL3sggeckiYIo-qINqnCLwQIE6tj5J07vjxcqNWLEZQvAp6oVL8oxl34l4o5f9ozkroi8P3QLvi7X7pMxWpf_ezHVKo9urJHK3pPIlJnYhfopMk4Zf4-BycHLU9MYUCrGrsl5NQu_DlDD9JmFJN4MtPtqeEOoi7HDWeqfygDoILZIMDD-rypxJno2xNNlC3jaDpER_nIdTTdvQMU3GexdY2ctHIxyElHZQnHs7IbYChnNzsehMcodET-o5G3W7vYqDrg8Gw977_Bl7dDtpC98VwMNT7erfffzd4M-wOnjvoRxwDg37xvq_3u-8GPejt9d8Ohh1E4CDF6az674D9hcBMfGaALM3J839DY5zs?type=png)](https://mermaid-js.github.io/mermaid-live-editor/edit#pako:eNqlWPFv2jgU_lci3y-txFWBwcbQtVLahh0SlB6BSdtyitzEKdGFJBcn21jT_u334qyE2A6GW_ghxvb3vef37M88npAbewSNkB_G39w1TjNteW1HdqTBQ_OHxxQnay2nJHWCKMkz7YuNNBv9XU0on3rMyc5s9IcfR5lGgx_ksqdfaSsY1SblqI3OaxDNtiFpQDU_CMMRzdL4HzL6Tdf1TtX-_VvgZeuRnnyXm6x7qw4f0-xfJ8GUnp013Xl7pdWDjheksJLz83Oe4IF4EiT0tkNS4pOURC6RAHdjPJxEnla1uGDjCIdbShyKNwkEqRFw4CBuFsQRy9KrF00En4YupMGoZmgWmyFJBU9xSjp47MvLy14u9tfIpak1QynBHuVjXbeSNHYJ5NCNk22Zry88S_Pbzfz-k3Nt3ja27SvHhqSPpHJHoJmZiw-mMzas5V9SLKSMrTmIHgWstbq2jNn9dHL3QW4XJ9UyRavGvbMwjVtLbhPMZS1IC8wtD2A9ksHucfJNwBaMBYJbc2neLJ3VbMJWbUhZ3DCnGUkF8M10ZS3NhRQDpyBONxjyHoSAbeVYmOP5YmYsnfFkClzOIc4kDgO6bqW6n08n1p8HGcocuHFESURzeR5u5neWeWetlPHc0agju-NUx_iAdz_XpfBwF_faP2Xoa_9UsfMDOPmKCI4nd8ZU4eVXnAY4AidxWJ4lZ5NnoPBfiUi4Wlrm4qN5FIufEvKAt0T0a7wwzWvjk3mcN2EMTKI4TOdA09SFugUx8YJSp_c14omn2Buz0bMMzfkiMHDjDRbJxaZdXl4V8C4ECW0XVxGzJ0Ht1piiMmwl6XDwCS1ksntQk0sGeTQVwWZAiVAfVHHmL7PuVeYLUa9V7v4_itehptuyQEsGGyBe5NXXQAP-Ux5aJakxuUXWj9b_BllTz1V634A2hPwkZLuEnyL3sggeckiYIo-qINqnCLwQIE6tj5J07vjxcqNWLEZQvAp6oVL8oxl34l4o5f9ozkroi8P3QLvi7X7pMxWpf_ezHVKo9urJHK3pPIlJnYhfopMk4Zf4-BycHLU9MYUCrGrsl5NQu_DlDD9JmFJN4MtPtqeEOoi7HDWeqfygDoILZIMDD-rypxJno2xNNlC3jaDpER_nIdTTdvQMU3GexdY2ctHIxyElHZQnHs7IbYChnNzsehMcodET-o5G3W7vYqDrg8Gw977_Bl7dDtpC98VwMNT7erfffzd4M-wOnjvoRxwDg37xvq_3u-8GPejt9d8Ohh1E4CDF6az674D9hcBMfGaALM3J839DY5zs)


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
