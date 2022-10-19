# umi-pipeline-nf Output
This document describes the output produced by the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

[![Workflow Flow Chart](https://mermaid.ink/img/pako:eNqFVl2PmkAU_SuGJzfZNjIzIvrQhFXcmvhBBR-2tZlMYVxNBSxg283u_veOChfuuJvyxJx75tx75g43PBthGkljYGz26Z9wK7KiFYzWSUs9h067PdyKJJH7j5ssjT1RbG9uypjZbrdgQb6tZhPuTTx3Opm7g-HCe-B37uh7GaY4PHOX9y4fO37wpWIwzPBXd74z89T7fcXoNvNZzUWvubC1VI7Hl64z8iuZfpNsdtDK1IpQ-QO82yRoA0UrhlZdLDZyA3cY8BN2Mu6AIrJiIi-mZmY4XfmBu4StyAtBXghuj9afpTteLGdOwMeTqRLkmjBBvgjyRTRf3mI68T9fKSBXBLkiNlr1rxs2XMx9d-6v4Ngp8kaRN4paQum7xw6quAEU2aPomlHrzQa8USBySJFDitrEkBVmvtOYuljtZBmyy1CnGLs-y_Fk7kyvC2bIJ0PtYj0sM1667p3z4Nabkb9O68OHTy8bkRe_-Ga3lzkPty_KdnVLz-EfMlIYqZrewEyzOqcTqMSrJGdOLLNHGfGzuiLbzaCaCdU4OHOLp4NscHpnMJMbmckkbEbsy27I3G9uB9S8GDs8Fds0OTkrZMYzKaIcscxSjOhAF6YEkgf4YvAY77j8W2QiLPglFSJ1SzFbByyYGEgdeKVHQv4HmH0YIE0l4BET6jyILJc83B9zdRR5XW7NJWWSrg5QGCYoCfBKW8TSAagO-axh-6rLdax_0aAdTZRSmCRN0Rom77emJpX3lVo6wGCWIHXgWeVFNzWA2jrQgxHSVIKNrAN1qgNIs1gUPEyTXCb5MUfE8k4y-LgIEgSYXp1mHSu_OdaDEYI0ALauNSB2aeLvcHOeFKeQbdwa6iOPxS5S_yDPJ-LaKLYylmtjoF4jkf1cG-vkVfHEsUj9pyQ0Bhuxz-WtcTxEopCjnXjMRAyojHZFms0ufzXnn5tb4yCSr2mqOEV2lK__AMtMa8c?type=png)](https://mermaid-js.github.io/mermaid-live-editor/edit#pako:eNqFVl2PmkAU_SuGJzfZNjIzIvrQhFXcmvhBBR-2tZlMYVxNBSxg283u_veOChfuuJvyxJx75tx75g43PBthGkljYGz26Z9wK7KiFYzWSUs9h067PdyKJJH7j5ssjT1RbG9uypjZbrdgQb6tZhPuTTx3Opm7g-HCe-B37uh7GaY4PHOX9y4fO37wpWIwzPBXd74z89T7fcXoNvNZzUWvubC1VI7Hl64z8iuZfpNsdtDK1IpQ-QO82yRoA0UrhlZdLDZyA3cY8BN2Mu6AIrJiIi-mZmY4XfmBu4StyAtBXghuj9afpTteLGdOwMeTqRLkmjBBvgjyRTRf3mI68T9fKSBXBLkiNlr1rxs2XMx9d-6v4Ngp8kaRN4paQum7xw6quAEU2aPomlHrzQa8USBySJFDitrEkBVmvtOYuljtZBmyy1CnGLs-y_Fk7kyvC2bIJ0PtYj0sM1667p3z4Nabkb9O68OHTy8bkRe_-Ga3lzkPty_KdnVLz-EfMlIYqZrewEyzOqcTqMSrJGdOLLNHGfGzuiLbzaCaCdU4OHOLp4NscHpnMJMbmckkbEbsy27I3G9uB9S8GDs8Fds0OTkrZMYzKaIcscxSjOhAF6YEkgf4YvAY77j8W2QiLPglFSJ1SzFbByyYGEgdeKVHQv4HmH0YIE0l4BET6jyILJc83B9zdRR5XW7NJWWSrg5QGCYoCfBKW8TSAagO-axh-6rLdax_0aAdTZRSmCRN0Rom77emJpX3lVo6wGCWIHXgWeVFNzWA2jrQgxHSVIKNrAN1qgNIs1gUPEyTXCb5MUfE8k4y-LgIEgSYXp1mHSu_OdaDEYI0ALauNSB2aeLvcHOeFKeQbdwa6iOPxS5S_yDPJ-LaKLYylmtjoF4jkf1cG-vkVfHEsUj9pyQ0Bhuxz-WtcTxEopCjnXjMRAyojHZFms0ufzXnn5tb4yCSr2mqOEV2lK__AMtMa8c)

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