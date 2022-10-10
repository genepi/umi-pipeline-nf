process MAP_1D {
    
    input:
        path fastq
        file reference 
    output:
        file "align/1d.bam", emit: bam_1d
        file "align/ad.bam.bai", emit: bai_1d

    """
        catfishq --max_n ${params.subset_reads} ${fastq} | minimap2 {params.minimap2_param} -t {threads} {input.REF} - | samtools sort -@ 5 -o {output.BAM} - && samtools index -@ {threads} {output.BAM}
}