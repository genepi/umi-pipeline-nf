process PARSE_BED {
    publishDir "${params.output}/bed", mode: 'copy', enabled: "${params.verbose}"

    input:
    path bed

    output:
    path "*.bed", emit: bed_files

    script:
    """
    while read -r line; do
        # Split BED line by tab
        IFS=\$'\\t' read -r -a fields <<< "\$line"
        target=\$(echo "\${fields[3]}" | tr -d '[:space:]')

        # Write to a per-target BED file
        echo "\$line" >> \${target}.bed
    done < ${bed}
    """
}
