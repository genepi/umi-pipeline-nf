#!/usr/bin/env python3
import pysam
import argparse
import sys
import os

def sort_and_index_bam(input_bam: str, n_threads: int = 1):
    """
    Sorts a BAM by coordinate and then creates its index (.bai).
    Replaces the input BAM with its sorted version.
    """
    base, ext = os.path.splitext(input_bam)
    sorted_bam = "{base}.sorted{ext}".format(base=base, ext=ext)

    # Sort BAM by coordinate
    pysam.sort(
        "-@", str(n_threads),
        "-o", sorted_bam,
        input_bam
    )
    os.replace(sorted_bam, input_bam)
    print("→ BAM sorted in-place: {bam}".format(bam=input_bam))

    # Index the sorted BAM
    pysam.index(input_bam)
    print("→ Index (.bai) created: {bam}.bai".format(bam=input_bam))


def parse_alignment(bam_file, reference_fasta, output_reference, output_bam, n_threads=1):
    # 1) Collect all cluster IDs from the read names
    cluster_ids = set()
    with pysam.AlignmentFile(bam_file, "rb") as in_bam:
        for read in in_bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            cid = "_".join(read.query_name.split('_')[:2])
            cluster_ids.add(cid)
    cluster_ids = sorted(cluster_ids)

    # 2) Build a brand-new header dict
    with pysam.AlignmentFile(bam_file, "rb") as in_bam:
        orig = in_bam.header.to_dict()
    new_header = {}

    if 'HD' in orig:
        new_header['HD'] = orig['HD']
    if 'PG' in orig:
        new_header['PG'] = orig['PG']

    new_header['SQ'] = [{'SN': cid, 'LN': 5104} for cid in cluster_ids] 

    # 3) Write unsorted output BAM with this header
    temp_bam = "{out}.unsorted".format(out=output_bam)
    out = pysam.AlignmentFile(temp_bam, "wb", header=new_header)

    # 4) Rewrite every mapped read using new header
    with pysam.AlignmentFile(bam_file, "rb") as in_bam:
        for read in in_bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            cid = "_".join(read.query_name.split('_')[:2])

            new_read = pysam.AlignedSegment(out.header)
            new_read.query_name = read.query_name
            new_read.query_sequence = read.query_sequence
            new_read.flag = read.flag
            new_read.reference_start = read.reference_start
            new_read.mapping_quality = read.mapping_quality
            new_read.cigar = read.cigar
            new_read.next_reference_id = read.next_reference_id
            new_read.next_reference_start = read.next_reference_start
            new_read.template_length = read.template_length
            new_read.query_qualities = read.query_qualities
            new_read.tags = read.tags

            new_read.reference_name = cid
            if new_read.next_reference_id != -1:
                new_read.next_reference_name = cid

            out.write(new_read)
    out.close()

    # 5) Sort and index the BAM in-place
    sort_and_index_bam(temp_bam, n_threads=n_threads)
    os.replace(temp_bam, output_bam)
    os.replace("{bam}.bai".format(bam=temp_bam), "{bam}.bai".format(bam=output_bam))

    # 6) Produce the multi-FASTA
    seq_parts = []
    with open(reference_fasta) as rf:
        for line in rf:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip())
    seq = "".join(seq_parts)

    with open(output_reference, "w") as fo:
        for cid in cluster_ids:
            fo.write(">{cid}\n{seq}\n".format(cid=cid, seq=seq))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Cluster-aware BAM reheader + multi-SEQUENCE FASTA generator"
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM")
    parser.add_argument("-r", "--reference", required=True, help="One-seq FASTA")
    parser.add_argument("-or", "--output_reference", required=True, help="Clustered FASTA")
    parser.add_argument("-ob", "--output_bam", required=True, help="Output BAM")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Threads for sorting")
    args = parser.parse_args()

    parse_alignment(
        bam_file=args.bam,
        reference_fasta=args.reference,
        output_reference=args.output_reference,
        output_bam=args.output_bam,
        n_threads=args.threads
    )
