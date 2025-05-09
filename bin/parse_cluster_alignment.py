import pysam
import argparse
import os

def sort_and_index_bam(input_bam, sorted_bam, n_threads=1):
    pysam.sort("-@", str(n_threads), "-o", sorted_bam, input_bam)
    pysam.index(sorted_bam)

def parse_fasta(reference_fasta):
    """
    Parses a multi-FASTA file into a dict {seq_name: seq_string}.
    """
    sequences = {}
    current_name = None
    current_seq = []

    with open(reference_fasta, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name is not None:
            sequences[current_name] = "".join(current_seq)

    return sequences

def parse_alignment(bam_file, reference_fasta, output_reference, output_bam):
    # Load all reference sequences
    reference_seqs = parse_fasta(reference_fasta)

    # First pass: collect cluster IDs and remember which original ref each maps to
    cluster_ids = set()
    cluster_to_ref = {}
    ref_lengths = {}

    with pysam.AlignmentFile(bam_file, "rb") as in_bam:
        # record original reference lengths
        for sq in in_bam.header["SQ"]:
            ref_lengths[sq["SN"]] = sq["LN"]

        for read in in_bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue
            cluster_id = "_".join(read.query_name.split('_')[:2])
            cluster_ids.add(cluster_id)
            # only set once per cluster
            if cluster_id not in cluster_to_ref:
                cluster_to_ref[cluster_id] = read.reference_name

    # Build new header SQ entries, and a mapping cluster -> new_name
    new_sq = []
    cluster_to_newname = {}
    for cluster in sorted(cluster_ids):
        oldref = cluster_to_ref[cluster]
        length = ref_lengths.get(oldref)
        if length is None:
            raise ValueError("Reference '{}' not found in BAM header".format(oldref))

        new_name = "{}_{}".format(oldref, cluster)
        cluster_to_newname[cluster] = new_name
        new_sq.append({"SN": new_name, "LN": length})

    new_header = {"HD": {"VN": "1.0"}, "SQ": new_sq}

    # Write an unsorted BAM with the new header
    base_out = output_bam.rsplit(".", 1)[0]
    temp_unsorted = "{}.unsorted.bam".format(base_out)
    out_bam = pysam.AlignmentFile(temp_unsorted, "wb", header=new_header)

    with pysam.AlignmentFile(bam_file, "rb") as in_bam:
        for read in in_bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            cluster_id = "_".join(read.query_name.split('_')[:2])
            new_name = cluster_to_newname[cluster_id]

            # build a fresh AlignedSegment tied to the new header
            new_read = pysam.AlignedSegment(out_bam.header)
            new_read.query_name        = read.query_name
            new_read.query_sequence    = read.query_sequence
            new_read.flag              = read.flag
            new_read.reference_start   = read.reference_start
            new_read.mapping_quality   = read.mapping_quality
            new_read.cigar             = read.cigar
            new_read.next_reference_id = read.next_reference_id
            new_read.next_reference_start = read.next_reference_start
            new_read.template_length   = read.template_length
            new_read.query_qualities   = read.query_qualities
            new_read.tags              = read.tags

            # assign the new combined reference name
            new_read.reference_name    = new_name

            out_bam.write(new_read)

    out_bam.close()

    # Sort & index, then clean up
    sort_and_index_bam(temp_unsorted, output_bam)
    os.remove(temp_unsorted)

    # Write out the new multi-sequence FASTA
    with open(output_reference, "w") as fasta_out:
        for cluster in sorted(cluster_ids):
            oldref   = cluster_to_ref[cluster]
            new_name = cluster_to_newname[cluster]
            seq = reference_seqs.get(oldref)
            if seq is None:
                raise ValueError("Reference '{}' not found in input FASTA".format(oldref))
            fasta_out.write(">{}\n{}\n".format(new_name, seq))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate cluster-specific reference FASTA and update BAM references to '<oldref>_<cluster>'."
    )
    parser.add_argument("-b", "--bam",               required=True, help="Input BAM file")
    parser.add_argument("-r", "--reference",         required=True, help="Multi-sequence reference FASTA")
    parser.add_argument("-or", "--output_reference", required=True, help="Output clustered reference FASTA")
    parser.add_argument("-ob", "--output_bam",       required=True, help="Output BAM with updated references")
    args = parser.parse_args()

    parse_alignment(
        args.bam,
        args.reference,
        args.output_reference,
        args.output_bam
    )
