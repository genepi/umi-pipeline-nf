import pysam
import argparse

def parse_alignment(bam_file, reference_fasta, output_reference, output_bam):
    # First pass: collect unique cluster IDs from the input BAM.
    cluster_ids = set()
    with pysam.AlignmentFile(bam_file, "rb") as in_bam:
        for read in in_bam:
            if read.is_unmapped:
                continue
            # Extract cluster id from the read name (adjust your split as needed)
            cluster_id = "_".join(read.query_name.split('_')[:2])
            cluster_ids.add(cluster_id)
    
    # Build a new header containing only the cluster references.
    # Here we set a dummy length (5104) for each new reference.
    new_sq = [{'SN': cluster, 'LN': 5104} for cluster in sorted(cluster_ids)]
    new_header = {'HD': {'VN': '1.0'}, 'SQ': new_sq}
    
    # Open output BAM with the new header.
    out_bam = pysam.AlignmentFile(output_bam, "wb", header=new_header)
    
    # Process the reads and update their reference id using the new header.
    with pysam.AlignmentFile(bam_file, "rb") as in_bam:
        for read in in_bam:
            if read.is_unmapped:
                continue
            cluster_id = "_".join(read.query_name.split('_')[:2])
            # Look up the reference id using the new header
            new_ref_id = out_bam.get_tid(cluster_id)
            if new_ref_id < 0:
                raise ValueError(f"Reference {cluster_id} not found in new header")
            read.reference_id = new_ref_id
            out_bam.write(read)
    out_bam.close()
    
    # Process the FASTA file to build the reference sequence string.
    with open(reference_fasta, "r") as ref_file:
        ref_seq = "".join(line.strip() for line in ref_file if not line.startswith(">"))
    
    # Write a new FASTA with one entry per cluster.
    with open(output_reference, "w") as fasta_out:
        for cluster in sorted(cluster_ids):
            fasta_out.write(f">{cluster}\n{ref_seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate cluster-specific reference FASTA and update BAM reference names."
    )
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-r", "--reference", required=True, help="Single reference FASTA file")
    parser.add_argument("-or", "--output_reference", required=True, help="Output clustered reference FASTA file")
    parser.add_argument("-ob", "--output_bam", required=True, help="Output BAM file with updated references")
    
    args = parser.parse_args()
    parse_alignment(args.bam, args.reference, args.output_reference, args.output_bam)
