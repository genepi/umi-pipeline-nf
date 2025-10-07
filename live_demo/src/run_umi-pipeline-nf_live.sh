#!/bin/bash

# Path to the input FASTQ file
INPUT_FASTQ="live_demo/data/input_fastq/input.fastq.gz"
# Base directory for barcode folders
NEXTFLOW_INPUT_DIR="live_demo/data/fastq_pass"
NEXTFLOW_OUTPUT_DIR="live_demo/output"
# Number of barcode directories
NUM_BARCODES=2

# Number of iterations and interval passed as arguments, with defaults
NUM_ITERATIONS=${1:-5}   # Default 5 iterations
INTERVAL=${2:-30}        # Default 60 seconds between iterations

create_outdir() {
    mkdir -p "${NEXTFLOW_INPUT_DIR}"
    for ((i=1; i<=NUM_BARCODES; i++)); do
        # Format number with leading zero if < 10
        BARCODE=$(printf "barcode%02d" "$i")
        mkdir -p "${NEXTFLOW_INPUT_DIR}/${BARCODE}"
    done
}

simulate_sequencing() {
    echo "Starting sequencing simulation for $NUM_ITERATIONS iterations, interval: $INTERVAL seconds"

    for ((iter=1; iter<=NUM_ITERATIONS; iter++)); do
        # Randomly choose a barcode directory
        BARCODE_DIR=$(printf "barcode%02d" $((RANDOM % NUM_BARCODES + 1)))

        # Output file path
        OUTPUT_FASTQ="${NEXTFLOW_INPUT_DIR}/${BARCODE_DIR}/input_$(date +%s).fastq.gz"

        # Copy the input FASTQ file to simulate sequencing output
        cp "$INPUT_FASTQ" "$OUTPUT_FASTQ"
        echo "[$(date)] Iteration $iter: Copied to ${OUTPUT_FASTQ}"

        # Wait for the specified interval before next iteration
        sleep "$INTERVAL"
    done

    echo "Sequencing simulation complete."
}

start_pipeline() {
    echo "Starting Nextflow pipeline..."
    nextflow run genepi/umi-pipeline-nf \
	-r Add_test_cases -latest \
	--live \
	--reference_based_polishing \
	--input "$NEXTFLOW_INPUT_DIR" \
	--output "$NEXTFLOW_OUTPUT_DIR" \
	-profile test,development
}

stop_pipeline() {
    echo "Stopping pipeline..."
    touch "$NEXTFLOW_OUTPUT_DIR/CONTINUE"
}

clean_up() {
	rm -r "${NEXTFLOW_INPUT_DIR}"
#	rm -r "${NEXTFLOW_OUTPUT_DIR}"
}

# Main workflow
create_outdir

# Start Nextflow pipeline in the background and save its PID
start_pipeline &
PIPELINE_PID=$!

# Run sequencing simulation in the foreground
simulate_sequencing

# Stop signal after sequencing ends
stop_pipeline

# Wait for Nextflow pipeline to complete
wait $PIPELINE_PID

# Now safe to clean up
clean_up
