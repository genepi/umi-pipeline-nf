#!/bin/bash

# ----------------------------------------
# Predefined settings
# ----------------------------------------
INPUT_FASTQ="live_demo/data/input_fastq/input.fastq.gz"
NEXTFLOW_INPUT_DIR="live_demo/data/fastq_pass"
NEXTFLOW_OUTPUT_DIR="live_demo/output"
NUM_BARCODES=2

# ----------------------------------------
# Default parameters
# ----------------------------------------
NUM_ITERATIONS=5
INTERVAL=30

# ----------------------------------------
# Parse arguments
# ----------------------------------------
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --iterations|-i) NUM_ITERATIONS="$2"; shift ;;
        --interval|-t) INTERVAL="$2"; shift ;;
        --help|-h)
            echo "Usage: $0 [--iterations <num>] [--interval <seconds>]"
            echo "  --iterations, -i    Number of sequencing iterations (default: 5)"
            echo "  --interval, -t      Time interval between iterations in seconds (default: 30)"
            exit 0
            ;;
        *)
            echo "Unknown parameter passed: $1"
            echo "Use --help for usage information."
            exit 1
            ;;
    esac
    shift
done

# ----------------------------------------
# Functions
# ----------------------------------------
create_outdir() {
    mkdir -p "${NEXTFLOW_INPUT_DIR}"
    for ((i=1; i<=NUM_BARCODES; i++)); do
        BARCODE=$(printf "barcode%02d" "$i")
        mkdir -p "${NEXTFLOW_INPUT_DIR}/${BARCODE}"
    done
}

simulate_sequencing() {
    echo "Starting sequencing simulation for $NUM_ITERATIONS iterations, interval: $INTERVAL seconds"

    for ((iter=1; iter<=NUM_ITERATIONS; iter++)); do
        BARCODE_DIR=$(printf "barcode%02d" $((RANDOM % NUM_BARCODES + 1)))
        OUTPUT_FASTQ="${NEXTFLOW_INPUT_DIR}/${BARCODE_DIR}/input_$(date +%s).fastq.gz"
        cp "$INPUT_FASTQ" "$OUTPUT_FASTQ"
        echo "[$(date)] Iteration $iter: Copied to ${OUTPUT_FASTQ}"
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
    echo "Cleaning up input directories..."
    rm -rf "${NEXTFLOW_INPUT_DIR}"
    # Uncomment the next line to remove output as well
    # rm -rf "${NEXTFLOW_OUTPUT_DIR}"
}

# ----------------------------------------
# Main workflow
# ----------------------------------------
create_outdir

start_pipeline &
PIPELINE_PID=$!

simulate_sequencing
stop_pipeline

wait $PIPELINE_PID
clean_up
