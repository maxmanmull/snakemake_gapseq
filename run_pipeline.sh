#!/bin/bash

# Fixed Pipeline Runner Script for non-SLURM systems

# Set defaults
CORES=32
JOBS=10
PROFILE=""
DRYRUN=false
FORCERUN=""
UNLOCK=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        -j|--jobs)
            JOBS="$2"
            shift 2
            ;;
        -n|--dryrun)
            DRYRUN=true
            shift
            ;;
        -f|--force)
            FORCERUN="$2"
            shift 2
            ;;
        --unlock)
            UNLOCK=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  -c, --cores NUM      Number of cores to use (default: 32)"
            echo "  -j, --jobs NUM       Number of parallel jobs (default: 10)"
            echo "  -n, --dryrun         Perform a dry run"
            echo "  -f, --force RULE     Force execution of specific rule"
            echo "  --unlock             Unlock the working directory"
            echo "  -h, --help           Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create necessary directories
mkdir -p logs
mkdir -p results

# Handle special cases
if [ "$UNLOCK" = true ]; then
    echo "Unlocking working directory..."
    snakemake --unlock
    exit 0
fi

# Build the basic Snakemake command (no cluster options)
SNAKEMAKE_CMD="snakemake \
    --use-singularity \
    --singularity-args '--bind /lustre:/lustre' \
    --cores $CORES \
    --jobs $JOBS \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
    --configfile config.yaml"

# Add dry run flag if requested
if [ "$DRYRUN" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --dry-run --printshellcmds"
fi

# Add force run if specified
if [ ! -z "$FORCERUN" ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --forcerun $FORCERUN"
fi

# Print the command
echo "Running Snakemake pipeline locally (no cluster)"
echo "Cores: $CORES, Jobs: $JOBS"
echo ""
echo "Command: $SNAKEMAKE_CMD"
echo ""

# Execute the command
eval $SNAKEMAKE_CMD

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "Pipeline completed successfully!"
else
    echo ""
    echo "Pipeline failed. Check the logs for details."
    exit 1
fi
