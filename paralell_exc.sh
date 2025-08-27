#!/bin/bash

# Advanced Parallel Pipeline Runner
# Optimized for maximum parallelization across samples and within samples

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Default settings for parallelization
TOTAL_CORES=$(nproc)  # Get total available cores
RECOMMENDED_JOBS=$((TOTAL_CORES / 4))  # Conservative estimate
MAX_JOBS=100  # Maximum Snakemake jobs
CONFIG="config_parallel.yaml"
PROFILE=""
DRYRUN=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--cores)
            TOTAL_CORES="$2"
            shift 2
            ;;
        -j|--jobs)
            MAX_JOBS="$2"
            shift 2
            ;;
        --auto)
            # Automatic optimization based on system
            print_info "Auto-detecting optimal settings..."
            TOTAL_CORES=$(nproc)
            TOTAL_MEM=$(free -g | awk '/^Mem:/{print $2}')
            
            # Calculate optimal job numbers based on resources
            if [ $TOTAL_MEM -gt 500 ]; then
                # High memory system
                MAX_JOBS=$((TOTAL_CORES / 2))
                print_info "High-memory system detected: ${TOTAL_MEM}GB RAM"
            elif [ $TOTAL_MEM -gt 200 ]; then
                # Medium memory system
                MAX_JOBS=$((TOTAL_CORES / 3))
                print_info "Medium-memory system detected: ${TOTAL_MEM}GB RAM"
            else
                # Low memory system
                MAX_JOBS=$((TOTAL_CORES / 4))
                print_warning "Low-memory system detected: ${TOTAL_MEM}GB RAM"
            fi
            shift
            ;;
        -n|--dryrun)
            DRYRUN=true
            shift
            ;;
        --config)
            CONFIG="$2"
            shift 2
            ;;
        -h|--help)
            echo "Advanced Parallel Pipeline Runner"
            echo "================================"
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  -c, --cores NUM    Total cores available (default: auto-detect)"
            echo "  -j, --jobs NUM     Maximum parallel jobs (default: cores/4)"
            echo "  --auto            Auto-optimize based on system resources"
            echo "  -n, --dryrun      Perform a dry run"
            echo "  --config FILE     Config file (default: config_parallel.yaml)"
            echo "  -h, --help        Show this help message"
            echo ""
            echo "Examples:"
            echo "  # Auto-optimize for your system:"
            echo "  $0 --auto"
            echo ""
            echo "  # Manual settings for large cluster:"
            echo "  $0 -c 128 -j 50"
            echo ""
            echo "  # Conservative settings:"
            echo "  $0 -c 32 -j 8"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Calculate parallelization strategy
print_info "==================================================="
print_info "PARALLELIZATION STRATEGY"
print_info "==================================================="
print_info "Total cores available: $TOTAL_CORES"
print_info "Maximum parallel jobs: $MAX_JOBS"
print_info "Configuration file: $CONFIG"

# Estimate parallelization capacity
SAMPLES_FILE="sample_sheet.csv"
if [ -f "$SAMPLES_FILE" ]; then
    NUM_SAMPLES=$(tail -n +2 "$SAMPLES_FILE" | wc -l)
    print_info "Total samples to process: $NUM_SAMPLES"
    
    # Count FASTQ vs FASTA samples
    NUM_FASTQ=$(grep -c "fastq" "$SAMPLES_FILE" || echo "0")
    NUM_FASTA=$(grep -c "fasta" "$SAMPLES_FILE" || echo "0")
    print_info "  - FASTQ samples (need assembly): $NUM_FASTQ"
    print_info "  - FASTA samples (pre-assembled): $NUM_FASTA"
    
    # Estimate parallel execution
    if [ $MAX_JOBS -ge $NUM_SAMPLES ]; then
        print_info "Can process ALL samples in parallel!"
    else
        BATCHES=$((NUM_SAMPLES / MAX_JOBS + 1))
        print_info "Will process samples in $BATCHES batches"
    fi
else
    print_warning "Sample sheet not found. Run prepare_metadata.py first."
fi

echo ""

# Check resource availability
print_info "Checking resource availability..."

# Check memory
TOTAL_MEM=$(free -g | awk '/^Mem:/{print $2}')
AVAIL_MEM=$(free -g | awk '/^Mem:/{print $7}')
print_info "Memory: ${AVAIL_MEM}GB available of ${TOTAL_MEM}GB total"

if [ $AVAIL_MEM -lt 100 ]; then
    print_warning "Low available memory. Consider reducing parallel jobs."
    print_warning "Recommended: Use -j $((MAX_JOBS / 2))"
fi

# Check disk space
DISK_AVAIL=$(df -BG . | awk 'NR==2 {print $4}' | sed 's/G//')
print_info "Disk space available: ${DISK_AVAIL}GB"

if [ $DISK_AVAIL -lt 500 ]; then
    print_warning "Low disk space. Pipeline may need 50-100GB per sample."
fi

echo ""

# Load required modules
print_info "Loading required modules..."
module load snakemake/7.32.4 2>/dev/null || print_warning "Snakemake module not found"
module load singularity/3.11.3 2>/dev/null || print_warning "Singularity module not found"

# Create necessary directories
mkdir -p logs
mkdir -p results

# Build resource groups for Snakemake
# This prevents too many memory-intensive jobs from running simultaneously
RESOURCE_GROUPS=""

# Limit assembly jobs (very memory intensive)
RESOURCE_GROUPS="$RESOURCE_GROUPS assembly:3"

# Limit GTDB-Tk jobs (memory intensive)
RESOURCE_GROUPS="$RESOURCE_GROUPS gtdbtk:5"

# Limit binning jobs
RESOURCE_GROUPS="$RESOURCE_GROUPS binning:5"

# Allow many MAG analyses
RESOURCE_GROUPS="$RESOURCE_GROUPS mag_analysis:20"

print_info "Resource groups: $RESOURCE_GROUPS"

# Build the Snakemake command
SNAKEMAKE_CMD="snakemake \
    --use-singularity \
    --singularity-args '--bind /lustre:/lustre' \
    --cores $TOTAL_CORES \
    --jobs $MAX_JOBS \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
    --configfile $CONFIG"

# Add resource groups to prevent memory overload
if [ ! -z "$RESOURCE_GROUPS" ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --resources $RESOURCE_GROUPS"
fi

# Add dry run flag if requested
if [ "$DRYRUN" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --dry-run --printshellcmds"
    print_info "DRY RUN MODE - No actual processing will occur"
fi

# For SLURM clusters, add cluster configuration
if command -v sbatch &> /dev/null; then
    print_info "SLURM detected. Adding cluster configuration..."
    
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD \
        --cluster 'sbatch \
            --partition={params.partition} \
            --time={params.time} \
            --mem={params.mem_mb} \
            --cpus-per-task={threads} \
            --job-name=smk_{rule}_{wildcards} \
            --output=logs/{rule}_{wildcards}_%j.out \
            --error=logs/{rule}_{wildcards}_%j.err' \
        --default-resources \
            partition=normal \
            mem_mb=8000 \
            time='04:00:00' \
        --cluster-cancel 'scancel'"
fi

# Print execution plan
echo ""
print_info "==================================================="
print_info "EXECUTION PLAN"
print_info "==================================================="

if [ "$DRYRUN" = true ]; then
    echo "Running: $SNAKEMAKE_CMD" | fold -s -w 80
    echo ""
    print_info "Analyzing workflow..."
    
    # Count total jobs
    TOTAL_JOBS=$(snakemake --dryrun --quiet 2>&1 | grep -c "^Job" || echo "0")
    print_info "Total jobs to execute: $TOTAL_JOBS"
    
    # Show job summary by rule
    print_info "Jobs by rule:"
    snakemake --dryrun --quiet 2>&1 | grep "^rule" | sort | uniq -c | sort -rn | head -20
else
    echo "Command: $SNAKEMAKE_CMD" | fold -s -w 80
fi

echo ""
print_info "==================================================="
print_info "Starting pipeline execution..."
print_info "==================================================="
echo ""

# Create a timestamp for this run
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="logs/pipeline_run_${TIMESTAMP}.log"

# Execute the pipeline
if [ "$DRYRUN" = false ]; then
    # Run with real-time monitoring
    eval $SNAKEMAKE_CMD 2>&1 | tee "$LOG_FILE"
    
    # Check exit status
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo ""
        print_info "==================================================="
        print_info "✓ Pipeline completed successfully!"
        print_info "==================================================="
        
        # Generate summary
        print_info "Generating results summary..."
        if [ -f "summarize_results.py" ]; then
            python summarize_results.py \
                --results-dir results \
                --output "results_${TIMESTAMP}.txt" \
                --json
            
            print_info "Results summary saved to: results_${TIMESTAMP}.txt"
        fi
        
        # Show statistics
        print_info "Pipeline statistics:"
        echo "  - Total runtime: $(grep "total runtime" "$LOG_FILE" 2>/dev/null || echo "See log file")"
        echo "  - Log file: $LOG_FILE"
        
    else
        echo ""
        print_error "==================================================="
        print_error "✗ Pipeline failed. Check logs for details."
        print_error "==================================================="
        print_error "Log file: $LOG_FILE"
        
        # Show recent errors
        print_error "Recent errors:"
        grep -E "Error|ERROR|Failed" "$LOG_FILE" | tail -5
        
        exit 1
    fi
else
    # Dry run mode
    eval $SNAKEMAKE_CMD
fi

echo ""
print_info "Done!"
