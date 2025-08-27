#!/bin/bash

# Quick setup script for local FASTA files
# This helps integrate your assembled genomes into the pipeline

echo "=================================================="
echo "Setting up Local FASTA Files for Pipeline"
echo "=================================================="
echo ""

# Check if input directory exists
if [ ! -d "input" ]; then
    echo "Error: input/ directory not found!"
    echo "Please create it and add your FASTA files there."
    exit 1
fi

# List FASTA files in input directory
echo "Found FASTA files in input/:"
ls -la input/*.fa input/*.fasta 2>/dev/null || echo "No FASTA files found"
echo ""

# Check for metadata file
METADATA_OPTIONS=(
    "input/meta_input_pipeline.csv"
    "../gapseq_pipeline_final/metadata.csv"
    "metadata.csv"
)

METADATA_FILE=""
for option in "${METADATA_OPTIONS[@]}"; do
    if [ -f "$option" ]; then
        METADATA_FILE="$option"
        echo "Using metadata file: $METADATA_FILE"
        break
    fi
done

if [ -z "$METADATA_FILE" ]; then
    echo "Error: No metadata file found!"
    echo "Please provide one of:"
    echo "  - input/meta_input_pipeline.csv"
    echo "  - metadata.csv"
    exit 1
fi

# Check if the FASTA file is already in metadata
echo ""
echo "Checking metadata for FASTA references..."

# Look for .fa or .fasta entries in metadata
grep -E "\.fa|\.fasta" "$METADATA_FILE" || echo "No FASTA files referenced in metadata"

# Ask if user wants to add FASTA file to metadata
echo ""
echo "Do you want to add your FASTA file to the metadata? (y/n)"
read -r response

if [[ "$response" == "y" ]]; then
    echo ""
    echo "Enter the sample name (filename with .fa extension):"
    echo "Example: Marinacidobacteraceae_bin.13_nanopore_reassembled.fa"
    read -r FASTA_NAME
    
    echo "Which media condition? (Taurine/Creatinine/Carnitine/Xylan/Chitin)"
    read -r MEDIA
    
    # Create a simple metadata entry
    echo ""
    echo "Add this line to your metadata CSV under the $MEDIA column:"
    echo "$FASTA_NAME"
    echo ""
    echo "Or create a simple metadata file:"
    
    cat > input/simple_metadata.csv << EOF
Comment,Taurine,Creatinine,Carnitine,Xylan,Chitin,Filetype
FASTA,$( [[ "$MEDIA" == "Taurine" ]] && echo "$FASTA_NAME" || echo ""),$( [[ "$MEDIA" == "Creatinine" ]] && echo "$FASTA_NAME" || echo ""),$( [[ "$MEDIA" == "Carnitine" ]] && echo "$FASTA_NAME" || echo ""),$( [[ "$MEDIA" == "Xylan" ]] && echo "$FASTA_NAME" || echo ""),$( [[ "$MEDIA" == "Chitin" ]] && echo "$FASTA_NAME" || echo ""),.fa
EOF
    
    echo "Created: input/simple_metadata.csv"
    METADATA_FILE="input/simple_metadata.csv"
fi

# Run the updated metadata parser
echo ""
echo "=================================================="
echo "Preparing metadata with local FASTA support..."
echo "=================================================="

if [ -f "prepare_metadata_v2.py" ]; then
    python prepare_metadata_v2.py "$METADATA_FILE" --input-dir input --validate
else
    echo "Warning: prepare_metadata_v2.py not found"
    echo "Using original parser (may not detect local files):"
    python prepare_metadata.py "$METADATA_FILE" --validate
fi

echo ""
echo "=================================================="
echo "Setup Complete!"
echo "=================================================="
echo ""
echo "Your FASTA file will be processed with:"
echo "  ✓ GTDB-Tk (taxonomy)"
echo "  ✓ Prodigal (gene prediction)"
echo "  ✓ Gapseq (metabolic modeling)"
echo ""
echo "It will SKIP:"
echo "  ✗ Trimming"
echo "  ✗ Assembly"
echo "  ✗ Binning"
echo ""
echo "To run the pipeline:"
echo "  ./run_pipeline.sh -c 32 -j 10"
echo ""
echo "To process ONLY your FASTA file:"
echo "  snakemake --use-singularity -c 8 results/09_gapseq/fasta/$FASTA_NAME/"
