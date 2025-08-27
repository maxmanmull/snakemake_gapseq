#!/bin/bash

# Container Checker and Setup Script
# This script checks for required containers and helps download missing ones

SINGULARITY_DIR="/lustre/BIF/nobackup/mulle088/singularity"

echo "=================================================="
echo "Checking Singularity Containers for Pipeline"
echo "=================================================="
echo ""

# Define required containers
declare -A CONTAINERS=(
    ["trimgalore"]="quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img"
    ["spades"]="quay.io-biocontainers-spades-3.15.5--h95f258a_0.img"
    ["quast"]="quay.io-biocontainers-quast-5.2.0--py39pl5321h2add14b_1.img"
    ["concoct"]="concoct_1.1.0--py312h71dcd68_7.sif"
    ["gtdbtk"]="quay.io-biocontainers-gtdbtk-2.4.1--pyhdfd78af_1.img"
    ["prodigal"]="prodigal_2.6.3--h779adbc_3.sif"
    ["gapseq"]="quay.io-biocontainers-gapseq-1.4.0--h9ee0642_1.img"
    ["multiqc"]="quay.io-biocontainers-multiqc-1.13--pyhdfd78af_0.img"
    ["checkm"]="checkm-genome_1.2.4--pyhdfd78af_2.sif"
)

# Check each container
MISSING_CONTAINERS=()
FOUND_CONTAINERS=()

for tool in "${!CONTAINERS[@]}"; do
    container_file="${SINGULARITY_DIR}/${CONTAINERS[$tool]}"
    if [ -f "$container_file" ]; then
        echo "✓ $tool: ${CONTAINERS[$tool]} found"
        FOUND_CONTAINERS+=("$tool")
    else
        echo "✗ $tool: ${CONTAINERS[$tool]} NOT FOUND"
        MISSING_CONTAINERS+=("$tool")
    fi
done

echo ""
echo "Summary:"
echo "  Found: ${#FOUND_CONTAINERS[@]} containers"
echo "  Missing: ${#MISSING_CONTAINERS[@]} containers"

# If CheckM is missing, provide download instructions
if [[ " ${MISSING_CONTAINERS[@]} " =~ " checkm " ]]; then
    echo ""
    echo "=================================================="
    echo "CheckM Container Setup Instructions"
    echo "=================================================="
    echo ""
    echo "The CheckM container is missing. To download it, run:"
    echo ""
    echo "cd $SINGULARITY_DIR"
    echo "singularity pull docker://quay.io/biocontainers/checkm-genome:1.2.4--pyhdfd78af_2"
    echo ""
    echo "This will create: checkm-genome_1.2.4--pyhdfd78af_2.sif"
    echo ""
fi

# Check for other missing containers
if [ ${#MISSING_CONTAINERS[@]} -gt 0 ]; then
    echo ""
    echo "=================================================="
    echo "Missing Container Download Commands"
    echo "=================================================="
    echo ""
    echo "To download missing containers, run these commands:"
    echo ""
    
    for tool in "${MISSING_CONTAINERS[@]}"; do
        case $tool in
            "checkm")
                # Already handled above
                ;;
            *)
                # For other containers, try to construct download command
                container_name="${CONTAINERS[$tool]}"
                if [[ $container_name == *.sif ]]; then
                    base_name="${container_name%.sif}"
                    echo "# $tool:"
                    echo "cd $SINGULARITY_DIR"
                    echo "singularity pull $container_name docker://quay.io/biocontainers/$base_name"
                    echo ""
                elif [[ $container_name == *.img ]]; then
                    base_name="${container_name%.img}"
                    echo "# $tool:"
                    echo "cd $SINGULARITY_DIR"
                    echo "singularity pull $container_name docker://$base_name"
                    echo ""
                fi
                ;;
        esac
    done
fi

# Check GTDB-Tk database
echo ""
echo "=================================================="
echo "Checking GTDB-Tk Database"
echo "=================================================="
echo ""

GTDBTK_DB="/lustre/BIF/nobackup/mulle088/gapseq_pipeline_final/gtdbtk_db/release226"

if [ -d "$GTDBTK_DB" ]; then
    echo "✓ GTDB-Tk database found at: $GTDBTK_DB"
    
    # Check for key database files
    if [ -d "$GTDBTK_DB/taxonomy" ] && [ -d "$GTDBTK_DB/fastani" ]; then
        echo "✓ Database structure appears complete"
    else
        echo "⚠ Database may be incomplete. Key directories missing."
        echo "  Please ensure the database is fully extracted."
    fi
else
    echo "✗ GTDB-Tk database NOT FOUND at: $GTDBTK_DB"
    echo ""
    echo "To download the GTDB-Tk database:"
    echo "1. Create directory: mkdir -p $GTDBTK_DB"
    echo "2. Download: wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package_data.tar.gz"
    echo "3. Extract: tar -xzf gtdbtk_package_data.tar.gz -C $GTDBTK_DB"
fi

# Check media files
echo ""
echo "=================================================="
echo "Checking Media Files"
echo "=================================================="
echo ""

MEDIA_DIR="/lustre/BIF/nobackup/mulle088/gapseq_pipeline_final/media"
MEDIA_FILES=("Taurine_media.csv" "Creatinine_media.csv" "Carnitine_media.csv" "Xylan_media.csv" "Chitin_media.csv")

MISSING_MEDIA=()
for media_file in "${MEDIA_FILES[@]}"; do
    if [ -f "$MEDIA_DIR/$media_file" ]; then
        echo "✓ $media_file found"
    else
        echo "✗ $media_file NOT FOUND"
        MISSING_MEDIA+=("$media_file")
    fi
done

if [ ${#MISSING_MEDIA[@]} -gt 0 ]; then
    echo ""
    echo "⚠ Missing media files. Please ensure all media definition files are in:"
    echo "  $MEDIA_DIR"
fi

# Final summary
echo ""
echo "=================================================="
echo "Setup Status Summary"
echo "=================================================="
echo ""

if [ ${#MISSING_CONTAINERS[@]} -eq 0 ] && [ ${#MISSING_MEDIA[@]} -eq 0 ] && [ -d "$GTDBTK_DB" ]; then
    echo "✓ All requirements satisfied! Pipeline is ready to run."
else
    echo "⚠ Some requirements are missing. Please address the issues above before running the pipeline."
    echo ""
    echo "Missing components:"
    [ ${#MISSING_CONTAINERS[@]} -gt 0 ] && echo "  - ${#MISSING_CONTAINERS[@]} container(s)"
    [ ${#MISSING_MEDIA[@]} -gt 0 ] && echo "  - ${#MISSING_MEDIA[@]} media file(s)"
    [ ! -d "$GTDBTK_DB" ] && echo "  - GTDB-Tk database"
fi

echo ""
echo "=================================================="
