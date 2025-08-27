#!/bin/bash

# Pipeline Visualization Script
# Generates visual representations of the pipeline structure

echo "=================================================="
echo "Pipeline Visualization Tool"
echo "=================================================="
echo ""

# Check if Snakemake and Graphviz are available
if ! command -v snakemake &> /dev/null; then
    echo "Error: Snakemake is not installed or not in PATH"
    echo "Please load the Snakemake module or install it"
    exit 1
fi

if ! command -v dot &> /dev/null; then
    echo "Warning: Graphviz (dot) is not installed"
    echo "Installing with conda/pip or loading module is recommended"
    echo "Trying to continue anyway..."
fi

# Function to generate DAG
generate_dag() {
    echo "Generating DAG (Directed Acyclic Graph) of all jobs..."
    
    if [ "$1" == "sample" ]; then
        # Generate DAG for a single sample
        echo "Enter sample ID (e.g., 12, 459, or Marinacidobacteraceae.fa):"
        read -r SAMPLE
        
        snakemake \
            --dag \
            results/08_gapseq/${SAMPLE}/${SAMPLE}_model.RDS \
            2>/dev/null | dot -Tpdf > pipeline_dag_${SAMPLE}.pdf
        
        if [ $? -eq 0 ]; then
            echo "✓ DAG for sample ${SAMPLE} saved to: pipeline_dag_${SAMPLE}.pdf"
        else
            echo "✗ Failed to generate DAG for sample ${SAMPLE}"
        fi
    else
        # Generate full DAG
        snakemake --dag 2>/dev/null | dot -Tpdf > pipeline_dag_full.pdf
        
        if [ $? -eq 0 ]; then
            echo "✓ Full DAG saved to: pipeline_dag_full.pdf"
        else
            echo "✗ Failed to generate full DAG"
        fi
    fi
}

# Function to generate rule graph
generate_rulegraph() {
    echo "Generating rule graph (simplified view)..."
    
    snakemake --rulegraph 2>/dev/null | dot -Tpdf > pipeline_rulegraph.pdf
    
    if [ $? -eq 0 ]; then
        echo "✓ Rule graph saved to: pipeline_rulegraph.pdf"
    else
        echo "✗ Failed to generate rule graph"
    fi
}

# Function to generate file graph
generate_filegraph() {
    echo "Generating file dependency graph..."
    
    if [ "$1" == "sample" ]; then
        echo "Enter sample ID:"
        read -r SAMPLE
        
        snakemake \
            --filegraph \
            results/08_gapseq/${SAMPLE}/${SAMPLE}_model.RDS \
            2>/dev/null | dot -Tpdf > pipeline_filegraph_${SAMPLE}.pdf
        
        if [ $? -eq 0 ]; then
            echo "✓ File graph for sample ${SAMPLE} saved to: pipeline_filegraph_${SAMPLE}.pdf"
        else
            echo "✗ Failed to generate file graph for sample ${SAMPLE}"
        fi
    else
        snakemake --filegraph 2>/dev/null | dot -Tpdf > pipeline_filegraph_full.pdf
        
        if [ $? -eq 0 ]; then
            echo "✓ Full file graph saved to: pipeline_filegraph_full.pdf"
        else
            echo "✗ Failed to generate full file graph"
        fi
    fi
}

# Function to show execution plan
show_execution_plan() {
    echo "Showing execution plan (dry run)..."
    echo ""
    
    if [ "$1" == "sample" ]; then
        echo "Enter sample ID:"
        read -r SAMPLE
        
        echo "Execution plan for sample ${SAMPLE}:"
        echo "======================================"
        snakemake \
            --dryrun \
            --printshellcmds \
            results/08_gapseq/${SAMPLE}/${SAMPLE}_model.RDS \
            2>/dev/null | grep -E "^rule |^    " | head -50
    else
        echo "Full execution plan (first 50 lines):"
        echo "======================================"
        snakemake --dryrun --printshellcmds 2>/dev/null | grep -E "^rule |^    " | head -50
    fi
}

# Function to show rule details
show_rule_details() {
    echo "Available rules in the pipeline:"
    echo "================================"
    
    grep "^rule " snakefile | sed 's/rule /  - /' | sed 's/://'
    
    echo ""
    echo "Enter rule name to see details (or 'skip' to continue):"
    read -r RULE
    
    if [ "$RULE" != "skip" ] && [ ! -z "$RULE" ]; then
        echo ""
        echo "Details for rule: $RULE"
        echo "========================"
        sed -n "/^rule $RULE:/,/^rule /p" snakefile | head -n -1
    fi
}

# Function to check parallelization potential
check_parallelization() {
    echo "Checking parallelization potential..."
    echo "====================================="
    
    echo ""
    echo "Rules that can run in parallel after genome preparation:"
    echo "  ✓ gtdbtk_classify - Taxonomic classification"
    echo "  ✓ prodigal_genes - Gene prediction"
    echo "  ✓ gapseq_find_pathways - Pathway identification"
    echo "  ✓ gapseq_find_transport - Transporter identification"
    
    echo ""
    echo "Sequential dependencies:"
    echo "  1. trim_reads → assemble_reads → map_reads_to_assembly"
    echo "  2. map_reads_to_assembly → calculate_coverage → binning_concoct"
    echo "  3. binning_concoct → checkm_bins"
    echo "  4. gapseq_find_pathways + gapseq_find_transport → gapseq_draft"
    echo "  5. gapseq_draft → gapseq_fill"
}

# Main menu
while true; do
    echo ""
    echo "=================================================="
    echo "Pipeline Visualization Options"
    echo "=================================================="
    echo "1. Generate full DAG (all samples)"
    echo "2. Generate DAG for specific sample"
    echo "3. Generate rule graph (simplified view)"
    echo "4. Generate file dependency graph"
    echo "5. Show execution plan (dry run)"
    echo "6. Show rule details"
    echo "7. Check parallelization potential"
    echo "8. Generate all visualizations"
    echo "9. Exit"
    echo ""
    echo "Choose an option (1-9):"
    read -r OPTION
    
    case $OPTION in
        1)
            generate_dag "full"
            ;;
        2)
            generate_dag "sample"
            ;;
        3)
            generate_rulegraph
            ;;
        4)
            echo "Generate for (a)ll samples or (s)pecific sample?"
            read -r CHOICE
            if [ "$CHOICE" == "s" ]; then
                generate_filegraph "sample"
            else
                generate_filegraph "full"
            fi
            ;;
        5)
            echo "Show plan for (a)ll samples or (s)pecific sample?"
            read -r CHOICE
            if [ "$CHOICE" == "s" ]; then
                show_execution_plan "sample"
            else
                show_execution_plan "full"
            fi
            ;;
        6)
            show_rule_details
            ;;
        7)
            check_parallelization
            ;;
        8)
            echo "Generating all visualizations..."
            generate_dag "full"
            generate_rulegraph
            generate_filegraph "full"
            check_parallelization
            echo ""
            echo "All visualizations completed!"
            ;;
        9)
            echo "Exiting..."
            exit 0
            ;;
        *)
            echo "Invalid option. Please choose 1-9."
            ;;
    esac
done
