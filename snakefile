import os
import pandas as pd
from pathlib import Path
import glob

# Configuration
configfile: "config.yaml"

# Load FASTA paths if available
if os.path.exists("config_paths.yaml"):
    configfile: "config_paths.yaml"

# Paths
BASE_DIR = config.get("base_dir", "/lustre/BIF/nobackup/mulle088")
SINGULARITY_DIR = config.get("singularity_dir", f"{BASE_DIR}/singularity")
DATA_DIR = config.get("data_dir", f"{BASE_DIR}/data_rp/90-1196727703/00_fastq")
MEDIA_DIR = config.get("media_dir", f"{BASE_DIR}/gapseq_pipeline_final/media")
INPUT_DIR = config.get("input_dir", "input")  # Local input directory
GTDBTK_DB = config.get("gtdbtk_db", f"{BASE_DIR}/gapseq_pipeline_final/gtdbtk_db/release226")
METADATA_FILE = config.get("metadata", f"{BASE_DIR}/gapseq_pipeline_final/metadata.csv")

# Get FASTA paths from config
FASTA_PATHS = config.get("fasta_paths", {})

# Output directories
RESULTS_DIR = config.get("results_dir", "results")
TRIMMED_DIR = f"{RESULTS_DIR}/01_trimmed"
ASSEMBLY_DIR = f"{RESULTS_DIR}/02_assembly"
QC_DIR = f"{RESULTS_DIR}/03_qc"
BINNING_DIR = f"{RESULTS_DIR}/04_binning"
BINS_DIR = f"{RESULTS_DIR}/05_bins"
CHECKM_DIR = f"{RESULTS_DIR}/06_checkm"
GTDBTK_DIR = f"{RESULTS_DIR}/07_gtdbtk"
PRODIGAL_DIR = f"{RESULTS_DIR}/08_prodigal"
GAPSEQ_DIR = f"{RESULTS_DIR}/09_gapseq"
FINAL_QC_DIR = f"{RESULTS_DIR}/10_final_qc"

# Parse metadata to get sample-media mapping
metadata_df = pd.read_csv(METADATA_FILE)
SAMPLE_MEDIA = {}
MEDIA_TYPES = ["Taurine", "Creatinine", "Carnitine", "Xylan", "Chitin"]

# Build sample-media mapping from metadata
for media in MEDIA_TYPES:
    if media in metadata_df.columns:
        samples = metadata_df[media].dropna().tolist()
        for sample in samples:
            if str(sample).endswith('.fa') or str(sample).endswith('.fasta'):
                # This is a pre-assembled genome
                SAMPLE_MEDIA[sample] = media
            elif str(sample).isdigit():
                # This is a fastq sample ID
                SAMPLE_MEDIA[sample] = media

# Get all fastq samples and fasta samples
FASTQ_SAMPLES = [s for s in SAMPLE_MEDIA.keys() if str(s).isdigit()]
FASTA_SAMPLES = [s for s in SAMPLE_MEDIA.keys() if s.endswith(('.fa', '.fasta'))]

# For binned samples, we'll create a wildcard for bins
# This will be dynamically determined after binning

# Containers
CONTAINERS = {
    "trimgalore": f"{SINGULARITY_DIR}/quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img",
    "spades": f"{SINGULARITY_DIR}/quay.io-biocontainers-spades-3.15.5--h95f258a_0.img",
    "bwa": f"{SINGULARITY_DIR}/bwa_0.7.17--hed695b0_7.sif",
    "samtools": f"{SINGULARITY_DIR}/samtools_1.16.1--h6899075_1.sif",
    "quast": f"{SINGULARITY_DIR}/quay.io-biocontainers-quast-5.2.0--py39pl5321h2add14b_1.img",
    "concoct": f"{SINGULARITY_DIR}/concoct_1.1.0--py312h71dcd68_7.sif",
    "checkm": f"{SINGULARITY_DIR}/checkm-genome_1.2.4--pyhdfd78af_2.sif",
    "gtdbtk": f"{SINGULARITY_DIR}/quay.io-biocontainers-gtdbtk-2.4.1--pyhdfd78af_1.img",
    "prodigal": f"{SINGULARITY_DIR}/prodigal_2.6.3--h779adbc_3.sif",
    "gapseq": f"{SINGULARITY_DIR}/quay.io-biocontainers-gapseq-1.4.0--h9ee0642_1.img",
    "multiqc": f"{SINGULARITY_DIR}/quay.io-biocontainers-multiqc-1.13--pyhdfd78af_0.img"
}

# Default rule
rule all:
    input:
        # Assembly QC for FASTQ samples
        expand(f"{QC_DIR}/assembly/{{sample}}/report.txt", sample=FASTQ_SAMPLES),
        # Binning checkpoint for FASTQ samples
        expand(f"{BINNING_DIR}/{{sample}}/bins_completed.txt", sample=FASTQ_SAMPLES),
        # CheckM for bins
        expand(f"{CHECKM_DIR}/{{sample}}/checkm_results.txt", sample=FASTQ_SAMPLES),
        # Batch GTDB-Tk for ALL genomes (MAGs + FASTA)
        f"{GTDBTK_DIR}/batch_results_parsed.txt",
        # For each FASTQ sample, process all bins
        expand(f"{BINS_DIR}/{{sample}}/all_bins_processed.txt", sample=FASTQ_SAMPLES),
        # For FASTA samples, process directly (no GTDB-Tk here since it's in batch)
        expand(f"{GAPSEQ_DIR}/fasta/{{sample}}/{{sample}}.RDS", sample=FASTA_SAMPLES),
        # Final QC report
        f"{FINAL_QC_DIR}/final_report.html"

# ================== MODULE 1: Quality Control and Trimming ==================

rule trim_reads:
    """
    Trim raw reads using Trim-Galore
    """
    input:
        r1 = f"{DATA_DIR}/{{sample}}_R1_001.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        r1 = f"{TRIMMED_DIR}/{{sample}}_R1_001_val_1.fq.gz",  
        r2 = f"{TRIMMED_DIR}/{{sample}}_R2_001_val_2.fq.gz",  
        report1 = f"{TRIMMED_DIR}/{{sample}}_R1_001.fastq.gz_trimming_report.txt",
        report2 = f"{TRIMMED_DIR}/{{sample}}_R2_001.fastq.gz_trimming_report.txt"
    params:
        quality = config.get("trim_quality", 20),
        min_length = config.get("trim_min_length", 50)
    threads: 1
    container: CONTAINERS["trimgalore"]
    shell:
        """
        trim_galore \
            --quality {params.quality} \
            --length {params.min_length} \
            --paired \
            --gzip \
            --cores {threads} \
            --output_dir {TRIMMED_DIR} \
            {input.r1} {input.r2}
        """
# ================== MODULE 2: Assembly ==================

rule assemble_reads:
    """
    Assemble trimmed reads using SPAdes
    """
    input:
        r1 = f"{TRIMMED_DIR}/{{sample}}_R1_001_val_1.fq.gz",
        r2 = f"{TRIMMED_DIR}/{{sample}}_R2_001_val_2.fq.gz"
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta",
        scaffolds = f"{ASSEMBLY_DIR}/{{sample}}/scaffolds.fasta"
    params:
        outdir = f"{ASSEMBLY_DIR}/{{sample}}",
        kmers = config.get("spades_kmers", "21,33,55,77")
    threads: 1
    resources:
        mem_mb=120000,
        time="24:00:00",
        assembly=1
    container: CONTAINERS["spades"]
    shell:
        """
        spades.py \
            --isolate \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            -k {params.kmers} \
            -t {threads} \
        """

# ================== MODULE 3: Assembly QC ==================

rule quast_assembly:
    """
    Quality assessment of assemblies using QUAST
    """
    input:
        contigs = rules.assemble_reads.output.contigs
    output:
        report = f"{QC_DIR}/assembly/{{sample}}/report.txt"
    params:
        outdir = f"{QC_DIR}/assembly/{{sample}}"
    threads: 1
    container: CONTAINERS["quast"]
    shell:
        """
        quast.py \
            {input.contigs} \
            -o {params.outdir} \
            --threads {threads} \
            --min-contig 500
        """

# ================== MODULE 4: Read Mapping for Binning ==================

rule map_reads_to_assembly:
    """
    Map reads back to assembly for coverage calculation
    """
    input:
        contigs = rules.assemble_reads.output.contigs,
        r1 = f"{TRIMMED_DIR}/{{sample}}_R1_001_val_1.fq.gz",
        r2 = f"{TRIMMED_DIR}/{{sample}}_R2_001_val_2.fq.gz"
    output:
        bam = f"{BINNING_DIR}/{{sample}}/{{sample}}_mapped.sorted.bam",
        bai = f"{BINNING_DIR}/{{sample}}/{{sample}}_mapped.sorted.bam.bai"
    threads: 1
    shell:
        """
        # Index contigs
        singularity exec {CONTAINERS[bwa]} bwa index {input.contigs}
        
        # Map reads and pipe to samtools
        singularity exec {CONTAINERS[bwa]} bwa mem -t {threads} {input.contigs} {input.r1} {input.r2} | \
        singularity exec {CONTAINERS[samtools]} samtools sort -@ {threads} -o {output.bam}
        
        # Index BAM
        singularity exec {CONTAINERS[samtools]} samtools index -@ {threads} {output.bam}
        """

rule calculate_coverage:
    """
    Calculate coverage depth for binning using CONCOCT's proper tools
    """
    input:
        contigs = rules.assemble_reads.output.contigs,
        bam = rules.map_reads_to_assembly.output.bam
    output:
        coverage = f"{BINNING_DIR}/{{sample}}/coverage_table.tsv"
    params:
        outdir = f"{BINNING_DIR}/{{sample}}"
    threads: 1
    container: CONTAINERS["concoct"]  # Container BEFORE shell
    shell:
        """
        # First cut up the contigs and create BED file
        cut_up_fasta.py \
            {input.contigs} \
            -c 10000 -o 0 \
            --merge_last \
            -b {params.outdir}/contigs_10K.bed \
            > {params.outdir}/contigs_10K.fa
        
        # Then use CONCOCT's coverage table tool with BED and BAM
        concoct_coverage_table.py \
            {params.outdir}/contigs_10K.bed \
            {input.bam} \
            > {output.coverage}
        """
# ================== MODULE 5: Binning ==================
rule binning_concoct:
    """
    Bin contigs using CONCOCT
    """
    input:
        contigs = rules.assemble_reads.output.contigs,
        coverage = rules.calculate_coverage.output.coverage
    output:
        bins_dir = directory(f"{BINNING_DIR}/{{sample}}/bins"),
        clustering = f"{BINNING_DIR}/{{sample}}/concoct_clustering.csv",
        completed = f"{BINNING_DIR}/{{sample}}/bins_completed.txt"
    params:
        outdir = f"{BINNING_DIR}/{{sample}}"
    threads: 1
    resources:
        mem_mb=32000,
        time="12:00:00",
        binning=1
    container: CONTAINERS["concoct"]  # Container BEFORE shell
    shell:
        """
        # Create output directories
        mkdir -p {params.outdir}/concoct_output
        mkdir -p {output.bins_dir}
        
        # The cut contigs and bed file already exist from calculate_coverage rule
        CUT_CONTIGS={params.outdir}/contigs_10K.fa
        
        # Check if we have enough contigs
        if [ -f "$CUT_CONTIGS" ]; then
            NUM_CONTIGS=$(grep -c "^>" $CUT_CONTIGS 2>/dev/null || echo "0")
        else
            NUM_CONTIGS=0
        fi
        
        echo "Cut contigs: $NUM_CONTIGS"
        
        if [ "$NUM_CONTIGS" -lt 4 ]; then
            echo "Too few contigs for binning"
            echo "contig_id,cluster_id" > {output.clustering}
            echo "no_bins_found" > {output.completed}
        else
            # Run CONCOCT
            concoct \
                --composition_file $CUT_CONTIGS \
                --coverage_file {input.coverage} \
                -b {params.outdir}/concoct_output/ \
                -t {threads} \
                --iterations 500 \
                --seed 42 || true
            
            # Process results if successful
            if [ -f {params.outdir}/concoct_output/clustering_gt1000.csv ]; then
                # Merge clustering
                merge_cutup_clustering.py \
                    {params.outdir}/concoct_output/clustering_gt1000.csv \
                    > {output.clustering}
                
                # Extract bins
                extract_fasta_bins.py \
                    {input.contigs} \
                    {output.clustering} \
                    --output_path {output.bins_dir}
                
                # List bins
                ls {output.bins_dir}/*.fa > {output.completed} 2>/dev/null || echo "no_bins_found" > {output.completed}
            else
                echo "CONCOCT produced no clustering"
                echo "contig_id,cluster_id" > {output.clustering}
                echo "no_bins_found" > {output.completed}
            fi
        fi
        
        # Ensure outputs exist
        touch {output.clustering}
        touch {output.completed}
        """


# ================== MODULE 6: Initial Bin Quality Check ==================
# Fixed CheckM rule with correct column parsing and no bc dependency

rule checkm_bins:
    """
    Check quality of bins using CheckM and filter based on completeness/contamination
    Fixed: correct column numbers and no bc command
    """
    input:
        bins_dir = rules.binning_concoct.output.bins_dir,
        completed = rules.binning_concoct.output.completed
    output:
        results = f"{CHECKM_DIR}/{{sample}}/checkm_results.txt",
        filtered_bins = f"{CHECKM_DIR}/{{sample}}/good_bins.txt"
    params:
        outdir = f"{CHECKM_DIR}/{{sample}}",
        completeness = config.get("checkm_completeness_threshold", 50),
        contamination = config.get("checkm_contamination_threshold", 10)
    threads: 8
    container: CONTAINERS["checkm"]
    shell:
        """
        # Check if we have any bins
        if grep -q "\.fa$" {input.completed}; then
            # Run CheckM lineage workflow
            checkm lineage_wf \
                -t {threads} \
                -x fa \
                {input.bins_dir} \
                {params.outdir} \
                > {output.results}
            
            # Create quality assessment table
            checkm qa \
                -o 2 \
                -f {params.outdir}/checkm_table.tsv \
                --tab_table \
                {params.outdir}/lineage.ms \
                {params.outdir}
            
            # Filter bins based on quality thresholds
            echo "Filtering bins with thresholds: Completeness>={params.completeness}%, Contamination<={params.contamination}%"
            
            # Parse CheckM table with CORRECT columns (6=completeness, 7=contamination)
            > {output.filtered_bins}
            tail -n +2 {params.outdir}/checkm_table.tsv | while IFS=$'\t' read -r bin_id col2 col3 col4 col5 completeness contamination rest; do
                
                # Use awk for floating point comparison (bc not available)
                PASS=$(awk -v comp="$completeness" -v cont="$contamination" -v min_comp="{params.completeness}" -v max_cont="{params.contamination}" \
                       'BEGIN {{if(comp >= min_comp && cont <= max_cont) print "YES"; else print "NO"}}')
                
                if [ "$PASS" == "YES" ]; then
                    echo "$bin_id" >> {output.filtered_bins}
                    echo "  ✓ PASS: $bin_id (C:${{completeness}}%, Cont:${{contamination}}%)"
                else
                    echo "  ✗ FAIL: $bin_id (C:${{completeness}}%, Cont:${{contamination}}%)"
                fi
            done
            
            # Report summary
            TOTAL_BINS=$(ls {input.bins_dir}/*.fa 2>/dev/null | wc -l)
            GOOD_BINS=$(cat {output.filtered_bins} 2>/dev/null | wc -l)
            echo ""
            echo "CheckM Summary for {wildcards.sample}:"
            echo "  Total bins: $TOTAL_BINS"
            echo "  Passed QC: $GOOD_BINS"
            echo "  Failed QC: $((TOTAL_BINS - GOOD_BINS))"
            
            # If no bins pass, suggest processing whole assembly
            if [ "$GOOD_BINS" -eq 0 ]; then
                echo ""
                echo "  No bins passed quality thresholds for {wildcards.sample}"
                echo "  Consider processing the whole assembly as a single genome"
            fi
            
        else
            echo "No bins found for {wildcards.sample}" > {output.results}
            touch {output.filtered_bins}
        fi
        """
# ================== MODULE 7: Process Individual Bins (Parallel) ==================

checkpoint get_bins:
    """
    Get list of bins for a sample to process - handles no bins gracefully
    """
    input:
        bins_dir = f"{BINNING_DIR}/{{sample}}/bins",
        filtered = f"{CHECKM_DIR}/{{sample}}/good_bins.txt"
    output:
        bins_list = f"{BINS_DIR}/{{sample}}/bins_to_process.txt"
    shell:
        """
        mkdir -p {BINS_DIR}/{wildcards.sample}
        
        # Check if we have any good quality bins
        if [ -s {input.filtered} ]; then
            # Process good quality bins
            for bin_id in $(cat {input.filtered}); do
                if [ -f {input.bins_dir}/${{bin_id}}.fa ]; then
                    echo "${{bin_id}}" >> {output.bins_list}
                    cp {input.bins_dir}/${{bin_id}}.fa {BINS_DIR}/{wildcards.sample}/
                fi
            done
        fi
        
        # If no bins were found, create a marker
        if [ ! -s {output.bins_list} ]; then
            echo "no_bins_found" > {output.bins_list}
        fi
        """

# Batch GTDB-Tk for all samples
rule prepare_gtdbtk_batch:
    """
    Prepare all genomes for batch GTDB-Tk processing - fixed version
    """
    input:
        # Make bins_lists optional - don't wait for all samples
        bins_lists = expand(f"{BINS_DIR}/{{sample}}/bins_to_process.txt", 
                          sample=FASTQ_SAMPLES),
        # Include FASTA samples
        fasta_files = [FASTA_PATHS.get(s, f"{INPUT_DIR}/{s}") for s in FASTA_SAMPLES]
    output:
        genome_dir = directory(f"{GTDBTK_DIR}/batch_genomes"),
        genome_list = f"{GTDBTK_DIR}/batch_genomes_list.txt"
    run:
        import os
        import shutil
        
        # Create directory for all genomes
        os.makedirs(output.genome_dir, exist_ok=True)
        genome_paths = []
        
        # Copy all MAG bins that exist
        for sample in FASTQ_SAMPLES:
            # Check both 05_bins and 04_binning for bins
            bins_processed = False
            
            # Try 05_bins first (quality filtered)
            bins_dir_filtered = f"{BINS_DIR}/{sample}"
            if os.path.exists(bins_dir_filtered):
                for bin_file in glob.glob(f"{bins_dir_filtered}/*.fa"):
                    if os.path.exists(bin_file):
                        bin_id = os.path.basename(bin_file).replace('.fa', '')
                        dst = f"{output.genome_dir}/{sample}_{bin_id}.fa"
                        shutil.copy(bin_file, dst)
                        genome_paths.append(dst)
                        bins_processed = True
            
            # If no quality bins, check raw bins
            if not bins_processed:
                bins_dir_raw = f"{BINNING_DIR}/{sample}/bins"
                if os.path.exists(bins_dir_raw):
                    for bin_file in glob.glob(f"{bins_dir_raw}/*.fa"):
                        if os.path.exists(bin_file):
                            bin_id = os.path.basename(bin_file).replace('.fa', '')
                            dst = f"{output.genome_dir}/{sample}_{bin_id}.fa"
                            shutil.copy(bin_file, dst)
                            genome_paths.append(dst)
        
        # Copy all FASTA samples
        for fasta_file in input.fasta_files:
            if os.path.exists(fasta_file):
                basename = os.path.basename(fasta_file)
                dst = f"{output.genome_dir}/{basename}"
                shutil.copy(fasta_file, dst)
                genome_paths.append(dst)
        
        # Write list of all genomes
        with open(output.genome_list, 'w') as f:
            if genome_paths:
                for path in genome_paths:
                    f.write(f"{path}\n")
            else:
                f.write("# No genomes found\n")
        
        print(f"Prepared {len(genome_paths)} genomes for GTDB-Tk")


rule gtdbtk_classify_batch:
    """
    Batch taxonomic classification using GTDB-Tk for ALL genomes at once
    """
    input:
        genome_dir = rules.prepare_gtdbtk_batch.output.genome_dir,
        genome_list = rules.prepare_gtdbtk_batch.output.genome_list
    output:
        summary = f"{GTDBTK_DIR}/batch_results/gtdbtk.bac120.summary.tsv",
        done = f"{GTDBTK_DIR}/batch_results/gtdbtk_complete.txt"
    params:
        outdir = f"{GTDBTK_DIR}/batch_results",
        db_path = GTDBTK_DB
    threads: 10
    resources:
        mem_mb=128000,
        time="24:00:00",
        gtdbtk=1
    container: CONTAINERS["gtdbtk"]
    shell:
        """
        export GTDBTK_DATA_PATH={params.db_path}
        
        # Count genomes
        NUM_GENOMES=$(ls {input.genome_dir}/*.fa 2>/dev/null | wc -l || echo "0")
        echo "Processing $NUM_GENOMES genomes in batch mode"
        
        if [ "$NUM_GENOMES" -gt 0 ]; then
            # Run GTDB-Tk WITHOUT the invalid --batch_file argument
            gtdbtk classify_wf \
                --genome_dir {input.genome_dir} \
                --extension fa \
                --out_dir {params.outdir} \
                --cpus {threads} \
                --pplacer_cpus {threads} \
                --skip_ani_screen
            
            echo "GTDB-Tk batch processing complete" > {output.done}
        else
            echo "No genomes to process" > {output.done}
            touch {output.summary}
        fi
        """
# Parse batch results back to individual samples
rule parse_gtdbtk_results:
    """
    Parse batch GTDB-Tk results and distribute to individual genome folders
    """
    input:
        summary = rules.gtdbtk_classify_batch.output.summary,
        done = rules.gtdbtk_classify_batch.output.done
    output:
        parsed = f"{GTDBTK_DIR}/batch_results_parsed.txt"
    run:
        import pandas as pd
        import os
        import shutil
        
        # Read the batch results
        if os.path.exists(input.summary) and os.path.getsize(input.summary) > 0:
            df = pd.read_csv(input.summary, sep='\t')
            
            # Parse and copy results to individual directories
            for _, row in df.iterrows():
                genome_id = row['user_genome']
                
                # Determine if it's a bin or FASTA sample
                if '_' in genome_id and genome_id.split('_')[0] in FASTQ_SAMPLES:
                    # It's a bin from a FASTQ sample
                    sample = genome_id.split('_')[0]
                    bin_id = '_'.join(genome_id.split('_')[1:]).replace('.fa', '')
                    
                    # Create directory and copy result
                    outdir = f"{GTDBTK_DIR}/bins/{sample}/{bin_id}"
                    os.makedirs(outdir, exist_ok=True)
                    
                    # Create individual summary file
                    row.to_frame().T.to_csv(
                        f"{outdir}/gtdbtk.bac120.summary.tsv",
                        sep='\t', index=False
                    )
                else:
                    # It's a FASTA sample
                    sample = genome_id.replace('.fa', '').replace('.fasta', '')
                    
                    # Create directory and copy result
                    outdir = f"{GTDBTK_DIR}/fasta/{sample}"
                    os.makedirs(outdir, exist_ok=True)
                    
                    # Create individual summary file
                    row.to_frame().T.to_csv(
                        f"{outdir}/gtdbtk.bac120.summary.tsv",
                        sep='\t', index=False
                    )
        
        with open(output.parsed, 'w') as f:
            f.write("Batch GTDB-Tk results parsed and distributed\n")

# Individual rules for other analyses (Prodigal, Gapseq) remain the same but without individual GTDB-Tk

rule prodigal_genes_bin:
    """
    Predict genes in individual bin using Prodigal
    """
    input:
        bin = f"{BINS_DIR}/{{sample}}/{{bin_id}}.fa"
    output:
        proteins = f"{PRODIGAL_DIR}/bins/{{sample}}/{{bin_id}}/proteins.faa",
        genes = f"{PRODIGAL_DIR}/bins/{{sample}}/{{bin_id}}/genes.fna",
        gff = f"{PRODIGAL_DIR}/bins/{{sample}}/{{bin_id}}/genes.gff"
    resources:
        mem_mb=4000,
        time="02:00:00",
        mag_analysis=1
    container: CONTAINERS["prodigal"]
    shell:
        """
        prodigal \
            -i {input.bin} \
            -a {output.proteins} \
            -d {output.genes} \
            -f gff \
            -o {output.gff} \
        """

rule gapseq_find_pathways_bin:
    """
    Find metabolic pathways in individual bin
    """
    input:
        bin = f"{BINS_DIR}/{{sample}}/{{bin_id}}.fa"
    output:
        pathways = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}/{{bin_id}}-all-Pathways.tbl",
        reactions = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}/{{bin_id}}-all-Reactions.tbl"
    params:
        outdir = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}"
    threads: 1
    container: CONTAINERS["gapseq"]
    shell:
        """
        cd {params.outdir}
        gapseq find -p all {input.bin}
        """

rule gapseq_find_transport_bin:
    """
    Find transporters in individual bin
    """
    input:
        bin = f"{BINS_DIR}/{{sample}}/{{bin_id}}.fa"
    output:
        transport = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}/{{bin_id}}-Transporter.tbl"
    params:
        outdir = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}"
    threads: 1
    container: CONTAINERS["gapseq"]
    shell:
        """
        cd {params.outdir}
        gapseq find-transport {input.bin}
        """

rule gapseq_draft_bin:
    """
    Create draft metabolic model for individual bin
    """
    input:
        bin = f"{BINS_DIR}/{{sample}}/{{bin_id}}.fa",
        pathways = rules.gapseq_find_pathways_bin.output.pathways,
        reactions = rules.gapseq_find_pathways_bin.output.reactions,
        transport = rules.gapseq_find_transport_bin.output.transport
    output:
        draft = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}/{{bin_id}}-draft.RDS",
        rxnWeights = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}/{{bin_id}}-rxnWeights.RDS",
        rxnXgenes = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}/{{bin_id}}-rxnXgenes.RDS"
    params:
        outdir = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}"
    container: CONTAINERS["gapseq"]
    shell:
        """
        cd {params.outdir}
        gapseq draft \
            -r {input.reactions} \
            -t {input.transport} \
            -p {input.pathways} \
            -c {input.bin}
        """

rule gapseq_fill_bin:
    """
    Fill metabolic model with media-specific gap-filling
    """
    input:
        draft = rules.gapseq_draft_bin.output.draft,
        rxnWeights = rules.gapseq_draft_bin.output.rxnWeights,
        rxnXgenes = rules.gapseq_draft_bin.output.rxnXgenes,
        media = lambda wildcards: f"{MEDIA_DIR}/{SAMPLE_MEDIA[wildcards.sample]}_media.csv"
    output:
        model = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}/{{bin_id}}.RDS"
    params:
        outdir = f"{GAPSEQ_DIR}/bins/{{sample}}/{{bin_id}}"
    container: CONTAINERS["gapseq"]
    shell:
        """
        cd {params.outdir}
        gapseq fill \
            -m {input.draft} \
            -c {input.rxnWeights} \
            -g {input.rxnXgenes} \
            -n {input.media}
        """
# Aggregate rule for all bins of a sample

def get_all_bin_outputs(wildcards):
    """
    Get all expected outputs for bins of a sample - handles no bins case
    """
    checkpoint_output = checkpoints.get_bins.get(**wildcards).output[0]
    
    # Read the bins list file
    with open(checkpoint_output) as f:
        content = f.read().strip()
    
    # Check if there are no bins
    if content == "no_bins_found" or not content:
        # Return empty list - no bin processing needed
        return []
    
    # Process each bin
    outputs = []
    with open(checkpoint_output) as f:
        bin_ids = [line.strip() for line in f if line.strip() and line.strip() != "no_bins_found"]
    
    for bin_id in bin_ids:
        outputs.extend([
            f"{PRODIGAL_DIR}/bins/{wildcards.sample}/{bin_id}/proteins.faa",
            f"{GAPSEQ_DIR}/bins/{wildcards.sample}/{bin_id}/{bin_id}_model.RDS"
        ])
    
    return outputs


rule process_all_bins:
    """
    Process all bins for a sample - handles no bins case
    """
    input:
        get_all_bin_outputs,
        gtdbtk_parsed = f"{GTDBTK_DIR}/batch_results_parsed.txt"
    output:
        f"{BINS_DIR}/{{sample}}/all_bins_processed.txt"
    shell:
        """
        echo "Bins processing complete for {wildcards.sample}" > {output}
        
        # Count processed bins
        if [ -d {GAPSEQ_DIR}/bins/{wildcards.sample} ]; then
            NUM_BINS=$(ls {GAPSEQ_DIR}/bins/{wildcards.sample}/*/{{*}}_model.RDS 2>/dev/null | wc -l || echo "0")
        else
            NUM_BINS=0
        fi
        
        echo "Processed bins: $NUM_BINS" >> {output}
        """
# ================== MODULE 8: Process FASTA samples (no individual GTDB-Tk) ==================

rule prodigal_genes_fasta:
    """
    Predict genes for pre-assembled genomes
    """
    input:
        genome = lambda wildcards: FASTA_PATHS.get(wildcards.sample, f"{INPUT_DIR}/{wildcards.sample}")
    output:
        proteins = f"{PRODIGAL_DIR}/fasta/{{sample}}/proteins.faa",
        genes = f"{PRODIGAL_DIR}/fasta/{{sample}}/genes.fna",
        gff = f"{PRODIGAL_DIR}/fasta/{{sample}}/genes.gff"
    container: CONTAINERS["prodigal"]
    shell:
        """
        prodigal \
            -i {input.genome} \
            -a {output.proteins} \
            -d {output.genes} \
            -f gff \
            -o {output.gff} \
            -p single
        """

rule gapseq_find_pathways_fasta:
    """
    Find pathways for pre-assembled genomes
    """
    input:
        genome = lambda wildcards: FASTA_PATHS.get(wildcards.sample, f"{INPUT_DIR}/{wildcards.sample}")
    output:
        pathways = f"{GAPSEQ_DIR}/fasta/{{sample}}/{{sample}}-all-Pathways.tbl",
        reactions = f"{GAPSEQ_DIR}/fasta/{{sample}}/{{sample}}-all-Reactions.tbl"
    params:
        outdir = f"{GAPSEQ_DIR}/fasta/{{sample}}"
    threads: 1
    container: CONTAINERS["gapseq"]
    shell:
        """
        cd {params.outdir}
        gapseq find -p all {input.genome}
        """

rule gapseq_find_transport_fasta:
    """
    Find transporters for pre-assembled genomes
    """
    input:
        genome = lambda wildcards: FASTA_PATHS.get(wildcards.sample, f"{INPUT_DIR}/{wildcards.sample}")
    output:
        transport = f"{GAPSEQ_DIR}/fasta/{{sample}}/{{sample}}-Transporter.tbl"
    params:
        outdir = f"{GAPSEQ_DIR}/fasta/{{sample}}"
    threads: 1
    container: CONTAINERS["gapseq"]
    shell:
        """
        cd {params.outdir}
        gapseq find-transport {input.genome}
        """

rule gapseq_draft_fasta:
    """
    Create draft model for pre-assembled genomes
    """
    input:
        genome = lambda wildcards: FASTA_PATHS.get(wildcards.sample, f"{INPUT_DIR}/{wildcards.sample}"),
        pathways = rules.gapseq_find_pathways_fasta.output.pathways,
        reactions = rules.gapseq_find_pathways_fasta.output.reactions,
        transport = rules.gapseq_find_transport_fasta.output.transport
    output:
        draft = f"{GAPSEQ_DIR}/fasta/{{sample}}/{{sample}}-draft.RDS",
        rxnWeights = f"{GAPSEQ_DIR}/fasta/{{sample}}/{{sample}}-rxnWeights.RDS",
        rxnXgenes = f"{GAPSEQ_DIR}/fasta/{{sample}}/{{sample}}-rxnXgenes.RDS"
    params:
        outdir = f"{GAPSEQ_DIR}/fasta/{{sample}}"
    container: CONTAINERS["gapseq"]
    shell:
        """
        cd {params.outdir}
        gapseq draft \
            -r {input.reactions} \
            -t {input.transport} \
            -p {input.pathways} \
            -c {input.genome}
        """

rule gapseq_fill_fasta:
    """
    Fill model for pre-assembled genomes
    """
    input:
        draft = rules.gapseq_draft_fasta.output.draft,
        rxnWeights = rules.gapseq_draft_fasta.output.rxnWeights,
        rxnXgenes = rules.gapseq_draft_fasta.output.rxnXgenes,
        media = lambda wildcards: f"{MEDIA_DIR}/{SAMPLE_MEDIA[wildcards.sample]}_media.csv"
    output:
        model = f"{GAPSEQ_DIR}/fasta/{{sample}}/{{sample}}.RDS"
    params:
        outdir = f"{GAPSEQ_DIR}/fasta/{{sample}}"
    container: CONTAINERS["gapseq"]
    shell:
        """
        cd {params.outdir}
        gapseq fill \
            -m {input.draft} \
            -c {input.rxnWeights} \
            -g {input.rxnXgenes} \
            -n {input.media}
        """
# ================== MODULE 9: Final Quality Control ==================

rule final_multiqc:
    """
    Generate final MultiQC report combining all results
    """
    input:
        expand(f"{TRIMMED_DIR}/{{sample}}_R1_001.fastq.gz_trimming_report.txt", sample=FASTQ_SAMPLES),
        expand(f"{QC_DIR}/assembly/{{sample}}/report.txt", sample=FASTQ_SAMPLES),
        expand(f"{CHECKM_DIR}/{{sample}}/checkm_results.txt", sample=FASTQ_SAMPLES)
    output:
        report = f"{FINAL_QC_DIR}/final_report.html"
    params:
        outdir = FINAL_QC_DIR
    container: CONTAINERS["multiqc"]
    shell:
        """
        multiqc \
            {TRIMMED_DIR} \
            {QC_DIR} \
            {CHECKM_DIR} \
            -o {params.outdir} \
            -n final_report.html \
            --force
        """
