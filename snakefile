import pandas as pd
from pathlib import Path
import glob
# Configuration
SINGULARITY_DIR = "/lustre/BIF/nobackup/mulle088/singularity"
DATA_DIR = "/lustre/BIF/nobackup/mulle088/data_rp/90-1196727703/00_fastq"
RESULTS_DIR = "results"

# Parse metadata to get FASTQ and FASTA samples
metadata = pd.read_csv("input/meta_input_pipeline.csv")
FASTQ_SAMPLES = []
FASTA_SAMPLES = []

# Extract unique sample IDs from all media columns
media_cols = ["Taurine", "Creatinine", "Carnitine", "Xylan", "Chitin"]
for col in media_cols:
    if col in metadata.columns:
        samples = metadata[col].dropna()
        for s in samples:
            s_str = str(s)
            if s_str.endswith('.fa'):
                if s_str not in FASTA_SAMPLES:
                    FASTA_SAMPLES.append(s_str)
            else:
                # Handle float conversion (459.0 -> 459)
                if '.' in s_str:
                    s_str = str(int(float(s_str)))
                if s_str not in FASTQ_SAMPLES:
                    FASTQ_SAMPLES.append(s_str)

print(f"Found {len(FASTQ_SAMPLES)} FASTQ samples: {FASTQ_SAMPLES}")
print(f"Found {len(FASTA_SAMPLES)} FASTA samples: {FASTA_SAMPLES}")

# Rule all - what we want to produce
rule all:
    input:
        expand(f"{RESULTS_DIR}/01_trimmed/{{sample}}_R1_001_val_1.fq.gz", sample=FASTQ_SAMPLES),
        expand(f"{RESULTS_DIR}/01_trimmed/{{sample}}_R2_001_val_2.fq.gz", sample=FASTQ_SAMPLES),
        expand(f"{RESULTS_DIR}/02_assembly/{{sample}}/contigs.fasta", sample=FASTQ_SAMPLES),
        expand(f"{RESULTS_DIR}/03_mapping/{{sample}}/{{sample}}.sorted.bam", sample=FASTQ_SAMPLES),
        expand(f"{RESULTS_DIR}/04_binning/{{sample}}/bins/done.txt", sample=FASTQ_SAMPLES),
        expand(f"{RESULTS_DIR}/05_checkm/{{sample}}/checkm_results.txt", sample=FASTQ_SAMPLES),
        expand(f"{RESULTS_DIR}/06_good_bins/{{sample}}/done.txt", sample=FASTQ_SAMPLES),
        f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/collected.txt",
        f"{RESULTS_DIR}/08_gtdbtk/done.txt",
        f"{RESULTS_DIR}/12_gapseq_models/multi_media_complete.txt",
        f"{RESULTS_DIR}/13_antismash/all_complete.txt",
        f"{RESULTS_DIR}/14_cazymes/all_done.txt",        
        f"{RESULTS_DIR}/15_qc_summary/mag_quality_master.tsv",
        f"{RESULTS_DIR}/15_qc_summary/mag_quality_summary.txt"
# Rule 1: Trim reads with Trim-Galore
rule trim_reads:
    input:
        r1 = f"{DATA_DIR}/{{sample}}_R1_001.fastq.gz",
        r2 = f"{DATA_DIR}/{{sample}}_R2_001.fastq.gz"
    output:
        r1 = f"{RESULTS_DIR}/01_trimmed/{{sample}}_R1_001_val_1.fq.gz",
        r2 = f"{RESULTS_DIR}/01_trimmed/{{sample}}_R2_001_val_2.fq.gz"
    params:
        outdir = f"{RESULTS_DIR}/01_trimmed"
    threads: 4
    container: f"{SINGULARITY_DIR}/quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img"
    shell:
        """
        mkdir -p {params.outdir}
        
        trim_galore \
            --quality 20 \
            --length 50 \
            --paired \
            --gzip \
            --cores {threads} \
            --output_dir {params.outdir} \
            {input.r1} {input.r2}
        """

# Rule 2: Assemble with SPAdes
rule assemble_spades:
    input:
        r1 = f"{RESULTS_DIR}/01_trimmed/{{sample}}_R1_001_val_1.fq.gz",
        r2 = f"{RESULTS_DIR}/01_trimmed/{{sample}}_R2_001_val_2.fq.gz"
    output:
        contigs = f"{RESULTS_DIR}/02_assembly/{{sample}}/contigs.fasta",
        scaffolds = f"{RESULTS_DIR}/02_assembly/{{sample}}/scaffolds.fasta"
    params:
        outdir = f"{RESULTS_DIR}/02_assembly/{{sample}}"
    threads: 16
    container: f"{SINGULARITY_DIR}/quay.io-biocontainers-spades-3.15.5--h95f258a_0.img"
    shell:
        """
        spades.py \
            --isolate \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            -k 21,33,55,77 \
            -t {threads} \
            -m 120
        """

# Rule 3: Map reads to assembly for binning
rule map_reads:
    input:
        contigs = f"{RESULTS_DIR}/02_assembly/{{sample}}/contigs.fasta",
        r1 = f"{RESULTS_DIR}/01_trimmed/{{sample}}_R1_001_val_1.fq.gz",
        r2 = f"{RESULTS_DIR}/01_trimmed/{{sample}}_R2_001_val_2.fq.gz"
    output:
        bam = f"{RESULTS_DIR}/03_mapping/{{sample}}/{{sample}}.sorted.bam",
        bai = f"{RESULTS_DIR}/03_mapping/{{sample}}/{{sample}}.sorted.bam.bai"
    params:
        outdir = f"{RESULTS_DIR}/03_mapping/{{sample}}",
        bwa_container = f"{SINGULARITY_DIR}/bwa_0.7.17--hed695b0_7.sif",
        samtools_container = f"{SINGULARITY_DIR}/samtools_1.16.1--h6899075_1.sif"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        
        # Index contigs
        singularity exec {params.bwa_container} \
            bwa index {input.contigs}
        
        # Map reads
        singularity exec {params.bwa_container} \
            bwa mem -t {threads} {input.contigs} {input.r1} {input.r2} | \
        singularity exec {params.samtools_container} \
            samtools sort -@ {threads} -o {output.bam}
        
        # Index BAM
        singularity exec {params.samtools_container} \
            samtools index -@ {threads} {output.bam}
        """

# Rule 4: Binning with CONCOCT
rule binning_concoct:
    input:
        contigs = f"{RESULTS_DIR}/02_assembly/{{sample}}/contigs.fasta",
        bam = f"{RESULTS_DIR}/03_mapping/{{sample}}/{{sample}}.sorted.bam"
    output:
        done = f"{RESULTS_DIR}/04_binning/{{sample}}/bins/done.txt"
    params:
        outdir = f"{RESULTS_DIR}/04_binning/{{sample}}",
        concoct_container = f"{SINGULARITY_DIR}/concoct_1.1.0--py312h71dcd68_7.sif"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}/bins
        
        # Cut up contigs
        singularity exec {params.concoct_container} \
            cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last \
            -b {params.outdir}/contigs_10K.bed > {params.outdir}/contigs_10K.fa
        
        # Generate coverage table
        singularity exec {params.concoct_container} \
            concoct_coverage_table.py {params.outdir}/contigs_10K.bed \
            {input.bam} > {params.outdir}/coverage_table.tsv
        
        # Run CONCOCT
        singularity exec {params.concoct_container} \
            concoct --composition_file {params.outdir}/contigs_10K.fa \
            --coverage_file {params.outdir}/coverage_table.tsv \
            -b {params.outdir}/ \
            -t {threads}
        
        # Merge clustering results
        singularity exec {params.concoct_container} \
            merge_cutup_clustering.py {params.outdir}/clustering_gt1000.csv \
            > {params.outdir}/clustering_merged.csv
        
        # Extract bins
        singularity exec {params.concoct_container} \
            extract_fasta_bins.py {input.contigs} \
            {params.outdir}/clustering_merged.csv \
            --output_path {params.outdir}/bins
        
        touch {output.done}
        """

# Rule 5: CheckM quality assessment
rule checkm_bins:
    input:
        done = f"{RESULTS_DIR}/04_binning/{{sample}}/bins/done.txt"
    output:
        results = f"{RESULTS_DIR}/05_checkm/{{sample}}/checkm_results.txt",
        table = f"{RESULTS_DIR}/05_checkm/{{sample}}/checkm_table.tsv"
    params:
        bins_dir = f"{RESULTS_DIR}/04_binning/{{sample}}/bins",
        outdir = f"{RESULTS_DIR}/05_checkm/{{sample}}",
        checkm_container = f"{SINGULARITY_DIR}/checkm-genome_1.2.4--pyhdfd78af_2.sif"
    threads: 8
    shell:
        """
        # Run CheckM lineage workflow
        singularity exec {params.checkm_container} \
            checkm lineage_wf \
            -t {threads} \
            -x fa \
            {params.bins_dir} \
            {params.outdir} \
            > {output.results}
        
        # Generate table output
        singularity exec {params.checkm_container} \
            checkm qa \
            -o 2 \
            -f {output.table} \
            --tab_table \
            {params.outdir}/lineage.ms \
            {params.outdir}
        """

# Rule 6: Filter bins >90% completeness
rule filter_good_bins:
    input:
        checkm_table = f"{RESULTS_DIR}/05_checkm/{{sample}}/checkm_table.tsv"
    output:
        done = f"{RESULTS_DIR}/06_good_bins/{{sample}}/done.txt",
        summary = f"{RESULTS_DIR}/06_good_bins/{{sample}}/summary.txt"
    params:
        bins_dir = f"{RESULTS_DIR}/04_binning/{{sample}}/bins",
        good_bins_dir = f"{RESULTS_DIR}/06_good_bins/{{sample}}",
        completeness_threshold = 90
    run:
        import os
        import shutil
        
        # Create output directory
        os.makedirs(params.good_bins_dir, exist_ok=True)
        
        good_bins = []
        total_bins = 0
        
        # Parse CheckM table
        with open(input.checkm_table, 'r') as f:
            lines = f.readlines()
            
        # Skip header line
        for line in lines[1:]:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    bin_id = parts[0]
                    completeness = float(parts[5])
                    contamination = float(parts[6])
                    
                    total_bins += 1
                    
                    # Check if completeness > threshold
                    if completeness > params.completeness_threshold:
                        # Try to copy the bin file
                        src_fa = f"{params.bins_dir}/{bin_id}.fa"
                        src_fasta = f"{params.bins_dir}/{bin_id}.fasta"
                        dst = f"{params.good_bins_dir}/{bin_id}.fa"
                        
                        if os.path.exists(src_fa):
                            shutil.copy(src_fa, dst)
                            good_bins.append(f"{bin_id}: Completeness={completeness:.1f}%, Contamination={contamination:.1f}%")
                        elif os.path.exists(src_fasta):
                            shutil.copy(src_fasta, dst)
                            good_bins.append(f"{bin_id}: Completeness={completeness:.1f}%, Contamination={contamination:.1f}%")
        
        # Write summary
        with open(output.summary, 'w') as f:
            for bin_info in good_bins:
                f.write(f"{bin_info}\n")
            f.write(f"\nTotal bins: {total_bins}\n")
            f.write(f"Good bins (>90% complete): {len(good_bins)}\n")
        
        # Create done file
        with open(output.done, 'w') as f:
            f.write("Filtering complete\n")

# Rule 7: Collect all genomes for GTDB-Tk
rule collect_all_genomes:
    input:
        good_bins_done = expand(f"{RESULTS_DIR}/06_good_bins/{{sample}}/done.txt", sample=FASTQ_SAMPLES)
    output:
        collected = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/collected.txt",
        summary = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/genome_list.txt"
    params:
        good_bins_dirs = expand(f"{RESULTS_DIR}/06_good_bins/{{sample}}", sample=FASTQ_SAMPLES),
        output_dir = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk",
        input_dir = "input"
    run:
        import os
        import shutil
        import glob
        
        # Create output directory
        os.makedirs(params.output_dir, exist_ok=True)
        
        genome_list = []
        
        # Collect all good bins from FASTQ samples
        for sample in FASTQ_SAMPLES:
            bins_dir = f"{RESULTS_DIR}/06_good_bins/{sample}"
            if os.path.exists(bins_dir):
                for bin_file in glob.glob(f"{bins_dir}/*.fa") + glob.glob(f"{bins_dir}/*.fasta"):
                    if os.path.exists(bin_file):
                        bin_name = os.path.basename(bin_file).replace('.fasta', '.fa')
                        # Name format: sample_binID.fa
                        new_name = f"{sample}_{bin_name}"
                        dst = os.path.join(params.output_dir, new_name)
                        shutil.copy(bin_file, dst)
                        genome_list.append(f"{new_name}\tMAG\t{sample}")
                        print(f"Copied {bin_name} from sample {sample}")
        
        # Copy FASTA samples
        for fasta in FASTA_SAMPLES:
            src = os.path.join(params.input_dir, fasta)
            if os.path.exists(src):
                # Keep original name but ensure .fa extension
                fasta_name = fasta if fasta.endswith('.fa') else fasta.replace('.fasta', '.fa')
                dst = os.path.join(params.output_dir, fasta_name)
                shutil.copy(src, dst)
                genome_list.append(f"{fasta_name}\tReference\tInput")
                print(f"Copied reference genome {fasta}")
        
        # Write summary
        with open(output.summary, 'w') as f:
            f.write("Genome\tType\tSource\n")
            for entry in genome_list:
                f.write(f"{entry}\n")
        
        # Write collected marker
        with open(output.collected, 'w') as f:
            f.write(f"Collected {len(genome_list)} genomes total\n")
            mag_count = sum(1 for g in genome_list if "MAG" in g)
            ref_count = sum(1 for g in genome_list if "Reference" in g)
            f.write(f"  MAGs: {mag_count}\n")
            f.write(f"  Reference genomes: {ref_count}\n")

# Rule 8: GTDB-Tk classification for all genomes
rule gtdbtk_classify:
    input:
        collected = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/collected.txt"
    output:
        summary = f"{RESULTS_DIR}/08_gtdbtk/gtdbtk.bac120.summary.tsv",
        done = f"{RESULTS_DIR}/08_gtdbtk/done.txt"
    params:
        genome_dir = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk",
        outdir = f"{RESULTS_DIR}/08_gtdbtk",
        gtdbtk_db = "/lustre/BIF/nobackup/mulle088/gapseq_pipeline_final/gtdbtk_db/release226",
        gtdbtk_container = f"{SINGULARITY_DIR}/quay.io-biocontainers-gtdbtk-2.4.1--pyhdfd78af_1.img"
    threads: 40
    shell:
        """
        export GTDBTK_DATA_PATH={params.gtdbtk_db}
        
        # Run GTDB-Tk classify workflow
        singularity exec -B {params.gtdbtk_db} {params.gtdbtk_container} \
            gtdbtk classify_wf \
            --genome_dir {params.genome_dir} \
            --extension fa \
            --out_dir {params.outdir} \
            --cpus {threads} \
            --skip_ani_screen \
            --pplacer_cpus 20
                    
        touch {output.done}
        """

# Create checkpoint to get genome list dynamically
checkpoint get_genome_list:
    input:
        collected = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/collected.txt"
    output:
        genome_list = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/genome_names.txt"
    shell:
        """
        ls {RESULTS_DIR}/07_all_genomes_for_gtdbtk/*.fa | xargs -n1 basename | sed 's/.fa$//' > {output.genome_list}
        """

# Rule 9: Gapseq find pathways - SIMPLIFIED
rule gapseq_pathways_single:
    input:
        genome = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/{{genome}}.fa"
    output:
        pathways = f"{RESULTS_DIR}/09_gapseq_pathways/{{genome}}-all-Pathways.tbl",
        reactions = f"{RESULTS_DIR}/09_gapseq_pathways/{{genome}}-all-Reactions.tbl"
    params:
        outdir = f"{RESULTS_DIR}/09_gapseq_pathways",
        gapseq_container = f"{SINGULARITY_DIR}/quay.io-biocontainers-gapseq-1.4.0--h9ee0642_1.img"
    threads: 4
    shell:
        """
        cd {params.outdir}
        singularity exec {params.gapseq_container} \
            gapseq find -p all ../../{input.genome}
        """

# Rule 10: Gapseq find transporters - SIMPLIFIED
rule gapseq_transporters_single:
    input:
        genome = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/{{genome}}.fa"
    output:
        transporters = f"{RESULTS_DIR}/10_gapseq_transporters/{{genome}}-Transporter.tbl"
    params:
        outdir = f"{RESULTS_DIR}/10_gapseq_transporters",
        gapseq_container = f"{SINGULARITY_DIR}/quay.io-biocontainers-gapseq-1.4.0--h9ee0642_1.img"
    threads: 4
    shell:
        """
        cd {params.outdir}
        singularity exec {params.gapseq_container} \
            gapseq find-transport ../../{input.genome}
        """

# Rule 11: Gapseq draft - SIMPLIFIED
rule gapseq_draft_single:
    input:
        genome = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/{{genome}}.fa",
        pathways = f"{RESULTS_DIR}/09_gapseq_pathways/{{genome}}-all-Pathways.tbl",
        reactions = f"{RESULTS_DIR}/09_gapseq_pathways/{{genome}}-all-Reactions.tbl",
        transporters = f"{RESULTS_DIR}/10_gapseq_transporters/{{genome}}-Transporter.tbl"
    output:
        draft = f"{RESULTS_DIR}/11_gapseq_draft/{{genome}}-draft.RDS"
    params:
        outdir = f"{RESULTS_DIR}/11_gapseq_draft",
        gapseq_container = f"{SINGULARITY_DIR}/quay.io-biocontainers-gapseq-1.4.0--h9ee0642_1.img"
    threads: 4
    shell:
        """
        cd {params.outdir}
        singularity exec {params.gapseq_container} \
            gapseq draft \
            -r ../../{input.reactions} \
            -t ../../{input.transporters} \
            -p ../../{input.pathways} \
            -c ../../{input.genome}
        """
# Aggregate function using checkpoint
def aggregate_gapseq_outputs(wildcards):
    checkpoint_output = checkpoints.get_genome_list.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f if line.strip()]
    
    return {
        'pathways': expand(f"{RESULTS_DIR}/09_gapseq_pathways/{{genome}}-all-Pathways.tbl", genome=genomes),
        'transporters': expand(f"{RESULTS_DIR}/10_gapseq_transporters/{{genome}}-Transporter.tbl", genome=genomes),
        'drafts': expand(f"{RESULTS_DIR}/11_gapseq_draft/{{genome}}-draft.RDS", genome=genomes)
    }


# Create a function to get all genome-media combinations
def get_genome_media_combinations():
    """
    Parse metadata to get all genome-media combinations
    Returns list of tuples: (genome, media)
    """
    metadata = pd.read_csv("input/meta_input_pipeline.csv")
    combinations = []
    
    for media in ["Taurine", "Creatinine", "Carnitine", "Xylan", "Chitin"]:
        if media in metadata.columns:
            samples = metadata[media].dropna().tolist()
            for sample in samples:
                sample_str = str(sample)
                
                # Handle FASTA samples
                if sample_str.endswith('.fa'):
                    combinations.append(("Marinacidobacteraceae", media))
                # Handle FASTQ samples (MAGs)
                else:
                    if '.' in sample_str:
                        sample_str = str(int(float(sample_str)))
                    
                    # Get all bins for this sample from the genome directory
                    genome_dir = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk"
                    if os.path.exists(genome_dir):
                        # Look for all genomes starting with this sample ID
                        for genome_file in glob.glob(f"{genome_dir}/{sample_str}_*.fa"):
                            genome_name = os.path.basename(genome_file).replace('.fa', '')
                            combinations.append((genome_name, media))
    
    return combinations

# Get all combinations
GENOME_MEDIA_COMBINATIONS = get_genome_media_combinations()

rule gapseq_fill_multi_media:
    input:
        draft = f"{RESULTS_DIR}/11_gapseq_draft/{{genome}}-draft.RDS",
        weights = f"{RESULTS_DIR}/11_gapseq_draft/{{genome}}-rxnWeights.RDS",
        genes = f"{RESULTS_DIR}/11_gapseq_draft/{{genome}}-rxnXgenes.RDS",
        media = "/lustre/BIF/nobackup/mulle088/snakemake/media/{media}_media.csv"
    output:
        model = f"{RESULTS_DIR}/12_gapseq_models/{{genome}}_{{media}}.RDS"
    params:
        outdir = f"{RESULTS_DIR}/12_gapseq_models",
        gapseq_container = f"{SINGULARITY_DIR}/quay.io-biocontainers-gapseq-1.4.0--h9ee0642_1.img"
    threads: 4
    shell:
        """
        cd {params.outdir}
        
        # Run gapseq fill (it will create {wildcards.genome}.RDS)
        singularity exec {params.gapseq_container} \
            gapseq fill \
            -m ../../{input.draft} \
            -c ../../{input.weights} \
            -g ../../{input.genes} \
            -n {input.media} \
            -o {wildcards.genome}_temp_{wildcards.media}
        
        # Rename the output to include media name
        if [ -f "{wildcards.genome}_temp_{wildcards.media}.RDS" ]; then
            mv {wildcards.genome}_temp_{wildcards.media}.RDS {wildcards.genome}_{wildcards.media}.RDS
        elif [ -f "{wildcards.genome}.RDS" ]; then
            mv {wildcards.genome}.RDS {wildcards.genome}_{wildcards.media}.RDS
        fi
        
        # Also handle the XML file if created
        if [ -f "{wildcards.genome}_temp_{wildcards.media}.xml" ]; then
            mv {wildcards.genome}_temp_{wildcards.media}.xml {wildcards.genome}_{wildcards.media}.xml
        elif [ -f "{wildcards.genome}.xml" ]; then
            mv {wildcards.genome}.xml {wildcards.genome}_{wildcards.media}.xml
        fi
        """

def get_all_genome_media_models(wildcards):
    """
    Dynamically determine all genome-media combinations that should exist
    """
    import pandas as pd
    import glob
    import os
    
    metadata = pd.read_csv("input/meta_input_pipeline.csv")
    expected_models = []
    
    # Get all existing genomes
    genome_dir = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk"
    if not os.path.exists(genome_dir):
        return []
    
    existing_genomes = [os.path.basename(f).replace('.fa', '') 
                       for f in glob.glob(f"{genome_dir}/*.fa")]
    
    # For each genome, determine which media it should have models for
    for genome in existing_genomes:
        genome_media = []
        
        if genome == "Marinacidobacteraceae":
            # FASTA sample - check all media columns
            for media in ["Taurine", "Creatinine", "Carnitine", "Xylan", "Chitin"]:
                if media in metadata.columns:
                    samples = metadata[media].dropna()
                    if any(str(s).endswith('.fa') for s in samples):
                        genome_media.append(media)
        else:
            # MAG sample - get sample ID and check media
            sample_id = genome.split('_')[0]
            
            for media in ["Taurine", "Creatinine", "Carnitine", "Xylan", "Chitin"]:
                if media in metadata.columns:
                    samples = metadata[media].dropna()
                    for s in samples:
                        s_str = str(s)
                        if '.' in s_str and not s_str.endswith('.fa'):
                            s_str = str(int(float(s_str)))
                        if s_str == sample_id:
                            genome_media.append(media)
                            break
        
        # Add expected models for this genome
        for media in genome_media:
            expected_models.append(f"{RESULTS_DIR}/12_gapseq_models/{genome}_{media}.RDS")
    
    return expected_models

# Update the aggregation rule to include all genome-media combinations
rule gapseq_all_models_done:
    input:
        models = [f"{RESULTS_DIR}/12_gapseq_models/{genome}_{media}.RDS" 
                  for genome, media in GENOME_MEDIA_COMBINATIONS]
    output:
        done = f"{RESULTS_DIR}/12_gapseq_models/all_models_complete.txt"
    shell:
        """
        echo "Created models for all genome-media combinations:" > {output.done}
        echo "Total models: {len(input.models)}" >> {output.done}
        
        # Count models per media type
        for media in Taurine Creatinine Carnitine Xylan Chitin; do
            count=$(ls {RESULTS_DIR}/12_gapseq_models/*_${{media}}.RDS 2>/dev/null | wc -l)
            echo "$media: $count models" >> {output.done}
        done
        """

rule gapseq_multi_media_complete:
    input:
        models = get_all_genome_media_models
    output:
        done = f"{RESULTS_DIR}/12_gapseq_models/multi_media_complete.txt",
        summary = f"{RESULTS_DIR}/12_gapseq_models/multi_media_summary.txt"
    run:
        import os
        
        # Count models by type
        media_counts = {"Taurine": 0, "Creatinine": 0, "Carnitine": 0, "Xylan": 0, "Chitin": 0}
        genome_counts = {}
        total_models = 0
        
        for model_path in input.models:
            if os.path.exists(model_path):
                total_models += 1
                # Extract genome and media from filename
                basename = os.path.basename(model_path).replace('.RDS', '')
                parts = basename.rsplit('_', 1)
                if len(parts) == 2:
                    genome, media = parts
                    
                    # Count by media
                    if media in media_counts:
                        media_counts[media] += 1
                    
                    # Count by genome
                    if genome not in genome_counts:
                        genome_counts[genome] = []
                    genome_counts[genome].append(media)
        
        # Write summary
        with open(output.summary, 'w') as f:
            f.write(f"Multi-Media Gap-Filling Summary\n")
            f.write(f"================================\n\n")
            f.write(f"Total models created: {total_models}\n")
            f.write(f"Total genomes: {len(genome_counts)}\n\n")
            
            f.write("Models per media type:\n")
            for media, count in sorted(media_counts.items()):
                f.write(f"  {media}: {count} models\n")
            
            f.write(f"\nGenomes with multiple media:\n")
            for genome, media_list in sorted(genome_counts.items()):
                if len(media_list) > 1:
                    f.write(f"  {genome}: {', '.join(media_list)}\n")
            
            f.write(f"\nAll genome-media combinations:\n")
            for genome, media_list in sorted(genome_counts.items()):
                for media in media_list:
                    f.write(f"  {genome}_{media}.RDS\n")
        
        # Write completion marker
        with open(output.done, 'w') as f:
            f.write(f"Completed {total_models} multi-media models\n")
            for model_path in sorted(input.models):
                if os.path.exists(model_path):
                    f.write(f"  ✓ {os.path.basename(model_path)}\n")


rule antismash_single:
    input:
        genome = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/{{genome}}.fa",
    output:
        done = f"{RESULTS_DIR}/13_antismash/{{genome}}/done.txt",
        json = f"{RESULTS_DIR}/13_antismash/{{genome}}/{{genome}}.json",
        gbk = f"{RESULTS_DIR}/13_antismash/{{genome}}/{{genome}}.gbk"
    params:
        outdir = f"{RESULTS_DIR}/13_antismash/{{genome}}",
        antismash_container = f"{SINGULARITY_DIR}/antismash-standalone-nonfree-6.0.0.img"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        # Run antiSMASH
        singularity exec {params.antismash_container} \
            antismash \
            --taxon bacteria \
            --genefinding-tool prodigal \
            --output-dir {params.outdir} \
            --output-basename {wildcards.genome} \
            --cpus {threads} \
            --minimal \
            --enable-genefunctions \
            --enable-lanthipeptides \
            --enable-lassopeptides \
            --enable-nrps-pks \
            --enable-sactipeptides \
            --enable-t2pks \
            --enable-thiopeptides \
            --enable-tta \
            {input.genome}
        touch {output.done}
        """

def aggregate_antismash_outputs(wildcards):
    checkpoint_output = checkpoints.get_genome_list.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f if line.strip()]
    return expand(f"{RESULTS_DIR}/13_antismash/{{genome}}/done.txt", genome=genomes)

rule antismash_summary:
    input:
        antismash_done = aggregate_antismash_outputs,
        # Ensure gapseq is complete before summarizing antiSMASH
        gapseq_done = f"{RESULTS_DIR}/12_gapseq_models/multi_media_complete.txt"
    output:
        summary = f"{RESULTS_DIR}/13_antismash/bgc_summary.txt",
        done = f"{RESULTS_DIR}/13_antismash/all_complete.txt"
    params:
        antismash_dir = f"{RESULTS_DIR}/13_antismash"
    run:
        import json
        import os
        
        bgc_summary = []
        genome_bgc_counts = {}
        bgc_types = {}
        
        # Parse each genome's antiSMASH results
        for done_file in input.antismash_done:
            genome_dir = os.path.dirname(done_file)
            genome_name = os.path.basename(genome_dir)
            json_file = os.path.join(genome_dir, f"{genome_name}.json")
            
            if os.path.exists(json_file):
                try:
                    with open(json_file, 'r') as f:
                        data = json.load(f)
                    
                    # Count BGCs for this genome
                    records = data.get('records', [])
                    total_bgcs = 0
                    genome_bgc_types = []
                    
                    for record in records:
                        areas = record.get('areas', [])
                        for area in areas:
                            total_bgcs += 1
                            products = area.get('products', [])
                            for product in products:
                                genome_bgc_types.append(product)
                                if product not in bgc_types:
                                    bgc_types[product] = 0
                                bgc_types[product] += 1
                    
                    genome_bgc_counts[genome_name] = {
                        'total': total_bgcs,
                        'types': genome_bgc_types
                    }
                    
                except Exception as e:
                    print(f"Error parsing {json_file}: {e}")
                    genome_bgc_counts[genome_name] = {'total': 0, 'types': []}
            else:
                genome_bgc_counts[genome_name] = {'total': 0, 'types': []}
        
        # Write summary
        with open(output.summary, 'w') as f:
            f.write("antiSMASH Biosynthetic Gene Cluster Summary\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Total genomes analyzed: {len(genome_bgc_counts)}\n")
            f.write(f"Total BGCs found: {sum(g['total'] for g in genome_bgc_counts.values())}\n\n")
            
            f.write("BGC types found:\n")
            for bgc_type, count in sorted(bgc_types.items(), key=lambda x: x[1], reverse=True):
                f.write(f"  {bgc_type}: {count}\n")
            
            f.write("\nBGCs per genome:\n")
            for genome, info in sorted(genome_bgc_counts.items(), key=lambda x: x[1]['total'], reverse=True):
                if info['total'] > 0:
                    f.write(f"  {genome}: {info['total']} BGCs\n")
                    type_counts = {}
                    for t in info['types']:
                        if t not in type_counts:
                            type_counts[t] = 0
                        type_counts[t] += 1
                    for t, c in sorted(type_counts.items()):
                        f.write(f"    - {t}: {c}\n")
                else:
                    f.write(f"  {genome}: No BGCs found\n")
        
        # Write completion marker
        with open(output.done, 'w') as f:
            f.write(f"antiSMASH analysis complete for {len(genome_bgc_counts)} genomes\n")


# Rule 14: Run dbCAN for CAZyme annotation
rule run_dbcan:
    input:
        genome = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/{{genome}}.fa"
    output:
        overview = f"{RESULTS_DIR}/14_cazymes/{{genome}}/overview.txt",
        done = f"{RESULTS_DIR}/14_cazymes/{{genome}}/done.txt"
    params:
        outdir = f"{RESULTS_DIR}/14_cazymes/{{genome}}",
        dbcan_container = f"{SINGULARITY_DIR}/run_dbcan_latest.sif",  # Updated container
        db_dir = "/lustre/BIF/nobackup/mulle088/databases/dbcan"
    threads: 8
    shell:
        """
        mkdir -p {params.outdir}
        
        # Run dbCAN with the latest container
        singularity exec {params.dbcan_container} \
            run_dbcan \
            {input.genome} \
            prok \
            --out_dir {params.outdir} \
            --db_dir {params.db_dir} \
            --use_signalP=FALSE \
            --dia_cpu {threads} \
            --hmm_cpu {threads}
        
        touch {output.done}
        """

# Aggregate CAZyme results for all genomes
def aggregate_cazyme_outputs(wildcards):
    checkpoint_output = checkpoints.get_genome_list.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f if line.strip()]
    return expand(f"{RESULTS_DIR}/14_cazymes/{{genome}}/done.txt", genome=genomes)

# Mark CAZyme analysis complete
rule cazymes_complete:
    input:
        aggregate_cazyme_outputs
    output:
        f"{RESULTS_DIR}/14_cazymes/all_done.txt"
    shell:
        """
        echo "CAZyme analysis complete" > {output}
        echo "Analyzed genomes:" >> {output}
        ls -d {RESULTS_DIR}/14_cazymes/*/ | wc -l >> {output}
        """

# Rule 15: Assess quality of MAGs (good bins) with QUAST
rule mag_assembly_stats:
    input:
        mag = f"{RESULTS_DIR}/07_all_genomes_for_gtdbtk/{{genome}}.fa"
    output:
        report = f"{RESULTS_DIR}/15_qc_summary/mag_stats/{{genome}}/report.tsv",
        stats = f"{RESULTS_DIR}/15_qc_summary/mag_stats/{{genome}}_stats.tsv"
    params:
        outdir = f"{RESULTS_DIR}/15_qc_summary/mag_stats/{{genome}}",
        quast_container = f"{SINGULARITY_DIR}/quay.io-biocontainers-quast-5.2.0--py39pl5321h2add14b_1.img"
    threads: 4
    shell:
        """
        # Run QUAST on MAG
        singularity exec {params.quast_container} \
            quast.py \
            {input.mag} \
            -o {params.outdir} \
            --threads {threads} \
            --min-contig 500 \
            --no-plots
        
        # Extract key MAG metrics
        echo -e "Genome\tTotal_contigs\tTotal_length\tN50\tN75\tN90\tLargest_contig\tGC_content" > {output.stats}
        
        awk -F'\t' '
            NR==1 {{for(i=1;i<=NF;i++) col[$i]=i}}
            /^# contigs \(>= 0 bp\)/ {{contigs=$2}}
            /^Total length \(>= 0 bp\)/ {{total=$2}}
            /^N50/ && !/N50[0-9]/ {{n50=$2}}
            /^N75/ {{n75=$2}}
            /^N90/ {{n90=$2}}
            /^Largest contig/ {{largest=$2}}
            /^GC \(%\)/ {{gc=$2}}
            END {{
                print "{wildcards.genome}", contigs, total, n50, n75, n90, largest, gc
            }}
        ' {output.report} | tr ' ' '\t' >> {output.stats}
        """

# Aggregate function for MAG stats
def aggregate_mag_stats(wildcards):
    checkpoint_output = checkpoints.get_genome_list.get(**wildcards).output[0]
    with open(checkpoint_output) as f:
        genomes = [line.strip() for line in f if line.strip()]
    return expand(f"{RESULTS_DIR}/15_qc_summary/mag_stats/{{genome}}_stats.tsv", genome=genomes)

# Rule 16: Comprehensive MAG quality report combining CheckM + QUAST + GTDB-Tk
rule comprehensive_mag_qc:
    input:
        mag_stats = aggregate_mag_stats,
        checkm_tables = expand(f"{RESULTS_DIR}/05_checkm/{{sample}}/checkm_table.tsv", sample=FASTQ_SAMPLES),
        gtdbtk = f"{RESULTS_DIR}/08_gtdbtk/gtdbtk.bac120.summary.tsv"
    output:
        master = f"{RESULTS_DIR}/15_qc_summary/mag_quality_master.tsv",
        summary = f"{RESULTS_DIR}/15_qc_summary/mag_quality_summary.txt",
        failed = f"{RESULTS_DIR}/15_qc_summary/failed_mags.tsv"
    run:
        import pandas as pd
        import os
        
        # Read all MAG stats
        mag_stats_dfs = []
        for stats_file in input.mag_stats:
            df = pd.read_csv(stats_file, sep='\t')
            mag_stats_dfs.append(df)
        mag_stats_df = pd.concat(mag_stats_dfs, ignore_index=True)
        
        # Read GTDB-Tk
        gtdbtk_df = pd.read_csv(input.gtdbtk, sep='\t')
        gtdbtk_dict = {row['user_genome']: {
            'classification': row.get('classification', 'Unknown'),
            'fastani_ani': row.get('fastani_ani', 'N/A'),
            'fastani_reference': row.get('fastani_reference', 'N/A')
        } for _, row in gtdbtk_df.iterrows()}
        
        # Compile CheckM data for all MAGs
        checkm_data = {}
        for sample in FASTQ_SAMPLES:
            checkm_file = f"{RESULTS_DIR}/05_checkm/{sample}/checkm_table.tsv"
            if os.path.exists(checkm_file):
                with open(checkm_file, 'r') as f:
                    for line in f.readlines()[1:]:
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 7:
                                bin_id = parts[0]
                                mag_name = f"{sample}_{bin_id}"
                                checkm_data[mag_name] = {
                                    'completeness': float(parts[5]),
                                    'contamination': float(parts[6]),
                                    'strain_heterogeneity': float(parts[7]) if len(parts) > 7 else 0
                                }
        
        # Add reference genome to CheckM data (100% complete, 0% contamination assumed)
        for fasta in FASTA_SAMPLES:
            genome_name = fasta.replace('.fa', '')
            checkm_data[genome_name] = {
                'completeness': 100.0,
                'contamination': 0.0,
                'strain_heterogeneity': 0.0
            }
        
        # Create master table
        master_data = []
        for _, row in mag_stats_df.iterrows():
            genome = row['Genome']
            
            # Get CheckM metrics
            checkm = checkm_data.get(genome, {})
            
            # Get GTDB-Tk classification
            gtdb = gtdbtk_dict.get(genome, {})
            
            # Quality tier determination
            completeness = checkm.get('completeness', 0)
            contamination = checkm.get('contamination', 100)
            
            if completeness >= 90 and contamination <= 5:
                quality_tier = "High"
            elif completeness >= 70 and contamination <= 10:
                quality_tier = "Medium"
            elif completeness >= 50 and contamination <= 10:
                quality_tier = "Low"
            else:
                quality_tier = "Failed"
            
            master_data.append({
                'Genome': genome,
                'Quality_Tier': quality_tier,
                'Completeness': completeness,
                'Contamination': contamination,
                'Strain_Het': checkm.get('strain_heterogeneity', 0),
                'Size_Mb': round(row['Total_length']/1000000, 2) if row['Total_length'] else 0,
                'Contigs': row['Total_contigs'],
                'N50': row['N50'],
                'GC%': row['GC_content'],
                'GTDB_Classification': gtdb.get('classification', 'Not classified'),
                'ANI%': gtdb.get('fastani_ani', 'N/A')
            })
        
        master_df = pd.DataFrame(master_data)
        master_df = master_df.sort_values(['Quality_Tier', 'Completeness'], ascending=[True, False])
        master_df.to_csv(output.master, sep='\t', index=False)
        
        # Failed MAGs
        failed_df = master_df[master_df['Quality_Tier'] == 'Failed']
        failed_df.to_csv(output.failed, sep='\t', index=False)
        
        # Summary report
        with open(output.summary, 'w') as f:
            f.write("MAG Quality Assessment Summary\n")
            f.write("=" * 60 + "\n\n")
            
            # Quality tiers
            tier_counts = master_df['Quality_Tier'].value_counts()
            f.write("Quality Tier Distribution:\n")
            for tier in ['High', 'Medium', 'Low', 'Failed']:
                count = tier_counts.get(tier, 0)
                pct = (count/len(master_df)*100) if len(master_df) > 0 else 0
                f.write(f"  {tier}: {count} ({pct:.1f}%)\n")
            
            f.write(f"\nTotal MAGs: {len(master_df)}\n")
            f.write(f"Average completeness: {master_df['Completeness'].mean():.1f}%\n")
            f.write(f"Average contamination: {master_df['Contamination'].mean():.1f}%\n")
            f.write(f"Average N50: {master_df['N50'].mean():.0f} bp\n")
            f.write(f"Average genome size: {master_df['Size_Mb'].mean():.2f} Mb\n")
            
            # Top MAGs
            f.write("\nTop 10 MAGs by completeness:\n")
            top_mags = master_df.nlargest(10, 'Completeness')[['Genome', 'Completeness', 'Contamination', 'Size_Mb', 'N50']]
            for _, row in top_mags.iterrows():
                f.write(f"  {row['Genome']}: {row['Completeness']:.1f}% complete, "
                       f"{row['Contamination']:.1f}% contam, {row['Size_Mb']:.2f} Mb, N50={row['N50']}\n")
            
            # Taxonomy distribution
            f.write("\nTaxonomic distribution (genus level):\n")
            master_df['Genus'] = master_df['GTDB_Classification'].apply(
                lambda x: x.split(';g__')[1].split(';')[0] if ';g__' in x else 'Unknown'
            )
            genus_counts = master_df[master_df['Genus'] != 'Unknown']['Genus'].value_counts().head(10)
            for genus, count in genus_counts.items():
                f.write(f"  {genus}: {count}\n")
