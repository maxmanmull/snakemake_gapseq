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
        f"{RESULTS_DIR}/12_gapseq_models/multi_media_complete.txt"

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

