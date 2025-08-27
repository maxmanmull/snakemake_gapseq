#!/usr/bin/env python3
"""
Results Summary Script for Metabolic Pipeline
Generates a comprehensive summary of all pipeline results
"""

import os
import sys
import glob
import pandas as pd
import argparse
from pathlib import Path
import json

def parse_checkm_results(checkm_file):
    """Parse CheckM results file"""
    bins_info = []
    
    if not os.path.exists(checkm_file):
        return pd.DataFrame()
    
    with open(checkm_file, 'r') as f:
        lines = f.readlines()
        
    # Find the data lines (skip headers)
    for line in lines:
        if line.strip() and not line.startswith('-') and 'Completeness' not in line:
            parts = line.split()
            if len(parts) >= 13:
                try:
                    bin_info = {
                        'bin_id': parts[0],
                        'completeness': float(parts[11]),
                        'contamination': float(parts[12])
                    }
                    bins_info.append(bin_info)
                except (ValueError, IndexError):
                    continue
    
    return pd.DataFrame(bins_info)

def parse_gtdbtk_summary(gtdbtk_file):
    """Parse GTDB-Tk summary file"""
    if not os.path.exists(gtdbtk_file):
        return None
    
    try:
        df = pd.read_csv(gtdbtk_file, sep='\t')
        if not df.empty:
            # Get the most relevant columns
            return {
                'classification': df.iloc[0].get('classification', 'Unknown'),
                'fastani_reference': df.iloc[0].get('fastani_reference', 'N/A'),
                'fastani_ani': df.iloc[0].get('fastani_ani', 'N/A')
            }
    except Exception as e:
        print(f"Error parsing GTDB-Tk file: {e}")
    
    return None

def count_gapseq_results(gapseq_dir, bin_id):
    """Count Gapseq results for a bin"""
    results = {
        'pathways': 0,
        'transporters': 0,
        'has_draft': False,
        'has_model': False
    }
    
    # Check for pathways
    pathways_file = f"{gapseq_dir}/{bin_id}-Pathways.tbl"
    if os.path.exists(pathways_file):
        try:
            df = pd.read_csv(pathways_file, sep='\t')
            results['pathways'] = len(df)
        except:
            pass
    
    # Check for transporters
    transport_file = f"{gapseq_dir}/{bin_id}-Transporter.tbl"
    if os.path.exists(transport_file):
        try:
            df = pd.read_csv(transport_file, sep='\t')
            results['transporters'] = len(df)
        except:
            pass
    
    # Check for models
    results['has_draft'] = os.path.exists(f"{gapseq_dir}/{bin_id}_draft.RDS")
    results['has_model'] = os.path.exists(f"{gapseq_dir}/{bin_id}_model.RDS")
    
    return results

def summarize_fastq_sample(sample_id, results_dir):
    """Summarize results for a FASTQ sample"""
    summary = {
        'sample_id': sample_id,
        'sample_type': 'FASTQ',
        'assembly_completed': False,
        'total_bins': 0,
        'quality_bins': 0,
        'processed_bins': []
    }
    
    # Check assembly
    assembly_file = f"{results_dir}/02_assembly/{sample_id}/contigs.fasta"
    summary['assembly_completed'] = os.path.exists(assembly_file)
    
    # Check bins
    bins_dir = f"{results_dir}/04_binning/{sample_id}/bins"
    if os.path.exists(bins_dir):
        all_bins = glob.glob(f"{bins_dir}/*.fa")
        summary['total_bins'] = len(all_bins)
    
    # Parse CheckM results
    checkm_file = f"{results_dir}/06_checkm/{sample_id}/checkm_results.txt"
    checkm_df = parse_checkm_results(checkm_file)
    
    # Filter quality bins
    if not checkm_df.empty:
        quality_bins = checkm_df[
            (checkm_df['completeness'] >= 50) & 
            (checkm_df['contamination'] <= 10)
        ]
        summary['quality_bins'] = len(quality_bins)
        
        # Check processing for each quality bin
        for _, bin_row in quality_bins.iterrows():
            bin_id = bin_row['bin_id']
            bin_info = {
                'bin_id': bin_id,
                'completeness': bin_row['completeness'],
                'contamination': bin_row['contamination'],
                'taxonomy': None,
                'genes_predicted': False,
                'pathways': 0,
                'transporters': 0,
                'model_created': False
            }
            
            # Check GTDB-Tk
            gtdbtk_file = f"{results_dir}/07_gtdbtk/bins/{sample_id}/{bin_id}/gtdbtk.bac120.summary.tsv"
            taxonomy = parse_gtdbtk_summary(gtdbtk_file)
            if taxonomy:
                bin_info['taxonomy'] = taxonomy['classification']
            
            # Check Prodigal
            proteins_file = f"{results_dir}/08_prodigal/bins/{sample_id}/{bin_id}/proteins.faa"
            bin_info['genes_predicted'] = os.path.exists(proteins_file)
            
            # Check Gapseq
            gapseq_dir = f"{results_dir}/09_gapseq/bins/{sample_id}/{bin_id}"
            gapseq_results = count_gapseq_results(gapseq_dir, bin_id)
            bin_info.update(gapseq_results)
            bin_info['model_created'] = gapseq_results['has_model']
            
            summary['processed_bins'].append(bin_info)
    
    return summary

def summarize_fasta_sample(sample_id, results_dir):
    """Summarize results for a FASTA sample"""
    summary = {
        'sample_id': sample_id,
        'sample_type': 'FASTA',
        'taxonomy': None,
        'genes_predicted': False,
        'pathways': 0,
        'transporters': 0,
        'model_created': False
    }
    
    # Check GTDB-Tk
    gtdbtk_file = f"{results_dir}/07_gtdbtk/fasta/{sample_id}/gtdbtk.bac120.summary.tsv"
    taxonomy = parse_gtdbtk_summary(gtdbtk_file)
    if taxonomy:
        summary['taxonomy'] = taxonomy['classification']
    
    # Check Prodigal
    proteins_file = f"{results_dir}/08_prodigal/fasta/{sample_id}/proteins.faa"
    summary['genes_predicted'] = os.path.exists(proteins_file)
    
    # Check Gapseq
    gapseq_dir = f"{results_dir}/09_gapseq/fasta/{sample_id}"
    gapseq_results = count_gapseq_results(gapseq_dir, sample_id)
    summary['pathways'] = gapseq_results['pathways']
    summary['transporters'] = gapseq_results['transporters']
    summary['model_created'] = gapseq_results['has_model']
    
    return summary

def generate_report(summaries, output_file):
    """Generate a formatted report from summaries"""
    
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("METABOLIC PIPELINE RESULTS SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        # Separate by type
        fastq_samples = [s for s in summaries if s['sample_type'] == 'FASTQ']
        fasta_samples = [s for s in summaries if s['sample_type'] == 'FASTA']
        
        # FASTQ samples summary
        if fastq_samples:
            f.write("FASTQ SAMPLES (Raw Reads → MAGs)\n")
            f.write("-"*40 + "\n\n")
            
            for sample in fastq_samples:
                f.write(f"Sample: {sample['sample_id']}\n")
                f.write(f"  Assembly: {'✓' if sample['assembly_completed'] else '✗'}\n")
                f.write(f"  Total bins: {sample['total_bins']}\n")
                f.write(f"  Quality bins (>50% complete, <10% contaminated): {sample['quality_bins']}\n")
                
                if sample['processed_bins']:
                    f.write(f"  Processed MAGs:\n")
                    for bin_info in sample['processed_bins']:
                        f.write(f"    {bin_info['bin_id']}:\n")
                        f.write(f"      Completeness: {bin_info['completeness']:.1f}%\n")
                        f.write(f"      Contamination: {bin_info['contamination']:.1f}%\n")
                        if bin_info['taxonomy']:
                            f.write(f"      Taxonomy: {bin_info['taxonomy']}\n")
                        f.write(f"      Genes predicted: {'✓' if bin_info['genes_predicted'] else '✗'}\n")
                        f.write(f"      Pathways found: {bin_info['pathways']}\n")
                        f.write(f"      Transporters found: {bin_info['transporters']}\n")
                        f.write(f"      Model created: {'✓' if bin_info['model_created'] else '✗'}\n")
                f.write("\n")
        
        # FASTA samples summary
        if fasta_samples:
            f.write("\nFASTA SAMPLES (Pre-assembled Genomes)\n")
            f.write("-"*40 + "\n\n")
            
            for sample in fasta_samples:
                f.write(f"Sample: {sample['sample_id']}\n")
                if sample['taxonomy']:
                    f.write(f"  Taxonomy: {sample['taxonomy']}\n")
                f.write(f"  Genes predicted: {'✓' if sample['genes_predicted'] else '✗'}\n")
                f.write(f"  Pathways found: {sample['pathways']}\n")
                f.write(f"  Transporters found: {sample['transporters']}\n")
                f.write(f"  Model created: {'✓' if sample['model_created'] else '✗'}\n\n")
        
        # Statistics
        f.write("\n" + "="*80 + "\n")
        f.write("OVERALL STATISTICS\n")
        f.write("="*80 + "\n\n")
        
        total_fastq = len(fastq_samples)
        total_fasta = len(fasta_samples)
        total_mags = sum(s['quality_bins'] for s in fastq_samples)
        successful_mags = sum(len([b for b in s['processed_bins'] if b['model_created']]) 
                              for s in fastq_samples)
        successful_fasta = sum(1 for s in fasta_samples if s['model_created'])
        
        f.write(f"Total FASTQ samples processed: {total_fastq}\n")
        f.write(f"Total FASTA samples processed: {total_fasta}\n")
        f.write(f"Total high-quality MAGs generated: {total_mags}\n")
        f.write(f"MAGs with completed models: {successful_mags}\n")
        f.write(f"FASTA genomes with completed models: {successful_fasta}\n")
        f.write(f"Total metabolic models created: {successful_mags + successful_fasta}\n")

def main():
    parser = argparse.ArgumentParser(
        description='Summarize results from the metabolic pipeline'
    )
    parser.add_argument(
        '--results-dir',
        default='results',
        help='Path to the results directory'
    )
    parser.add_argument(
        '--metadata',
        default='../gapseq_pipeline_final/metadata.csv',
        help='Path to metadata file'
    )
    parser.add_argument(
        '--output',
        default='pipeline_summary.txt',
        help='Output summary file'
    )
    parser.add_argument(
        '--json',
        action='store_true',
        help='Also output JSON format'
    )
    
    args = parser.parse_args()
    
    if not os.path.exists(args.results_dir):
        print(f"Error: Results directory {args.results_dir} not found")
        sys.exit(1)
    
    # Parse metadata to get sample list
    metadata_df = pd.read_csv(args.metadata)
    media_types = ["Taurine", "Creatinine", "Carnitine", "Xylan", "Chitin"]
    
    all_samples = []
    for media in media_types:
        if media in metadata_df.columns:
            samples = metadata_df[media].dropna().tolist()
            for sample in samples:
                sample_str = str(sample)
                if sample_str.endswith('.fa') or sample_str.endswith('.fasta'):
                    all_samples.append(('FASTA', sample_str))
                elif sample_str.replace('.', '').isdigit() or sample_str.isdigit():
                    if '.' in sample_str:
                        sample_id = str(int(float(sample_str)))
                    else:
                        sample_id = sample_str
                    all_samples.append(('FASTQ', sample_id))
    
    # Process each sample
    print(f"Processing {len(all_samples)} samples...")
    summaries = []
    
    for sample_type, sample_id in all_samples:
        print(f"  Processing {sample_type} sample: {sample_id}")
        
        if sample_type == 'FASTQ':
            summary = summarize_fastq_sample(sample_id, args.results_dir)
        else:
            summary = summarize_fasta_sample(sample_id, args.results_dir)
        
        summaries.append(summary)
    
    # Generate report
    print(f"\nGenerating report: {args.output}")
    generate_report(summaries, args.output)
    
    # Optionally save JSON
    if args.json:
        json_output = args.output.replace('.txt', '.json')
        with open(json_output, 'w') as f:
            json.dump(summaries, f, indent=2)
        print(f"JSON summary saved to: {json_output}")
    
    print(f"\nSummary complete! View the report: {args.output}")

if __name__ == "__main__":
    main()
