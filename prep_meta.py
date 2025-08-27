#!/usr/bin/env python3
"""
Enhanced Metadata Parser that handles samples with multiple media types
Creates separate pipeline entries for each sample-media combination
"""

import pandas as pd
import argparse
import os
from pathlib import Path

def parse_metadata_multi(metadata_file, input_dir="input"):
    """
    Parse metadata file allowing samples to have multiple media types
    Each sample-media combination becomes a separate pipeline entry
    """
    # Read the metadata file
    df = pd.read_csv(metadata_file)
    
    # Media types from columns
    media_types = ["Taurine", "Creatinine", "Carnitine", "Xylan", "Chitin"]
    
    # Initialize the mapping dictionaries
    sample_entries = []  # List of all sample-media combinations
    
    # Process each media type column
    for media in media_types:
        if media in df.columns:
            # Get all non-null values from this column
            samples = df[media].dropna().tolist()
            
            for sample in samples:
                sample_str = str(sample)
                
                # Clean up float representations (e.g., 459.0 → 459)
                if '.' in sample_str and sample_str.replace('.', '').replace('0', '').isdigit():
                    sample_str = str(int(float(sample)))
                
                # Check if it's a FASTA file
                if sample_str.endswith('.fa') or sample_str.endswith('.fasta'):
                    # For FASTA files, create unique ID with media suffix
                    sample_id = f"{sample_str}_{media}"
                    
                    # Look for the actual file
                    possible_paths = [
                        f"{input_dir}/{sample_str}",
                        f"{input_dir}/Marinacidobacteraceae_bin.13_nanopore_reassembled.fa",  # Check actual file
                        f"../gapseq_pipeline_final/{sample_str}",
                    ]
                    
                    actual_path = None
                    for path in possible_paths:
                        if os.path.exists(path):
                            actual_path = os.path.abspath(path)
                            break
                    
                    if actual_path:
                        sample_entries.append({
                            'sample_id': sample_id,
                            'original_name': sample_str,
                            'sample_type': 'fasta',
                            'media': media,
                            'path': actual_path
                        })
                    else:
                        print(f"Warning: FASTA file not found: {sample_str}")
                        
                # Check if it's a sample number (fastq)
                elif sample_str.isdigit():
                    # For FASTQ, create unique ID with media suffix if multiple media
                    sample_id = f"{sample_str}_{media}"
                    
                    sample_entries.append({
                        'sample_id': sample_id,
                        'original_name': sample_str,
                        'sample_type': 'fastq',
                        'media': media,
                        'path': 'fastq'
                    })
    
    return sample_entries

def create_sample_sheet_multi(sample_entries, data_dir):
    """
    Create sample sheet with unique entries for each sample-media combination
    """
    rows = []
    
    for entry in sample_entries:
        if entry['sample_type'] == 'fastq':
            # Check for FASTQ files
            r1_file = f"{data_dir}/{entry['original_name']}_R1_001.fastq.gz"
            r2_file = f"{data_dir}/{entry['original_name']}_R2_001.fastq.gz"
            
            if os.path.exists(r1_file) and os.path.exists(r2_file):
                rows.append({
                    'sample_id': entry['sample_id'],
                    'original_id': entry['original_name'],
                    'sample_type': 'fastq',
                    'media': entry['media'],
                    'r1': r1_file,
                    'r2': r2_file,
                    'fasta_path': 'NA'
                })
            else:
                print(f"Warning: FASTQ files not found for {entry['original_name']}")
        else:
            # FASTA entry
            rows.append({
                'sample_id': entry['sample_id'],
                'original_id': entry['original_name'],
                'sample_type': 'fasta',
                'media': entry['media'],
                'r1': 'NA',
                'r2': 'NA',
                'fasta_path': entry['path']
            })
    
    return pd.DataFrame(rows)

def print_summary_multi(df):
    """
    Print summary of the sample sheet with multiple media handling
    """
    print("\n" + "="*60)
    print("SAMPLE SHEET SUMMARY (Multi-Media Support)")
    print("="*60)
    
    # Count unique samples
    unique_samples = df['original_id'].nunique()
    total_entries = len(df)
    
    print(f"\nTotal unique samples: {unique_samples}")
    print(f"Total pipeline entries: {total_entries}")
    
    if total_entries > unique_samples:
        print(f"→ {total_entries - unique_samples} samples will be processed with multiple media")
    
    # Show samples with multiple media
    multi_media = df.groupby('original_id')['media'].apply(list).reset_index()
    multi_media = multi_media[multi_media['media'].apply(len) > 1]
    
    if not multi_media.empty:
        print("\nSamples with multiple media conditions:")
        for _, row in multi_media.iterrows():
            print(f"  {row['original_id']}: {', '.join(row['media'])}")
    
    # Count by type
    print(f"\nBy sample type:")
    print(f"  FASTQ entries: {len(df[df['sample_type'] == 'fastq'])}")
    print(f"  FASTA entries: {len(df[df['sample_type'] == 'fasta'])}")
    
    # Count by media
    print(f"\nBy media type:")
    for media in df['media'].unique():
        count = len(df[df['media'] == media])
        print(f"  {media}: {count} entries")
    
    print("\nSample details:")
    for _, row in df.iterrows():
        if row['original_id'] != row['sample_id']:
            print(f"  {row['sample_id']} (from {row['original_id']}, {row['sample_type']}): {row['media']} media")
        else:
            print(f"  {row['sample_id']} ({row['sample_type']}): {row['media']} media")

def main():
    parser = argparse.ArgumentParser(
        description='Prepare metadata with support for multiple media per sample'
    )
    parser.add_argument(
        'metadata',
        help='Path to the metadata CSV file'
    )
    parser.add_argument(
        '--data-dir',
        default='/lustre/BIF/nobackup/mulle088/data_rp/90-1196727703/00_fastq',
        help='Directory containing FASTQ files'
    )
    parser.add_argument(
        '--input-dir',
        default='input',
        help='Directory containing local FASTA files'
    )
    parser.add_argument(
        '--output',
        default='sample_sheet_multi.csv',
        help='Output sample sheet file'
    )
    
    args = parser.parse_args()
    
    # Parse metadata
    print(f"Parsing metadata from: {args.metadata}")
    print(f"Looking for FASTA files in: {args.input_dir}")
    print(f"Looking for FASTQ files in: {args.data_dir}")
    print()
    
    sample_entries = parse_metadata_multi(args.metadata, args.input_dir)
    
    if not sample_entries:
        print("No samples found in metadata!")
        return
    
    # Create sample sheet
    df = create_sample_sheet_multi(sample_entries, args.data_dir)
    
    if df.empty:
        print("No valid samples found!")
        return
    
    # Save sample sheet
    df.to_csv(args.output, index=False)
    print(f"\nSample sheet saved to: {args.output}")
    
    # Create config for FASTA paths
    fasta_entries = df[df['sample_type'] == 'fasta']
    if not fasta_entries.empty:
        with open('config_paths_multi.yaml', 'w') as f:
            f.write("# FASTA file paths for pipeline (multi-media version)\n")
            f.write("fasta_paths:\n")
            for _, row in fasta_entries.iterrows():
                f.write(f"  {row['sample_id']}: \"{row['fasta_path']}\"\n")
        print(f"FASTA paths config saved to: config_paths_multi.yaml")
    
    # Print summary
    print_summary_multi(df)
    
    print("\n" + "="*60)
    print("✓ Sample sheet created with multi-media support!")
    print("Each sample-media combination will be processed separately.")
    print("Use this with the pipeline for complete analysis.")

if __name__ == "__main__":
    main()
