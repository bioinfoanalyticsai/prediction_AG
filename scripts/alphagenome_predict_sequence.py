#!/usr/bin/env python3
"""
AlphaGenome sequence prediction - predict signals for entire sequences.
Uses dna_model.predict_sequence to generate position-level predictions.
"""

import argparse
import pandas as pd
import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client
import warnings
import sys

warnings.filterwarnings('ignore')

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Predict AlphaGenome signals for DNA sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:

  # Predict for a direct DNA sequence
  python alphagenome_predict_sequence.py \
    --sequence "GATTACAGATTACA" \
    --api-key YOUR_KEY --tissues brain

  # Predict for genomic region (using Interval)
  python alphagenome_predict_sequence.py \
    --chrom chr6 --start 49163874 --end 49165042 \
    --api-key YOUR_KEY --tissues intestine

  # Predict for genomic region from FASTA file
  python alphagenome_predict_sequence.py \
    --chrom chr6 --start 49163874 --end 49165042 \
    --fasta mm10.fa.gz \
    --api-key YOUR_KEY --tissues brain

  # Predict for multiple regions from BED
  python alphagenome_predict_sequence.py \
    --bed-file regions.bed \
    --api-key YOUR_KEY --tissues intestine --output sequence_predictions/
        """
    )
    parser.add_argument('--sequence', type=str, help='Direct DNA sequence (e.g., "GATTACA")')
    parser.add_argument('--chrom', type=str, help='Chromosome (e.g., chr6)')
    parser.add_argument('--start', type=int, help='Start position (0-based)')
    parser.add_argument('--end', type=int, help='End position')
    parser.add_argument('--fasta', type=str, help='FASTA file for fetching sequences')
    parser.add_argument('--bed-file', type=str, help='BED file with regions')
    parser.add_argument('--api-key', type=str, required=True, help='AlphaGenome API key')
    parser.add_argument('--tissues', type=str, nargs='+', default=['brain'],
                       help='Tissues to predict')
    parser.add_argument('--output', type=str, nargs='+',
                       choices=['rna_seq', 'atac', 'dnase', 'cage', 'chip_tf', 
                               'chip_histone', 'procap', 'all'],
                       default=['all'],
                       help='Output types to save')
    parser.add_argument('--output-dir', type=str, default='./sequence_predictions/',
                       help='Output directory for results')
    parser.add_argument('--save-csv', action='store_true',
                       help='Save predictions as CSV (position x tracks)')
    parser.add_argument('--save-npy', action='store_true',
                       help='Save predictions as numpy arrays')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')
    
    return parser.parse_args()

TISSUE_MAPPING = {
    'brain': 'UBERON:0000955',
    'lung': 'UBERON:0002048',
    'colon': 'UBERON:0001155',
    'intestine': 'UBERON:0001157',
    'kidney': 'UBERON:0002113',
}

MM10_CHROM_SIZES = {
    'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
    'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
    'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
    'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
    'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299,
    'chrY': 91744698, 'chrM': 16571,
}

def tissue_to_uberon(tissues):
    """Convert tissue names to UBERON codes."""
    uberon_terms = []
    for tissue in tissues:
        if tissue.lower() in TISSUE_MAPPING:
            uberon_terms.append(TISSUE_MAPPING[tissue.lower()])
        else:
            print(f"⚠ Unknown tissue '{tissue}'. Supported: {', '.join(TISSUE_MAPPING.keys())}")
    return uberon_terms

def get_output_types(output_arg):
    """Convert output argument to AlphaGenome OutputType."""
    if 'all' in output_arg:
        return [
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.ATAC,
            dna_client.OutputType.DNASE,
            dna_client.OutputType.CAGE,
            dna_client.OutputType.CHIP_HISTONE,
            dna_client.OutputType.CHIP_TF,
            dna_client.OutputType.PROCAP,
        ]
    
    output_map = {
        'rna_seq': dna_client.OutputType.RNA_SEQ,
        'atac': dna_client.OutputType.ATAC,
        'dnase': dna_client.OutputType.DNASE,
        'cage': dna_client.OutputType.CAGE,
        'chip_tf': dna_client.OutputType.CHIP_TF,
        'chip_histone': dna_client.OutputType.CHIP_HISTONE,
        'procap': dna_client.OutputType.PROCAP,
    }
    
    return [output_map[o] for o in output_arg if o in output_map]

def fetch_sequence_from_fasta(fasta_file, chrom, start, end):
    """
    Fetch DNA sequence from FASTA file.
    Returns: DNA sequence string
    """
    try:
        from pysam import FastaFile
        with FastaFile(fasta_file) as f:
            sequence = f.fetch(chrom, start, end).upper()
        return sequence
    except ImportError:
        print("⚠ pysam not available. Use --sequence directly or install: pip install pysam")
        return None
    except Exception as e:
        print(f"✗ Error fetching sequence: {e}")
        return None

def predict_sequence_dna(dna_model, sequence, uberon_terms, output_types, verbose=False):
    """
    Predict signals for a DNA sequence string using predict_sequence.
    
    The sequence will be padded with 'N's to valid lengths.
    Valid lengths: 16KB, 100KB, 500KB, 1MB
    
    Returns: Output object with track data for each output type
    """
    if verbose:
        print(f"Input sequence length: {len(sequence)} bp")
    
    try:
        # Determine appropriate sequence length
        valid_lengths = {
            'SEQUENCE_LENGTH_16KB': dna_client.SEQUENCE_LENGTH_16KB,
            'SEQUENCE_LENGTH_100KB': dna_client.SEQUENCE_LENGTH_100KB,
            'SEQUENCE_LENGTH_500KB': dna_client.SEQUENCE_LENGTH_500KB,
            'SEQUENCE_LENGTH_1MB': dna_client.SEQUENCE_LENGTH_1MB,
        }
        
        # Use smallest valid length that fits the sequence
        seq_length = dna_client.SEQUENCE_LENGTH_1MB
        for length_name in ['SEQUENCE_LENGTH_16KB', 'SEQUENCE_LENGTH_100KB', 'SEQUENCE_LENGTH_500KB', 'SEQUENCE_LENGTH_1MB']:
            length_val = valid_lengths[length_name]
            if len(sequence) <= length_val:
                seq_length = length_val
                if verbose:
                    print(f"Using sequence length: {length_name}")
                break
        
        # Pad sequence to valid length
        padded_sequence = sequence.center(seq_length, 'N')
        
        if verbose:
            print(f"Padded sequence length: {len(padded_sequence)} bp")
            print(f"Predicting signals...")
        
        # Use predict_sequence with DNA string
        output = dna_model.predict_sequence(
            sequence=padded_sequence,
            requested_outputs=output_types,
            ontology_terms=uberon_terms,
        )
        
        if verbose:
            print(f"✓ Prediction complete")
        
        return output
    
    except Exception as e:
        print(f"✗ Error predicting sequence: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return None

def predict_sequence_region(dna_model, chrom, start, end, uberon_terms, output_types, fasta_file=None, verbose=False):
    """
    Predict signals for a genomic region.
    
    Can use either:
    1. Fetch from FASTA file (if provided)
    2. Use Interval directly (fallback)
    
    Returns: Output object with track data for each output type
    """
    sequence = None
    
    # Try to fetch from FASTA first
    if fasta_file:
        if verbose:
            print(f"Fetching sequence from FASTA: {chrom}:{start}-{end}")
        sequence = fetch_sequence_from_fasta(fasta_file, chrom, start, end)
    
    if sequence:
        # Use DNA sequence method
        if verbose:
            print(f"Using predict_sequence with DNA string")
        return predict_sequence_dna(dna_model, sequence, uberon_terms, output_types, verbose)
    
    else:
        # Fallback to predict_interval
        if verbose:
            print(f"Using predict_interval (interval-based)")
        
        try:
            interval = genome.Interval(
                chromosome=chrom,
                start=start,
                end=end
            )
            
            if verbose:
                print(f"Predicting on interval ({end - start} bp)...")
            
            output = dna_model.predict_interval(
                interval=interval,
                requested_outputs=output_types,
                ontology_terms=uberon_terms,
                organism=dna_client.Organism.MUS_MUSCULUS
            )
            
            if verbose:
                print(f"✓ Prediction complete")
            
            return output
        
        except Exception as e:
            print(f"✗ Error predicting sequence: {e}")
            if verbose:
                import traceback
                traceback.print_exc()
            return None

def extract_sequence_predictions(output, output_types):
    """
    Extract predictions from Output object.
    Returns: Dictionary mapping track types to prediction arrays and metadata.
    
    Note: predict_interval returns aggregated values (mean, max, etc.) per track,
    not position-level data. The script saves these summary values.
    For true position-level predictions, position_index represents regions,
    not individual bases.
    """
    predictions = {}
    
    try:
        tracks_map = {
            'rna_seq': output.rna_seq,
            'atac': output.atac,
            'dnase': output.dnase,
            'cage': output.cage,
            'chip_histone': output.chip_histone,
            'chip_tf': output.chip_tf,
            'procap': output.procap,
        }
        
        for track_name, track_data in tracks_map.items():
            if track_data is None:
                continue
            
            if hasattr(track_data, 'values') and track_data.values is not None:
                values = track_data.values
                
                # values shape: (n_positions, n_tracks) or similar
                if values.size > 0:
                    predictions[track_name] = {
                        'values': values,  # Prediction array
                        'shape': values.shape,
                    }
                    
                    if hasattr(track_data, 'metadata') and track_data.metadata is not None:
                        predictions[track_name]['metadata'] = track_data.metadata
                        predictions[track_name]['n_tracks'] = len(track_data.metadata)
                    else:
                        if len(values.shape) == 2:
                            predictions[track_name]['n_tracks'] = values.shape[1]
                        else:
                            predictions[track_name]['n_tracks'] = 1
        
        return predictions
    
    except Exception as e:
        print(f"✗ Error extracting predictions: {e}")
        return {}

def save_sequence_predictions(predictions, chrom, start, end, output_dir, save_csv=False, save_npy=False, verbose=False):
    """Save sequence predictions to files."""
    import os
    from pathlib import Path
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    region_name = f"{chrom}_{start}_{end}"
    
    for track_name, track_data in predictions.items():
        values = track_data['values']
        n_positions, n_tracks = values.shape
        
        if verbose:
            print(f"\n  {track_name}:")
            print(f"    Shape: {values.shape} ({n_positions} positions x {n_tracks} tracks)")
        
        # Save as NPY (numpy binary)
        if save_npy:
            npy_path = output_path / f"{region_name}_{track_name}.npy"
            np.save(npy_path, values)
            if verbose:
                print(f"    Saved: {npy_path}")
        
        # Save as CSV
        if save_csv:
            track_labels = []
            if 'metadata' in track_data:
                track_labels = [f"{name}" for name in track_data['metadata']['name'].tolist()]
            else:
                track_labels = [f"track_{i}" for i in range(n_tracks)]
            
            df = pd.DataFrame(values, columns=track_labels)
            df.insert(0, 'position', np.arange(start, end))
            
            csv_path = output_path / f"{region_name}_{track_name}.csv"
            df.to_csv(csv_path, index=False)
            if verbose:
                print(f"    Saved: {csv_path}")
    
    return output_path

def print_sequence_summary(predictions, region_name, region_size):
    """Print summary of sequence predictions."""
    print("\n" + "="*80)
    print(f"Sequence Prediction Summary: {region_name}")
    print("="*80 + "\n")
    
    print(f"Region size: {region_size:,} bp")
    print(f"\nPrediction Tracks:\n")
    
    for track_name, track_data in predictions.items():
        values = track_data['values']
        n_positions, n_tracks = values.shape
        
        print(f"  {track_name}:")
        print(f"    Positions: {n_positions:,}")
        print(f"    Tracks (tissues/cell types): {n_tracks}")
        print(f"    Data shape: {values.shape}")
        print(f"    Value range: [{np.min(values):.6f}, {np.max(values):.6f}]")
        print(f"    Mean value: {np.mean(values):.6f}")
        print(f"    Std dev: {np.std(values):.6f}")
        
        if 'metadata' in track_data:
            track_names = track_data['metadata']['name'].tolist()[:3]
            print(f"    Sample tracks: {', '.join(str(t) for t in track_names)}")
            if n_tracks > 3:
                print(f"               ... and {n_tracks - 3} more")
        
        print()

def main():
    args = parse_arguments()
    
    print("\n" + "="*80)
    print("AlphaGenome Sequence Prediction")
    print("="*80 + "\n")
    
    # Get regions/sequences to predict
    regions = []
    
    if args.sequence:
        # Direct DNA sequence
        print(f"✓ Using direct DNA sequence: {args.sequence[:50]}{'...' if len(args.sequence) > 50 else ''}")
        regions.append(('sequence', args.sequence, None))
    
    elif args.bed_file:
        # Load from BED
        try:
            bed_df = pd.read_csv(args.bed_file, sep='\t', header=None, comment='#')
            if bed_df.shape[1] < 3:
                raise ValueError("BED file must have at least 3 columns")
            bed_df.columns = ['chrom', 'start', 'end'] + [f'col_{i}' for i in range(3, bed_df.shape[1])]
            bed_df['start'] = pd.to_numeric(bed_df['start'])
            bed_df['end'] = pd.to_numeric(bed_df['end'])
            
            for _, row in bed_df.iterrows():
                regions.append(('interval', row['chrom'], (int(row['start']), int(row['end']))))
            
            print(f"✓ Loaded {len(regions)} regions from {args.bed_file}\n")
        
        except Exception as e:
            print(f"✗ Error reading BED file: {e}")
            sys.exit(1)
    
    elif args.chrom and args.start is not None and args.end is not None:
        # Single region from coordinates
        regions.append(('interval', args.chrom, (args.start, args.end)))
        print(f"✓ Predicting region: {args.chrom}:{args.start}-{args.end}\n")
    
    else:
        print("✗ Must provide either:")
        print("  - --sequence (direct DNA sequence)")
        print("  - --bed-file (genomic regions)")
        print("  - --chrom/--start/--end (single region)")
        sys.exit(1)
    
    # Get tissues
    uberon_terms = tissue_to_uberon(args.tissues)
    if not uberon_terms:
        print("✗ No valid tissues!")
        sys.exit(1)
    
    # Get output types
    output_types = get_output_types(args.output)
    
    # Initialize AlphaGenome
    print("Initializing AlphaGenome model...")
    dna_model = dna_client.create(args.api_key)
    print("✓ Model loaded\n")
    
    print(f"Configuration:")
    print(f"  Tissues: {', '.join(args.tissues)}")
    print(f"  Output types: {len(output_types)}")
    print(f"  Regions/Sequences: {len(regions)}")
    print(f"  Save CSV: {args.save_csv}")
    print(f"  Save NPY: {args.save_npy}\n" + "-"*80 + "\n")
    
    # Process each region
    for idx, region_data in enumerate(regions):
        region_type = region_data[0]
        
        if region_type == 'sequence':
            _, sequence, _ = region_data
            region_name = f"custom_{idx}"
            print(f"[{idx + 1}/{len(regions)}] DNA sequence ({len(sequence)} bp)")
            
            try:
                # Predict for direct sequence
                output = predict_sequence_dna(
                    dna_model, sequence, uberon_terms, output_types,
                    verbose=args.verbose
                )
                
                if output is None:
                    print("  ✗ Failed\n")
                    continue
            
            except Exception as e:
                print(f"  ✗ Error: {e}\n")
                continue
        
        else:  # interval
            _, chrom, (start, end) = region_data
            region_name = f"{chrom}_{start}_{end}"
            print(f"[{idx + 1}/{len(regions)}] {chrom}:{start}-{end}")
            
            try:
                # Check bounds
                if not is_valid_region(chrom, start, end, dna_client.SEQUENCE_LENGTH_1MB):
                    print("  ✗ out of bounds\n")
                    continue
                
                # Predict for interval/genomic region
                output = predict_sequence_region(
                    dna_model, chrom, start, end, uberon_terms, output_types,
                    fasta_file=args.fasta, verbose=args.verbose
                )
                
                if output is None:
                    print("  ✗ Failed\n")
                    continue
            
            except Exception as e:
                print(f"  ✗ Error: {e}\n")
                continue
        
        # Extract predictions
        predictions = extract_sequence_predictions(output, output_types)
        
        if not predictions:
            print("  ✗ No predictions extracted\n")
            continue
        
        # Print summary
        if region_type == 'sequence':
            print_sequence_summary(predictions, region_name, len(sequence))
        else:
            print_sequence_summary(predictions, region_name, end - start)
        
        # Save predictions
        if args.save_csv or args.save_npy:
            if region_type == 'sequence':
                save_sequence_predictions(
                    predictions, 'custom', idx, None, args.output_dir,
                    save_csv=args.save_csv, save_npy=args.save_npy,
                    verbose=args.verbose
                )
            else:
                save_sequence_predictions(
                    predictions, chrom, start, end, args.output_dir,
                    save_csv=args.save_csv, save_npy=args.save_npy,
                    verbose=args.verbose
                )
    
    print("\n" + "="*80)
    print("Sequence Prediction Complete")
    print("="*80)
    
    if args.save_csv or args.save_npy:
        print(f"\n✓ Predictions saved to: {args.output_dir}")
        print(f"\nUsage tips:")
        print(f"  - Load NPY: data = np.load('file.npy')")
        print(f"  - Load CSV: df = pd.read_csv('file.csv')")
        print(f"  - Analyze in Python")
        print(f"  - Visualize with: python plot_sequence_tracks.py --npy-dir {args.output_dir}")

if __name__ == '__main__':
    main()
