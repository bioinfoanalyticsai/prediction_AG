import argparse
import pandas as pd
import json
import numpy as np
from alphagenome.data import genome
from alphagenome.models import dna_client
import warnings
import sys
import os
from typing import Dict, List, Tuple

warnings.filterwarnings('ignore')

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='GTF annotation + AlphaGenome pipeline with nearest gene support',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:

  # Find overlapping genes, or nearest if no overlap
  python alphagenome_gtf_pipeline_nearest.py regions.bed mm10.gtf \
    --api-key YOUR_KEY --tissues intestine

  # Only overlapping genes (no nearest)
  python alphagenome_gtf_pipeline_nearest.py regions.bed mm10.gtf \
    --api-key YOUR_KEY --tissues brain --overlap-only

  # With buffer for extended search
  python alphagenome_gtf_pipeline_nearest.py regions.bed mm10.gtf \
    --api-key YOUR_KEY --tissues lung --buffer 10000
        """
    )
    parser.add_argument('bed_file', type=str, help='Input BED file')
    parser.add_argument('gtf_file', type=str, help='GTF file with gene annotations')
    parser.add_argument('--api-key', type=str, required=True, help='AlphaGenome API key')
    parser.add_argument('--tissues', type=str, nargs='+', default=['brain'],
                       help='Tissues: brain, lung, colon, intestine, kidney')
    parser.add_argument('--output', type=str, default='alphagenome_gtf_results_nearest.csv',
                       help='Output CSV')
    parser.add_argument('--sequence-length', type=str, default='SEQUENCE_LENGTH_100KB',
                       choices=['SEQUENCE_LENGTH_16KB', 'SEQUENCE_LENGTH_100KB',
                               'SEQUENCE_LENGTH_500KB', 'SEQUENCE_LENGTH_1MB'])
    parser.add_argument('--feature-type', type=str, default='gene',
                       choices=['gene', 'exon', 'transcript'])
    parser.add_argument('--overlap', type=str, default='any',
                       choices=['any', 'strict', 'partial'])
    parser.add_argument('--buffer', type=int, default=0, help='Buffer region (bp)')
    parser.add_argument('--overlap-only', action='store_true',
                       help='Only report overlapping genes (no nearest)')
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

def parse_gtf(gtf_file, feature_type='gene'):
    """Parse GTF file and return genes indexed by chromosome."""
    genes = {}
    print(f"Parsing GTF file: {gtf_file}")
    
    try:
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                if fields[2] != feature_type:
                    continue
                
                chrom = fields[0]
                start = int(fields[3]) - 1
                end = int(fields[4])
                strand = fields[6]
                
                attributes = parse_gtf_attributes(fields[8])
                
                if chrom not in genes:
                    genes[chrom] = []
                
                gene_info = {
                    'gene_id': attributes.get('gene_id', 'unknown'),
                    'gene_name': attributes.get('gene_name', attributes.get('gene_id', 'unknown')),
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'biotype': attributes.get('gene_biotype', 'unknown'),
                }
                
                genes[chrom].append(gene_info)
    except Exception as e:
        print(f"✗ Error parsing GTF: {e}")
        sys.exit(1)
    
    # Sort genes by position
    for chrom in genes:
        genes[chrom].sort(key=lambda x: x['start'])
    
    total_genes = sum(len(g) for g in genes.values())
    print(f"✓ Found {total_genes} {feature_type} features\n")
    return genes

def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string."""
    attributes = {}
    pairs = attr_string.strip().rstrip(';').split(';')
    for pair in pairs:
        pair = pair.strip()
        if ' ' in pair:
            key, value = pair.split(' ', 1)
            attributes[key.strip()] = value.strip().strip('"')
    return attributes

def calculate_distance(region_start, region_end, gene_start, gene_end):
    """Calculate distance between region and gene."""
    if region_start <= gene_end and region_end >= gene_start:
        return 0  # Overlapping
    elif region_end < gene_start:
        return gene_start - region_end
    else:
        return region_start - gene_end

def find_overlapping_gene(chrom, start, end, genes, overlap_type='any', buffer=0):
    """Find overlapping gene."""
    if chrom not in genes:
        return None
    
    start_buffered = start - buffer
    end_buffered = end + buffer
    
    for gene in genes[chrom]:
        if gene['end'] < start_buffered or gene['start'] > end_buffered:
            continue
        
        # Check overlap
        if overlap_type == 'any':
            if start_buffered < gene['end'] and end_buffered > gene['start']:
                return gene
        elif overlap_type == 'strict':
            if gene['start'] <= start_buffered and end_buffered <= gene['end']:
                return gene
        elif overlap_type == 'partial':
            overlap_start = max(start_buffered, gene['start'])
            overlap_end = min(end_buffered, gene['end'])
            overlap_length = max(0, overlap_end - overlap_start)
            region_length = end_buffered - start_buffered
            if overlap_length >= (region_length * 0.5):
                return gene
    
    return None

def find_nearest_gene(chrom, start, end, genes):
    """
    Find the nearest gene to the region (if not overlapping).
    Returns (gene_info, distance, direction)
    """
    if chrom not in genes:
        return None, None, None
    
    nearest_gene = None
    nearest_distance = float('inf')
    direction = None
    
    for gene in genes[chrom]:
        distance = calculate_distance(start, end, gene['start'], gene['end'])
        
        # Skip overlapping genes (distance 0)
        if distance == 0:
            continue
        
        if distance < nearest_distance:
            nearest_distance = distance
            nearest_gene = gene
            
            # Determine direction
            if end < gene['start']:
                direction = 'upstream'
            else:
                direction = 'downstream'
    
    return nearest_gene, nearest_distance, direction

def extract_predictions(output):
    """Extract predictions from AlphaGenome Output object."""
    pred_dict = {}
    
    try:
        tracks = {
            'rna_seq': output.rna_seq,
            'atac': output.atac,
            'dnase': output.dnase,
            'cage': output.cage,
            'chip_histone': output.chip_histone,
            'chip_tf': output.chip_tf,
        }
        
        for track_name, track_data in tracks.items():
            if track_data is None:
                pred_dict[f'{track_name}_available'] = False
            else:
                if hasattr(track_data, 'values') and track_data.values is not None:
                    values = track_data.values
                    
                    if values.size > 0 and values.shape[1] > 0:
                        flat_values = values.flatten()
                        pred_dict[f'{track_name}_mean'] = float(np.mean(flat_values))
                        pred_dict[f'{track_name}_max'] = float(np.max(flat_values))
                        pred_dict[f'{track_name}_std'] = float(np.std(flat_values))
                        pred_dict[f'{track_name}_available'] = True
                    else:
                        pred_dict[f'{track_name}_available'] = False
        
        return pred_dict
    except Exception as e:
        return {'error': str(e)}

def is_valid_region(chrom, start, end, context_size):
    """Check if region is within chromosome bounds."""
    if chrom not in MM10_CHROM_SIZES:
        return False
    chrom_size = MM10_CHROM_SIZES[chrom]
    padded_start = max(0, start - context_size // 2)
    padded_end = min(chrom_size, end + context_size // 2)
    return not (padded_start >= chrom_size or padded_end <= 0)

def tissue_to_uberon(tissues):
    """Convert tissue names to UBERON codes."""
    uberon_terms = []
    for tissue in tissues:
        if tissue.lower() in TISSUE_MAPPING:
            uberon_terms.append(TISSUE_MAPPING[tissue.lower()])
        else:
            print(f"⚠ Unknown tissue '{tissue}'. Supported: {', '.join(TISSUE_MAPPING.keys())}")
    return uberon_terms

def main():
    args = parse_arguments()
    
    print("\n" + "="*80)
    print("AlphaGenome GTF Annotation + Prediction Pipeline (with Nearest Gene)")
    print("="*80 + "\n")
    
    # Step 1: Parse GTF
    genes = parse_gtf(args.gtf_file, feature_type=args.feature_type)
    
    # Step 2: Read BED
    print(f"Reading BED file: {args.bed_file}")
    try:
        bed_df = pd.read_csv(args.bed_file, sep='\t', header=None, comment='#')
        if bed_df.shape[1] < 3:
            raise ValueError("BED file must have at least 3 columns")
        bed_df.columns = ['chrom', 'start', 'end'] + [f'col_{i}' for i in range(3, bed_df.shape[1])]
        bed_df['start'] = pd.to_numeric(bed_df['start'])
        bed_df['end'] = pd.to_numeric(bed_df['end'])
        print(f"✓ Loaded {len(bed_df)} regions\n")
    except Exception as e:
        print(f"✗ Error reading BED: {e}")
        sys.exit(1)
    
    # Step 3: Initialize AlphaGenome
    uberon_terms = tissue_to_uberon(args.tissues)
    if not uberon_terms:
        print("✗ No valid tissues!")
        return
    
    print("Initializing AlphaGenome model...")
    dna_model = dna_client.create(args.api_key)
    print("✓ Model loaded\n")
    
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[args.sequence_length]
    seq_length_bp = int(args.sequence_length.split('_')[-1].replace('KB', '000').replace('MB', '000000'))
    
    print("Configuration:")
    print(f"  Tissues: {', '.join(args.tissues)}")
    print(f"  Sequence length: {args.sequence_length}")
    print(f"  Overlap type: {args.overlap}")
    print(f"  Buffer: {args.buffer} bp")
    print(f"  Find nearest: {not args.overlap_only}\n" + "-"*80 + "\n")
    
    # Step 4: Run predictions with gene annotation
    results = []
    total_regions = len(bed_df)
    stats = {'overlapping': 0, 'nearest': 0, 'none': 0, 'success': 0, 'error': 0}
    
    for idx, row in bed_df.iterrows():
        region_size = row['end'] - row['start']
        
        # Find gene information
        overlapping_gene = find_overlapping_gene(
            row['chrom'], row['start'], row['end'], genes,
            overlap_type=args.overlap, buffer=args.buffer
        )
        
        if overlapping_gene:
            gene = overlapping_gene
            annotation_type = 'overlapping'
            gene_distance = 0
            direction = 'N/A'
            stats['overlapping'] += 1
        elif not args.overlap_only:
            nearest_gene, distance, direction = find_nearest_gene(
                row['chrom'], row['start'], row['end'], genes
            )
            if nearest_gene:
                gene = nearest_gene
                annotation_type = f'nearest_{direction}'
                gene_distance = distance
                stats['nearest'] += 1
            else:
                gene = None
                annotation_type = 'none'
                gene_distance = None
                direction = None
                stats['none'] += 1
        else:
            gene = None
            annotation_type = 'none'
            gene_distance = None
            direction = None
            stats['none'] += 1
        
        # Build gene info display
        gene_info = f"({gene['gene_name']})" if gene else "(no gene)"
        if annotation_type.startswith('nearest'):
            gene_info += f" {direction} {gene_distance}bp"
        
        region_str = f"{row['chrom']}:{row['start']}-{row['end']} {gene_info}"
        print(f"[{idx + 1:>3}/{total_regions}] {region_str:<65} ... ", end='', flush=True)
        
        # Check bounds
        if not is_valid_region(row['chrom'], row['start'], row['end'], seq_length_bp):
            print("✗ out of bounds")
            result_row = {
                'chrom': row['chrom'],
                'start': row['start'],
                'end': row['end'],
                'region_size': region_size,
                'status': 'skipped',
                'annotation_type': annotation_type,
            }
            
            if gene:
                result_row['gene_name'] = gene['gene_name']
                result_row['gene_id'] = gene['gene_id']
                result_row['biotype'] = gene['biotype']
                result_row['strand'] = gene['strand']
                result_row['gene_distance'] = gene_distance if gene_distance else ''
                result_row['direction'] = direction if direction else ''
            
            results.append(result_row)
            continue
        
        try:
            interval = genome.Interval(
                chromosome=row['chrom'],
                start=int(row['start']),
                end=int(row['end'])
            )
            interval = interval.resize(sequence_length)
            
            output = dna_model.predict_interval(
                interval=interval,
                requested_outputs=[
                    dna_client.OutputType.RNA_SEQ,
                    dna_client.OutputType.ATAC,
                    dna_client.OutputType.DNASE,
                    dna_client.OutputType.CAGE,
                    dna_client.OutputType.CHIP_HISTONE,
                    dna_client.OutputType.CHIP_TF,
                    dna_client.OutputType.SPLICE_SITES,
                    dna_client.OutputType.SPLICE_SITE_USAGE,
                    dna_client.OutputType.SPLICE_JUNCTIONS,
                    dna_client.OutputType.CONTACT_MAPS,
                    dna_client.OutputType.PROCAP,
                ],
                ontology_terms=uberon_terms,
                organism=dna_client.Organism.MUS_MUSCULUS
            )
            
            predictions = extract_predictions(output)
            
            result_row = {
                'chrom': row['chrom'],
                'start': row['start'],
                'end': row['end'],
                'region_size': region_size,
                'status': 'success',
                'annotation_type': annotation_type,
            }
            
            # Add gene information
            if gene:
                result_row['gene_name'] = gene['gene_name']
                result_row['gene_id'] = gene['gene_id']
                result_row['biotype'] = gene['biotype']
                result_row['strand'] = gene['strand']
                result_row['gene_start'] = gene['start']
                result_row['gene_end'] = gene['end']
                result_row['gene_distance'] = gene_distance if gene_distance else ''
                result_row['direction'] = direction if direction else ''
            else:
                result_row['gene_name'] = 'no_gene'
                result_row['gene_id'] = ''
                result_row['biotype'] = 'none'
            
            # Add predictions
            result_row.update(predictions)
            results.append(result_row)
            stats['success'] += 1
            
            print("✓")
            
        except Exception as e:
            error_msg = str(e)[:50]
            print(f"✗ {error_msg}")
            result_row = {
                'chrom': row['chrom'],
                'start': row['start'],
                'end': row['end'],
                'region_size': region_size,
                'status': 'error',
                'annotation_type': annotation_type,
            }
            
            if gene:
                result_row['gene_name'] = gene['gene_name']
                result_row['gene_id'] = gene['gene_id']
                result_row['biotype'] = gene['biotype']
            
            results.append(result_row)
            stats['error'] += 1
    
    # Save results
    print("\n" + "-"*80)
    df_results = pd.DataFrame(results)
    df_results.to_csv(args.output, index=False)
    
    print(f"\n✓ DONE! Results saved to: {args.output}")
    print(f"\nSummary:")
    print(f"  Total regions: {total_regions}")
    print(f"  Overlapping genes: {stats['overlapping']}")
    print(f"  Nearest genes: {stats['nearest']}")
    print(f"  No genes found: {stats['none']}")
    print(f"  Successfully predicted: {stats['success']}")
    print(f"  Errors: {stats['error']}")
    
    # Show sample results
    if stats['success'] > 0:
        print(f"\nSample results (first 3 successful):")
        success_df = df_results[df_results['status'] == 'success'].head(3)
        cols_to_show = ['chrom', 'start', 'end', 'gene_name', 'annotation_type', 'gene_distance', 'direction', 'rna_seq_mean', 'dnase_mean', 'chip_tf_mean', 'splice_sites', 'splice_site_usage', 'splice_junction', 'contact_maps', 'procap']
        display_df = success_df[[c for c in cols_to_show if c in success_df.columns]]
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        print(display_df.to_string())

if __name__ == '__main__':
    main()
