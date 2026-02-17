#!/usr/bin/env python3
"""
Generate publication-quality plots from AlphaGenome predictions.
Creates visualizations of RNA-seq, chromatin accessibility, and ChIP-seq predictions.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 6)
plt.rcParams['font.size'] = 10

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate plots from AlphaGenome predictions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:

  # Generate all plots
  python plot_alphagenome_results.py results.csv --output plots/

  # Specific plots only
  python plot_alphagenome_results.py results.csv --plots rna_seq atac --output plots/

  # With gene filtering
  python plot_alphagenome_results.py results.csv --output plots/ --exclude-intergenic

  # Compare by biotype
  python plot_alphagenome_results.py results.csv --output plots/ --by-biotype
        """
    )
    parser.add_argument('csv_file', type=str, help='AlphaGenome results CSV')
    parser.add_argument('--output', type=str, default='./alphagenome_plots/',
                       help='Output directory for plots')
    parser.add_argument('--plots', type=str, nargs='+',
                       choices=['rna_seq', 'atac', 'dnase', 'cage', 'chip_tf', 'all'],
                       default=['all'],
                       help='Which plots to generate')
    parser.add_argument('--exclude-intergenic', action='store_true',
                       help='Exclude intergenic regions')
    parser.add_argument('--by-biotype', action='store_true',
                       help='Color plots by gene biotype')
    parser.add_argument('--by-annotation', action='store_true',
                       help='Color plots by annotation type (overlapping/nearest)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for saved plots')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')
    
    return parser.parse_args()

def load_results(csv_file):
    """Load results CSV."""
    try:
        df = pd.read_csv(csv_file)
        print(f"✓ Loaded {len(df)} regions from {csv_file}")
        return df
    except Exception as e:
        print(f"✗ Error loading CSV: {e}")
        sys.exit(1)

def filter_data(df, exclude_intergenic=False):
    """Filter data based on options."""
    original_len = len(df)
    
    if exclude_intergenic:
        df = df[df['gene_name'] != 'no_gene']
        filtered_len = len(df)
        print(f"  Excluded {original_len - filtered_len} intergenic regions")
    
    return df

def plot_rna_seq_distribution(df, output_dir, by_biotype=False, by_annotation=False, dpi=300):
    """Plot RNA-seq prediction distribution."""
    if 'rna_seq_mean' not in df.columns:
        print("⚠ RNA-seq data not available")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('RNA-seq Expression Predictions', fontsize=16, fontweight='bold')
    
    # Histogram
    ax = axes[0, 0]
    ax.hist(df['rna_seq_mean'].dropna(), bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax.set_xlabel('Mean Expression Level')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of RNA-seq Predictions')
    ax.set_yscale('log')
    
    # Scatter: mean vs max
    ax = axes[0, 1]
    if 'rna_seq_max' in df.columns:
        scatter = ax.scatter(df['rna_seq_mean'], df['rna_seq_max'], 
                            alpha=0.6, s=50, c=df['rna_seq_mean'], cmap='viridis')
        ax.set_xlabel('Mean Expression')
        ax.set_ylabel('Max Expression')
        ax.set_title('Mean vs Max Expression')
        plt.colorbar(scatter, ax=ax, label='Mean Expression')
    
    # Box plot by biotype
    ax = axes[1, 0]
    if by_biotype and 'biotype' in df.columns:
        biotypes = df[df['biotype'] != 'intergenic']['biotype'].unique()[:5]
        data_for_box = [df[df['biotype'] == bt]['rna_seq_mean'].dropna() for bt in biotypes]
        bp = ax.boxplot(data_for_box, labels=biotypes, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
        ax.set_xlabel('Gene Biotype')
        ax.set_ylabel('Expression Level')
        ax.set_title('RNA-seq by Gene Type')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    else:
        # Overall distribution
        ax.hist(df['rna_seq_mean'].dropna(), bins=40, color='lightblue', alpha=0.7, edgecolor='black')
        ax.set_xlabel('Expression Level')
        ax.set_ylabel('Frequency')
        ax.set_title('Expression Level Distribution')
    
    # Top expressed regions
    ax = axes[1, 1]
    top_n = 10
    top_genes = df.nlargest(top_n, 'rna_seq_mean')[['gene_name', 'rna_seq_mean']]
    if len(top_genes) > 0:
        genes_display = [str(g)[:20] for g in top_genes['gene_name']]
        ax.barh(genes_display, top_genes['rna_seq_mean'].values, color='steelblue')
        ax.set_xlabel('Mean Expression')
        ax.set_title(f'Top {top_n} Expressed Genes')
        ax.invert_yaxis()
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'rna_seq_distribution.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()

def plot_chromatin_accessibility(df, output_dir, dpi=300):
    """Plot chromatin accessibility (ATAC, DNase)."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Chromatin Accessibility Predictions', fontsize=16, fontweight='bold')
    
    # DNase distribution
    if 'dnase_mean' in df.columns:
        ax = axes[0, 0]
        ax.hist(df['dnase_mean'].dropna(), bins=50, color='coral', alpha=0.7, edgecolor='black')
        ax.set_xlabel('DNase Signal')
        ax.set_ylabel('Frequency')
        ax.set_title('DNase Hypersensitivity Distribution')
        ax.set_yscale('log')
    
    # ATAC distribution
    if 'atac_mean' in df.columns:
        ax = axes[0, 1]
        atac_data = df['atac_mean'].dropna()
        if len(atac_data) > 0:
            ax.hist(atac_data, bins=50, color='mediumseagreen', alpha=0.7, edgecolor='black')
            ax.set_xlabel('ATAC Signal')
            ax.set_ylabel('Frequency')
            ax.set_title('ATAC Accessibility Distribution')
        else:
            ax.text(0.5, 0.5, 'No ATAC data available', ha='center', va='center')
            ax.set_xticks([])
            ax.set_yticks([])
    
    # DNase vs RNA-seq correlation
    if 'dnase_mean' in df.columns and 'rna_seq_mean' in df.columns:
        ax = axes[1, 0]
        valid_data = df[['dnase_mean', 'rna_seq_mean']].dropna()
        if len(valid_data) > 0:
            scatter = ax.scatter(valid_data['dnase_mean'], valid_data['rna_seq_mean'],
                                alpha=0.6, s=50, c=valid_data['dnase_mean'], cmap='RdYlBu_r')
            ax.set_xlabel('DNase Accessibility')
            ax.set_ylabel('RNA-seq Expression')
            ax.set_title('Chromatin Accessibility vs Expression')
            plt.colorbar(scatter, ax=ax, label='DNase')
            
            # Add correlation
            corr = valid_data.corr().iloc[0, 1]
            ax.text(0.05, 0.95, f'Corr: {corr:.3f}', transform=ax.transAxes,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Top accessible regions
    if 'dnase_mean' in df.columns:
        ax = axes[1, 1]
        top_n = 10
        top_accessible = df.nlargest(top_n, 'dnase_mean')[['gene_name', 'dnase_mean']]
        genes_display = [str(g)[:20] for g in top_accessible['gene_name']]
        ax.barh(genes_display, top_accessible['dnase_mean'].values, color='coral')
        ax.set_xlabel('DNase Signal')
        ax.set_title(f'Top {top_n} Accessible Regions')
        ax.invert_yaxis()
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'chromatin_accessibility.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()

def plot_chip_seq_tf_binding(df, output_dir, dpi=300):
    """Plot ChIP-seq transcription factor binding."""
    if 'chip_tf_mean' not in df.columns:
        print("⚠ ChIP-seq data not available")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Transcription Factor Binding (ChIP-seq)', fontsize=16, fontweight='bold')
    
    # Histogram
    ax = axes[0, 0]
    ax.hist(df['chip_tf_mean'].dropna(), bins=50, color='mediumpurple', alpha=0.7, edgecolor='black')
    ax.set_xlabel('TF Binding Signal')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of TF Binding')
    ax.set_yscale('log')
    
    # Mean vs Max
    ax = axes[0, 1]
    if 'chip_tf_max' in df.columns:
        scatter = ax.scatter(df['chip_tf_mean'], df['chip_tf_max'],
                            alpha=0.6, s=50, c=df['chip_tf_mean'], cmap='plasma')
        ax.set_xlabel('Mean TF Binding')
        ax.set_ylabel('Max TF Binding')
        ax.set_title('Mean vs Max TF Binding')
        plt.colorbar(scatter, ax=ax, label='Mean Binding')
    
    # TF vs Expression
    ax = axes[1, 0]
    if 'rna_seq_mean' in df.columns:
        valid_data = df[['chip_tf_mean', 'rna_seq_mean']].dropna()
        if len(valid_data) > 0:
            scatter = ax.scatter(valid_data['chip_tf_mean'], valid_data['rna_seq_mean'],
                                alpha=0.6, s=50, c=valid_data['chip_tf_mean'], cmap='viridis')
            ax.set_xlabel('TF Binding Signal')
            ax.set_ylabel('Gene Expression')
            ax.set_title('TF Binding vs Gene Expression')
            plt.colorbar(scatter, ax=ax, label='TF Binding')
    
    # Top TF binding sites
    ax = axes[1, 1]
    top_n = 10
    top_tf = df.nlargest(top_n, 'chip_tf_mean')[['gene_name', 'chip_tf_mean']]
    genes_display = [str(g)[:20] for g in top_tf['gene_name']]
    ax.barh(genes_display, top_tf['chip_tf_mean'].values, color='mediumpurple')
    ax.set_xlabel('TF Binding Signal')
    ax.set_title(f'Top {top_n} TF Binding Sites')
    ax.invert_yaxis()
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'chip_seq_tf_binding.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()

def plot_all_predictions_heatmap(df, output_dir, dpi=300):
    """Create heatmap of all prediction types."""
    # Select prediction columns
    pred_cols = [col for col in df.columns if '_mean' in col and col != 'gene_distance']
    
    if len(pred_cols) == 0:
        print("⚠ No prediction data for heatmap")
        return
    
    # Normalize data
    df_subset = df[['gene_name'] + pred_cols].copy()
    df_subset = df_subset.dropna(subset=pred_cols, how='all')
    
    if len(df_subset) == 0:
        return
    
    # Take top genes by variance
    top_n = min(20, len(df_subset))
    top_indices = df_subset[pred_cols].var(axis=1).nlargest(top_n).index
    df_heatmap = df_subset.loc[top_indices].copy()
    
    # Set gene names as index
    df_heatmap = df_heatmap.set_index('gene_name')
    
    # Normalize (0-1) for better visualization
    df_heatmap_norm = (df_heatmap - df_heatmap.min()) / (df_heatmap.max() - df_heatmap.min())
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 10))
    
    sns.heatmap(df_heatmap_norm, cmap='RdYlBu_r', cbar_kws={'label': 'Normalized Signal'},
                ax=ax, annot=False, linewidths=0.5)
    
    ax.set_title('All Prediction Types Heatmap (Top Variable Regions)', 
                fontweight='bold', fontsize=14)
    ax.set_xlabel('Prediction Type')
    ax.set_ylabel('Region (Gene)')
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    output_path = Path(output_dir) / 'all_predictions_heatmap.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()

def plot_annotation_summary(df, output_dir, dpi=300):
    """Plot summary of annotation types."""
    if 'annotation_type' not in df.columns:
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle('Annotation Summary', fontsize=14, fontweight='bold')
    
    # Annotation type counts
    ax = axes[0]
    ann_counts = df['annotation_type'].value_counts()
    colors = ['#2ecc71', '#3498db', '#e74c3c', '#95a5a6'][:len(ann_counts)]
    ann_counts.plot(kind='bar', ax=ax, color=colors)
    ax.set_title('Regions by Annotation Type')
    ax.set_xlabel('Annotation Type')
    ax.set_ylabel('Count')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    # Distance distribution for nearest genes
    if 'gene_distance' in df.columns:
        ax = axes[1]
        nearest_mask = df['annotation_type'].str.contains('nearest', na=False)
        distances = df[nearest_mask]['gene_distance'].dropna()
        if len(distances) > 0:
            ax.hist(distances, bins=50, color='skyblue', alpha=0.7, edgecolor='black')
            ax.set_xlabel('Distance to Gene (bp)')
            ax.set_ylabel('Frequency')
            ax.set_title('Distance Distribution for Nearest Genes')
            ax.set_yscale('log')
        else:
            ax.text(0.5, 0.5, 'No nearest gene data', ha='center', va='center')
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'annotation_summary.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()

def plot_biotype_comparison(df, output_dir, dpi=300):
    """Plot predictions by gene biotype."""
    if 'biotype' not in df.columns:
        return
    
    biotypes = df[df['biotype'] != 'intergenic']['biotype'].unique()[:6]
    if len(biotypes) == 0:
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Predictions by Gene Biotype', fontsize=14, fontweight='bold')
    
    # RNA-seq by biotype
    if 'rna_seq_mean' in df.columns:
        ax = axes[0, 0]
        data = [df[df['biotype'] == bt]['rna_seq_mean'].dropna() for bt in biotypes]
        bp = ax.boxplot(data, labels=biotypes, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
        ax.set_ylabel('Expression Level')
        ax.set_title('RNA-seq by Biotype')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # DNase by biotype
    if 'dnase_mean' in df.columns:
        ax = axes[0, 1]
        data = [df[df['biotype'] == bt]['dnase_mean'].dropna() for bt in biotypes]
        bp = ax.boxplot(data, labels=biotypes, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightcoral')
        ax.set_ylabel('Accessibility')
        ax.set_title('DNase by Biotype')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # ChIP-seq by biotype
    if 'chip_tf_mean' in df.columns:
        ax = axes[1, 0]
        data = [df[df['biotype'] == bt]['chip_tf_mean'].dropna() for bt in biotypes]
        bp = ax.boxplot(data, labels=biotypes, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightgreen')
        ax.set_ylabel('TF Binding')
        ax.set_title('TF Binding by Biotype')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Biotype counts
    ax = axes[1, 1]
    biotype_counts = df[df['biotype'] != 'intergenic']['biotype'].value_counts()
    ax.barh(biotype_counts.index, biotype_counts.values, color='steelblue')
    ax.set_xlabel('Count')
    ax.set_title('Gene Count by Biotype')
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'biotype_comparison.png'
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"✓ Saved: {output_path}")
    plt.close()

def main():
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "="*80)
    print("AlphaGenome Results Visualization")
    print("="*80 + "\n")
    
    # Load data
    df = load_results(args.csv_file)
    
    # Filter data
    df = filter_data(df, exclude_intergenic=args.exclude_intergenic)
    print(f"✓ {len(df)} regions to plot\n")
    
    # Determine which plots to generate
    plots_to_gen = args.plots
    if 'all' in plots_to_gen:
        plots_to_gen = ['rna_seq', 'atac', 'dnase', 'cage', 'chip_tf']
    
    print("Generating plots...\n")
    
    # Generate plots
    if 'rna_seq' in plots_to_gen:
        print("  RNA-seq distribution...", end=' ')
        plot_rna_seq_distribution(df, output_dir, by_biotype=args.by_biotype, 
                                 by_annotation=args.by_annotation, dpi=args.dpi)
    
    if 'atac' in plots_to_gen or 'dnase' in plots_to_gen:
        print("  Chromatin accessibility...", end=' ')
        plot_chromatin_accessibility(df, output_dir, dpi=args.dpi)
    
    if 'chip_tf' in plots_to_gen:
        print("  ChIP-seq TF binding...", end=' ')
        plot_chip_seq_tf_binding(df, output_dir, dpi=args.dpi)
    
    # Additional plots
    print("\n  All predictions heatmap...", end=' ')
    plot_all_predictions_heatmap(df, output_dir, dpi=args.dpi)
    
    print("\n  Annotation summary...", end=' ')
    plot_annotation_summary(df, output_dir, dpi=args.dpi)
    
    if args.by_biotype:
        print("\n  Biotype comparison...", end=' ')
        plot_biotype_comparison(df, output_dir, dpi=args.dpi)
    
    print("\n\n" + "="*80)
    print("Summary Statistics")
    print("="*80 + "\n")
    
    # Print summary stats
    if 'rna_seq_mean' in df.columns:
        print("RNA-seq Expression:")
        print(f"  Mean: {df['rna_seq_mean'].mean():.6f}")
        print(f"  Median: {df['rna_seq_mean'].median():.6f}")
        print(f"  Std Dev: {df['rna_seq_mean'].std():.6f}")
    
    if 'dnase_mean' in df.columns:
        print("\nDNase Accessibility:")
        print(f"  Mean: {df['dnase_mean'].mean():.6f}")
        print(f"  Median: {df['dnase_mean'].median():.6f}")
        print(f"  Std Dev: {df['dnase_mean'].std():.6f}")
    
    if 'chip_tf_mean' in df.columns:
        print("\nTF Binding (ChIP-seq):")
        print(f"  Mean: {df['chip_tf_mean'].mean():.6f}")
        print(f"  Median: {df['chip_tf_mean'].median():.6f}")
        print(f"  Std Dev: {df['chip_tf_mean'].std():.6f}")
    
    print(f"\n✓ All plots saved to: {output_dir}")

if __name__ == '__main__':
    main()
