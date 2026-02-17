#!/usr/bin/env python3
"""
Simple script to open, view, and analyze AlphaGenome NPY files.

Usage:
    python view_npy.py your_file.npy
    python view_npy.py custom_0_rna_seq.npy
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path

def view_npy(filename, tissue_names=None):
    """Load and display NPY file contents."""
    
    # Load file
    try:
        data = np.load(filename)
    except FileNotFoundError:
        print(f"✗ File not found: {filename}")
        sys.exit(1)
    except Exception as e:
        print(f"✗ Error loading file: {e}")
        sys.exit(1)
    
    # Header
    print("\n" + "="*80)
    print(f"NPY File Viewer")
    print("="*80)
    print(f"\nFile: {filename}")
    
    # Basic info
    print(f"\n📊 SHAPE & SIZE")
    print(f"  Dimensions: {data.shape}")
    print(f"  Positions: {data.shape[0]:,}")
    print(f"  Tissues/Tracks: {data.shape[1] if len(data.shape) > 1 else 1}")
    print(f"  Total values: {data.size:,}")
    print(f"  Memory: {data.nbytes / (1024*1024):.2f} MB")
    
    # Data type
    print(f"\n💾 DATA TYPE")
    print(f"  Type: {data.dtype}")
    print(f"  Format: NumPy array")
    
    # Statistics
    print(f"\n📈 STATISTICS")
    print(f"  Min: {data.min():.6e}")
    print(f"  Max: {data.max():.6e}")
    print(f"  Mean: {data.mean():.6e}")
    print(f"  Median: {np.median(data):.6e}")
    print(f"  Std Dev: {data.std():.6e}")
    
    # Per-position mean
    if len(data.shape) > 1:
        mean_across_tissues = data.mean(axis=1)
        print(f"\n  Mean signal across tissues:")
        print(f"    Min: {mean_across_tissues.min():.6e}")
        print(f"    Max: {mean_across_tissues.max():.6e}")
        print(f"    Mean: {mean_across_tissues.mean():.6e}")
    
    # Per-tissue stats
    if len(data.shape) > 1 and data.shape[1] > 1:
        print(f"\n🔬 PER-TISSUE STATISTICS")
        if tissue_names is None:
            tissue_names = [f"Tissue {i}" for i in range(data.shape[1])]
        
        for i in range(data.shape[1]):
            tissue_data = data[:, i]
            tissue_name = tissue_names[i] if i < len(tissue_names) else f"Tissue {i}"
            print(f"  {tissue_name}:")
            print(f"    Mean: {tissue_data.mean():.6f}")
            print(f"    Max: {tissue_data.max():.6f}")
            print(f"    Min: {tissue_data.min():.6f}")
            print(f"    Std: {tissue_data.std():.6f}")
    
    # Value ranges
    print(f"\n📊 VALUE DISTRIBUTION")
    if len(data.shape) > 1:
        mean_signal = data.mean(axis=1)
    else:
        mean_signal = data
    
    print(f"  < 0.001: {(mean_signal < 0.001).sum():,} positions ({(mean_signal < 0.001).sum()/len(mean_signal)*100:.1f}%)")
    print(f"  0.001-0.01: {((mean_signal >= 0.001) & (mean_signal < 0.01)).sum():,} positions")
    print(f"  0.01-0.1: {((mean_signal >= 0.01) & (mean_signal < 0.1)).sum():,} positions")
    print(f"  > 0.1: {(mean_signal > 0.1).sum():,} positions")
    
    # Peak detection
    threshold = np.percentile(mean_signal, 95)
    peak_count = (mean_signal > threshold).sum()
    print(f"\n🔝 PEAK REGIONS (Top 5%)")
    print(f"  Count: {peak_count:,} positions")
    if peak_count > 0:
        peak_positions = np.where(mean_signal > threshold)[0]
        print(f"  Top 10 peak positions: {peak_positions[:10]}")
        print(f"  Peak values: {mean_signal[peak_positions[:10]]}")
    
    # First/last rows
    print(f"\n📄 FIRST 10 ROWS")
    if len(data.shape) > 1:
        print(f"  Position | {' | '.join([f'Tissue {i}' for i in range(data.shape[1])])}")
        print(f"  " + "-"*50)
        for i in range(min(10, len(data))):
            row_str = f"  {i:>8} | " + " | ".join([f"{val:.6f}" for val in data[i]])
            print(row_str)
    else:
        print(f"  {data[:10]}")
    
    # Create visualization
    print(f"\n📊 GENERATING VISUALIZATIONS...")
    
    fig = plt.figure(figsize=(16, 10))
    
    # Plot 1: Signal profile
    ax1 = plt.subplot(2, 3, 1)
    if len(data.shape) > 1:
        mean_signal = data.mean(axis=1)
    else:
        mean_signal = data
    ax1.plot(mean_signal, linewidth=0.8, color='steelblue')
    ax1.set_title('Signal Profile', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Position (bp)')
    ax1.set_ylabel('Mean Signal')
    ax1.grid(alpha=0.3)
    
    # Plot 2: Per-tissue comparison
    if len(data.shape) > 1 and data.shape[1] > 1:
        ax2 = plt.subplot(2, 3, 2)
        colors = plt.cm.tab10(range(data.shape[1]))
        for i in range(data.shape[1]):
            tissue_name = tissue_names[i] if tissue_names and i < len(tissue_names) else f"Tissue {i}"
            ax2.plot(data[:, i], label=tissue_name, alpha=0.6, color=colors[i])
        ax2.set_title('Tissue-Specific Signals', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Position (bp)')
        ax2.set_ylabel('Signal')
        ax2.legend(loc='best', fontsize=8)
        ax2.grid(alpha=0.3)
    
    # Plot 3: Heatmap
    if len(data.shape) > 1:
        ax3 = plt.subplot(2, 3, 3)
        im = ax3.imshow(data.T, aspect='auto', cmap='RdYlBu_r', origin='lower')
        ax3.set_title('Heatmap (Positions × Tissues)', fontsize=12, fontweight='bold')
        ax3.set_xlabel('Position (bp)')
        ax3.set_ylabel('Tissue')
        plt.colorbar(im, ax=ax3, label='Signal')
    
    # Plot 4: Distribution
    ax4 = plt.subplot(2, 3, 4)
    ax4.hist(mean_signal, bins=100, color='steelblue', alpha=0.7, edgecolor='black')
    ax4.set_title('Value Distribution', fontsize=12, fontweight='bold')
    ax4.set_xlabel('Signal Value')
    ax4.set_ylabel('Frequency (log scale)')
    ax4.set_yscale('log')
    ax4.grid(alpha=0.3)
    
    # Plot 5: Box plot by percentile
    ax5 = plt.subplot(2, 3, 5)
    percentiles = [10, 25, 50, 75, 90, 95, 99]
    percentile_values = [np.percentile(mean_signal, p) for p in percentiles]
    ax5.bar([str(p) for p in percentiles], percentile_values, color='coral', alpha=0.7, edgecolor='black')
    ax5.set_title('Percentile Values', fontsize=12, fontweight='bold')
    ax5.set_xlabel('Percentile')
    ax5.set_ylabel('Signal Value')
    ax5.grid(alpha=0.3, axis='y')
    
    # Plot 6: Cumulative distribution
    ax6 = plt.subplot(2, 3, 6)
    sorted_signal = np.sort(mean_signal)
    cumsum = np.arange(1, len(sorted_signal) + 1) / len(sorted_signal)
    ax6.plot(sorted_signal, cumsum, color='green', linewidth=2)
    ax6.set_title('Cumulative Distribution', fontsize=12, fontweight='bold')
    ax6.set_xlabel('Signal Value')
    ax6.set_ylabel('Cumulative Probability')
    ax6.grid(alpha=0.3)
    
    plt.suptitle(f'NPY File Analysis: {Path(filename).name}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save plot
    output_png = filename.replace('.npy', '_visualization.png')
    plt.savefig(output_png, dpi=150, bbox_inches='tight')
    print(f"✓ Saved visualization: {output_png}")
    
    # Export to CSV
    if len(data.shape) > 1:
        if tissue_names is None:
            tissue_names = [f"Tissue {i}" for i in range(data.shape[1])]
        tissue_names = tissue_names[:data.shape[1]]
        
        df = pd.DataFrame(data, columns=tissue_names)
        df.insert(0, 'Position', range(len(data)))
    else:
        df = pd.DataFrame({'Value': data})
    
    output_csv = filename.replace('.npy', '.csv')
    df.to_csv(output_csv, index=False)
    print(f"✓ Exported to CSV: {output_csv}")
    
    # Print sample of dataframe
    print(f"\n📋 CSV PREVIEW (first 5 rows)")
    print(df.head().to_string(index=False))
    
    print(f"\n✓ Done!")
    print(f"\nFiles created:")
    print(f"  - {output_png}")
    print(f"  - {output_csv}")
    
    return data, df

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python view_npy.py <filename.npy>")
        print("\nExample:")
        print("  python view_npy.py custom_0_rna_seq.npy")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    # Try to guess tissue names from filename
    tissue_names = None
    if 'brain' in filename.lower():
        tissue_names = ['Brain']
    elif 'lung' in filename.lower():
        tissue_names = ['Lung']
    
    # View the file
    data, df = view_npy(filename, tissue_names=tissue_names)
