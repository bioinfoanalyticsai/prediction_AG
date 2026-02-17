# AlphaGenome DNA Sequence Predictions

## Overview

AlphaGenome's `predict_sequence()` method now works with actual DNA sequences!

**Three ways to provide sequences:**
1. ✅ Direct DNA strings (ATCG...)
2. ✅ Genomic coordinates (chrom:start-end)
3. ✅ FASTA files

---

## Quick Start

### Method 1: Direct DNA Sequence

```bash
python alphagenome_predict_sequence.py \
  --sequence "GATTACAGATTACA" \
  --api-key YOUR_KEY \
  --tissues brain lung \
  --save-npy
```

The sequence is automatically padded with 'N's to valid AlphaGenome lengths (16KB, 100KB, 500KB, 1MB).

### Method 2: Genomic Coordinates

```bash
# Uses predict_interval (fallback)
python alphagenome_predict_sequence.py \
  --chrom chr6 --start 49163874 --end 49165042 \
  --api-key YOUR_KEY \
  --tissues brain
```

### Method 3: FASTA File

```bash
# Fetches sequence from FASTA
python alphagenome_predict_sequence.py \
  --chrom chr6 --start 49163874 --end 49165042 \
  --fasta mm10.fa.gz \
  --api-key YOUR_KEY \
  --tissues brain
```

---

## Full Examples

### Example 1: Simple Sequence

```bash
python alphagenome_predict_sequence.py \
  --sequence "GATTACA" \
  --api-key KEY \
  --tissues lung brain \
  --output rna_seq cage dnase atac \
  --save-npy --save-csv
```

**Output:**
```
sequence_predictions/
├── custom_0_rna_seq.npy
├── custom_0_rna_seq.csv
├── custom_0_cage.npy
├── custom_0_cage.csv
├── custom_0_dnase.npy
├── custom_0_dnase.csv
└── custom_0_atac.npy
```

### Example 2: Real Sequence with Tissues

```bash
# Predict for a known regulatory sequence
python alphagenome_predict_sequence.py \
  --sequence "TTAAACAGATCTGAGGTCTCTGAGGCTCCAGAGGATGTAAGAGCCTTCCGCGTTACTGTGCAATGCCGGCCTAAGAGCCGTGTGTAGCAAAAGAGCGCTCCG" \
  --api-key KEY \
  --tissues brain lung intestine \
  --save-npy

# Check predictions
python plot_sequence_tracks.py \
  --npy-file custom_0_rna_seq.npy \
  --chrom custom --start 0 --end 100
```

### Example 3: Batch Processing

```bash
# Predict from BED file
python alphagenome_predict_sequence.py \
  --bed-file regions.bed \
  --fasta mm10.fa.gz \
  --api-key KEY \
  --tissues brain \
  --save-npy
```

---

## Input Specifications

### Sequence Requirements

- **Format:** ACGT or N (upper case)
- **Length:** Any length (automatically padded)
- **Valid:** A, C, G, T, N (case-insensitive, converted to uppercase)

### Padding

Sequences are automatically padded to valid AlphaGenome lengths:

| Length Name | Size |
|------------|------|
| SEQUENCE_LENGTH_16KB | 16,384 bp |
| SEQUENCE_LENGTH_100KB | 102,400 bp |
| SEQUENCE_LENGTH_500KB | 512,000 bp |
| SEQUENCE_LENGTH_1MB | 1,048,576 bp |

**Example:**
```python
# Input
sequence = "GATTACA"  # 7 bp

# Auto-padded (example to 100KB)
padded = "N...N" + "GATTACA" + "N...N"  # Total: 102,400 bp
```

---

## Code Example

### Direct Python Usage

```python
from alphagenome.models import dna_client

# Initialize model
dna_model = dna_client.create(api_key='YOUR_KEY')

# Your sequence
sequence = 'GATTACA'.center(dna_client.SEQUENCE_LENGTH_1MB, 'N')

# Predict
output = dna_model.predict_sequence(
    sequence=sequence,
    requested_outputs=[
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.DNASE,
        dna_client.OutputType.ATAC,
        dna_client.OutputType.CAGE,
    ],
    ontology_terms=[
        'UBERON:0002048',  # Lung
        'UBERON:0000955',  # Brain
    ],
)

# Access results
print(f'RNA-seq shape: {output.rna_seq.values.shape}')
print(f'DNase shape: {output.dnase.values.shape}')
print(f'ATAC shape: {output.atac.values.shape}')
print(f'CAGE shape: {output.cage.values.shape}')
```

---

## Output Format

### NumPy Files

**Shape:** (n_positions, n_tissues)

```python
import numpy as np

# Load
data = np.load('custom_0_rna_seq.npy')
print(data.shape)  # e.g., (1048576, 2) for 1MB + 2 tissues

# Access
tissue_1_data = data[:, 0]  # All positions, tissue 1
position_1000 = data[1000, :]  # Position 1000, all tissues
```

### CSV Files

**Columns:** position, track_0, track_1, ...

```
position,track_0,track_1
0,0.001,0.002
1,0.001,0.003
2,0.002,0.002
...
```

---

## Analyzing Predictions

### Find Peaks

```python
import numpy as np

data = np.load('rna_seq.npy')
mean_signal = np.mean(data, axis=1)

# Find peak positions
peak_threshold = np.percentile(mean_signal, 90)
peak_positions = np.where(mean_signal > peak_threshold)[0]

print(f"Peak positions: {peak_positions}")
```

### Compare Tissues

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.load('rna_seq.npy')
tissue_1 = data[:, 0]
tissue_2 = data[:, 1]

plt.plot(tissue_1, label='Tissue 1')
plt.plot(tissue_2, label='Tissue 2')
plt.legend()
plt.show()
```

### Identify Regulatory Elements

```python
import numpy as np

rna = np.load('rna_seq.npy').mean(axis=1)
dnase = np.load('dnase.npy').mean(axis=1)
atac = np.load('atac.npy').mean(axis=1)

# Find accessible regions (high DNase + ATAC)
accessible = (dnase > np.percentile(dnase, 75)) & \
             (atac > np.percentile(atac, 75))

# Find promoters (high accessibility + expression)
promoters = accessible & (rna > np.percentile(rna, 75))

print(f"Accessible positions: {np.sum(accessible)}")
print(f"Promoter positions: {np.sum(promoters)}")
```

---

## Common Use Cases

### 1. Score Variants

```bash
# Get wild-type signal
python alphagenome_predict_sequence.py \
  --sequence "GATTACAGATTACA" \
  --api-key KEY --tissues brain --save-npy

# Get mutant signal (with SNP)
python alphagenome_predict_sequence.py \
  --sequence "GATTACAGATTACA".replace('G', 'A', 1) \
  --api-key KEY --tissues brain --save-npy

# Compare impacts
import numpy as np
wt = np.load('custom_0_rna_seq.npy').mean(axis=1)
mut = np.load('custom_1_rna_seq.npy').mean(axis=1)
delta = np.mean(np.abs(wt - mut))
print(f"Average impact: {delta}")
```

### 2. Analyze Enhancers

```bash
# Get enhancer sequence from FASTA
python alphagenome_predict_sequence.py \
  --chrom chr6 --start 49163874 --end 49165042 \
  --fasta mm10.fa.gz \
  --api-key KEY \
  --tissues brain lung intestine \
  --save-npy

# Visualize
python plot_sequence_tracks.py \
  --npy-dir sequence_predictions/ \
  --output enhancer_plots/
```

### 3. Predict Regulatory Activity

```bash
# Sequence of interest
python alphagenome_predict_sequence.py \
  --sequence "YOUR_SEQUENCE_HERE" \
  --api-key KEY \
  --tissues brain lung colon intestine kidney \
  --output rna_seq dnase cage \
  --save-npy

# Python analysis
import numpy as np

for tissue in ['brain', 'lung', 'colon', 'intestine', 'kidney']:
    data = np.load(f'custom_0_rna_seq_{tissue}.npy')
    print(f"{tissue}: mean expression = {data.mean():.6f}")
```

---

## Troubleshooting

### Issue: "Sequence invalid"

**Check:**
- Only ACGT and N characters
- Correct case (converted automatically but verify)
- No special characters or spaces

```python
# Clean sequence
sequence = "gaTTacA"  # Automatically uppercased
sequence = sequence.upper().replace('-', 'N')  # Remove gaps
```

### Issue: Results seem wrong

**Check:**
- Verify tissue ontology terms are correct
- Confirm output types match what you expect
- Compare with genomic coordinates (should be similar if sequence matches)

### Issue: Very long sequences slow

**Solutions:**
- Use smaller sequences (< 500KB)
- Process in batches
- Cache results

---

## Performance

| Sequence Length | Est. Time | Memory |
|-----------------|-----------|--------|
| 100 bp | < 1 sec | ~50 MB |
| 1 KB | 1-2 sec | ~100 MB |
| 100 KB | 5-10 sec | ~500 MB |
| 1 MB | 30-60 sec | ~2 GB |

---

## Tips & Tricks

1. **Center your sequence:**
   ```python
   padded = sequence.center(dna_client.SEQUENCE_LENGTH_1MB, 'N')
   ```

2. **Use consistent tissue names:**
   ```python
   tissues = ['brain', 'lung', 'intestine']  # All lowercase
   ```

3. **Save both NPY and CSV:**
   ```bash
   --save-npy --save-csv  # Double save for flexibility
   ```

4. **Analyze in batches:**
   ```bash
   for seq in sequences:
     python alphagenome_predict_sequence.py \
       --sequence "$seq" --api-key KEY --save-npy
   done
   ```

---

## Next Steps

1. **Run predictions:**
   ```bash
   python alphagenome_predict_sequence.py --sequence "GATTACA" --api-key KEY --tissues brain --save-npy
   ```

2. **Analyze results:**
   ```python
   import numpy as np
   data = np.load('custom_0_rna_seq.npy')
   print(data.shape)
   ```

3. **Visualize:**
   ```bash
   python plot_sequence_tracks.py --npy-dir sequence_predictions/ --output plots/
   ```

4. **Integrate into workflows:**
   - Variant effect prediction
   - Enhancer discovery
   - Regulatory variant scoring

---

## Files Updated

✅ `alphagenome_predict_sequence.py` - Supports DNA sequences  
✅ `plot_sequence_tracks.py` - Visualize results

---

Happy predicting! 🧬
