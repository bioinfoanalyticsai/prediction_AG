# View NPY Files - Quick Start

## Simplest Way: Use the Script

```bash
python view_npy.py custom_0_rna_seq.npy
```

**Output:**
- 📊 Detailed statistics
- 📈 6 visualization plots
- 📋 CSV export
- ✓ Summary report

---

## Manual Python (3 Lines)

```python
import numpy as np

data = np.load('custom_0_rna_seq.npy')
print(data.shape)      # See dimensions
print(data)            # View all values
```

---

## Common Tasks

### View Dimensions
```python
import numpy as np
data = np.load('custom_0_rna_seq.npy')
print(data.shape)  # (1048576, 2) = positions × tissues
```

### View Statistics
```python
import numpy as np
data = np.load('custom_0_rna_seq.npy')
print(f"Min: {data.min()}")
print(f"Max: {data.max()}")
print(f"Mean: {data.mean()}")
```

### Get Specific Values
```python
import numpy as np
data = np.load('custom_0_rna_seq.npy')

# Position 1000, tissue 0 (brain)
value = data[1000, 0]

# All positions, tissue 0
brain_signal = data[:, 0]

# Position 500, all tissues
all_tissues = data[500, :]
```

### Plot Signal
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.load('custom_0_rna_seq.npy')
mean = data.mean(axis=1)  # Average across tissues

plt.plot(mean)
plt.xlabel('Position (bp)')
plt.ylabel('Signal')
plt.title('RNA-seq Profile')
plt.savefig('plot.png')
plt.show()
```

### Export to CSV
```python
import numpy as np
import pandas as pd

data = np.load('custom_0_rna_seq.npy')
df = pd.DataFrame(data, columns=['Brain', 'Lung'])
df.to_csv('output.csv', index=False)
```

### Find Peaks
```python
import numpy as np

data = np.load('custom_0_rna_seq.npy')
mean = data.mean(axis=1)

# Find top 5% positions
threshold = np.percentile(mean, 95)
peaks = np.where(mean > threshold)[0]

print(f"Peak positions: {peaks}")
```

---

## File Structure

```
custom_0_rna_seq.npy
├─ Shape: (1048576, 2)
│  ├─ 1,048,576 rows = genomic positions
│  └─ 2 columns = tissues (brain, lung)
├─ Data type: float32
└─ Memory: ~8 MB
```

---

## Access Patterns

```python
data = np.load('file.npy')

# By position
data[0]        # Position 0, all tissues
data[100]      # Position 100, all tissues
data[0:10]     # Positions 0-9, all tissues

# By tissue  
data[:, 0]     # All positions, tissue 0 (brain)
data[:, 1]     # All positions, tissue 1 (lung)

# Specific cell
data[1000, 0]  # Position 1000, tissue 0
```

---

## Comparison: Multiple Files

```python
import numpy as np

# Load multiple track types
rna_seq = np.load('custom_0_rna_seq.npy').mean(axis=1)
dnase = np.load('custom_0_dnase.npy').mean(axis=1)
atac = np.load('custom_0_atac.npy').mean(axis=1)

# Find regions high in all
high_all = (rna_seq > 0.01) & (dnase > 0.05) & (atac > 0.02)
positions = np.where(high_all)[0]

print(f"Regions high in all: {len(positions)}")
```

---

## Visualization

### Line Plot
```python
import matplotlib.pyplot as plt
import numpy as np

data = np.load('file.npy')
plt.plot(data.mean(axis=1))
plt.savefig('plot.png')
```

### Heatmap
```python
import matplotlib.pyplot as plt
import numpy as np

data = np.load('file.npy')
plt.imshow(data.T, aspect='auto', cmap='RdYlBu_r')
plt.colorbar()
plt.savefig('heatmap.png')
```

### Compare Tissues
```python
import matplotlib.pyplot as plt
import numpy as np

data = np.load('file.npy')
plt.plot(data[:, 0], label='Tissue 0')
plt.plot(data[:, 1], label='Tissue 1')
plt.legend()
plt.savefig('comparison.png')
```

---

## Script Usage

```bash
# View and analyze
python view_npy.py custom_0_rna_seq.npy

# Creates:
# - custom_0_rna_seq_visualization.png (6 plots)
# - custom_0_rna_seq.csv (spreadsheet)
# - Console output (statistics)
```

---

## Files Created

1. **PNG file** - 6 visualization plots
2. **CSV file** - Spreadsheet format (open in Excel)
3. **Console output** - Statistics and summary

---

## Quick Summary

| Task | Command |
|------|---------|
| View everything | `python view_npy.py file.npy` |
| Load in Python | `data = np.load('file.npy')` |
| See shape | `print(data.shape)` |
| Get stats | `print(data.mean())` |
| Plot | `plt.plot(data.mean(axis=1))` |
| Export CSV | `pd.DataFrame(data).to_csv('out.csv')` |

---

**Start with: `python view_npy.py your_file.npy`** 🎉
