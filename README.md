# prediction_AG
Prediction examples using AlphaGenome

# To run the prediction with AlphaGenome with intervals:

## Step 1: Prepare Your Data

### Get GTF File
```bash
# Download GTF (GENCODE)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.basic.annotation.gtf.gz
gunzip gencode.vM38.basic.annotation.gtf.gz

```

### Prepare BED File OR Use you bed file in this format
Format: chrom, start, end (tab-separated)
```
chr6	49163874	49165042
chr10	5743395	5746168
chr2	18720549	18722409
```

---

## Step 2: Run AlphaGenome Analysis

### Option A: With Nearest Gene Annotation (RECOMMENDED)
```bash
python alphagenome_gtf_pipeline_interval.py your_regions.bed gencode.vM38.basic.annotation.gtf \
  --api-key YOUR_API_KEY \
  --tissues intestine \
  --buffer 10000 \
  --output results.csv
```

### Option B: Simple GTF Annotation + AlphaGenome
```bash
python alphagenome_gtf_pipeline.py your_regions.bed gencode.vM38.basic.annotation.gtf \
  --api-key YOUR_API_KEY \
  --tissues brain \
  --output results.csv
```

### Option C: No GTF (Just Predictions)
```bash
python alphagenome_final.py your_regions.bed \
  --api-key YOUR_API_KEY \
  --tissues brain \
  --output results.csv
```

---

## Step 3: Generate Plots

### Basic Plots
```bash
python plot_alphagenome_results.py results.csv --output plots/
```

This generates:
- ✓ RNA-seq distribution
- ✓ Chromatin accessibility (ATAC, DNase)
- ✓ ChIP-seq (TF binding)
- ✓ All predictions heatmap
- ✓ Annotation summary
- ✓ Statistical summary

### Enhanced Plots with Biotype Analysis
```bash
python plot_alphagenome_results.py results.csv \
  --output plots/ \
  --exclude-intergenic \
  --by-biotype
```

### Publication Quality
```bash
python plot_alphagenome_results.py results.csv \
  --output publication_figs/ \
  --by-biotype \
  --exclude-intergenic \
  --dpi 600
```


#To run AlphaGenome DNA Sequence Predictions

**Three ways to provide sequences:**
1. ✅ Direct DNA strings (ATCG...)
2. ✅ Genomic coordinates (chrom:start-end)
3. ✅ FASTA files

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

#For more information on analysis interpretation please look into:
#DNA_SEQUENCE_GUIDE.md
#README_alphagenome_interval.md

#To make the plots from NPY file:
#use VIEW_NPY_QUICK_START.md
