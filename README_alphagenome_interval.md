# Complete AlphaGenome Analysis Workflow with Visualization

## End-to-End Pipeline

This guide shows the complete workflow from raw BED file to publication-ready plots.

---

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

---

## Step 4: Interpret Results

### Check Summary Statistics
```
RNA-seq Expression:
  Mean: 0.015919
  Median: 0.005123
  Std Dev: 0.045321

DNase Accessibility:
  Mean: 0.052443
  Median: 0.021345
  Std Dev: 0.134567

TF Binding (ChIP-seq):
  Mean: 78.234
  Median: 35.678
  Std Dev: 156.789
```

### Examine Plots
1. **rna_seq_distribution.png** - Gene expression levels
2. **chromatin_accessibility.png** - Open chromatin regions
3. **chip_seq_tf_binding.png** - Transcription factor binding
4. **all_predictions_heatmap.png** - All data combined
5. **annotation_summary.png** - Gene overlap types
6. **biotype_comparison.png** - Predictions by gene type

---

## Complete Workflow Example: B Cell Analysis

```bash
# Step 1: Get reference files
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.refGene.gtf.gz
gunzip mm10.refGene.gtf.gz

# Step 2: Verify setup
python debug_gtf.py mm10.gtf --bed-file bcell_regions.bed

# Step 3: Run analysis for B cells (use intestine tissue = GALT)
python alphagenome_gtf_pipeline_nearest.py bcell_regions.bed mm10.gtf \
  --api-key YOUR_API_KEY \
  --tissues intestine \
  --buffer 5000 \
  --output bcell_predictions.csv

# Step 4: Generate plots
python plot_alphagenome_results.py bcell_predictions.csv \
  --output bcell_plots/ \
  --exclude-intergenic \
  --by-biotype

# Step 5: Check results
head -5 bcell_predictions.csv
ls -lh bcell_plots/
```

---

## What Each Plot Shows

### RNA-seq Distribution
Shows predicted expression levels for your regions:
- **Histogram**: Most genes have low expression (expected)
- **Top genes**: These should include known B cell genes (CD19, Ighg1, etc.)
- **Mean vs Max**: Variation in expression across region

### Chromatin Accessibility
Shows open/closed chromatin:
- **DNase/ATAC**: Most active regions accessible
- **Correlation with RNA**: Good = open chromatin → expressed
- **Top regions**: Most accessible sites in your data

### ChIP-seq TF Binding
Shows transcription factor (CTCF) binding:
- **Distribution**: Wide range of binding strengths
- **Top sites**: Strongest CTCF binding locations
- **TF vs Expression**: Relationship between binding and expression

### Heatmap
All prediction types together:
- **Pattern**: Coordinated activity of related genes
- **Clusters**: Genes with similar activity profiles
- **Completeness**: Shows which predictions are available

### Annotation Summary
Overview of gene annotations:
- **Bar chart**: Overlapping vs nearest genes
- **Distance plot**: How far nearest genes are
- **Interpretation**: Region-gene relationships

---

## Interpreting Your B Cell Results

### Expected Findings:
```
✓ Top expressed genes: CD19, Ighg1, BCR-related genes
✓ High chromatin accessibility at B cell genes
✓ CTCF binding at immunoglobulin loci
✓ Good correlation between accessibility and expression
✓ Overlapping genes match B cell biology
```

### What the Plots Should Show:
1. **RNA-seq**: Peak at B cell genes
2. **DNase**: Accessible at B cell promoters
3. **ChIP-seq**: CTCF at switching regions
4. **Heatmap**: B cell gene cluster
5. **Annotation**: Good overlapping genes

---

## Quality Control Checks

### Before Proceeding:
```bash
# 1. Verify GTF has genes
grep -c $'\tgene\t' mm10.gtf
# Should output: > 0

# 2. Check regions have annotations
grep -v "^chrom" results.csv | grep -c "no_gene"
# Should output: < 10% of total

# 3. Verify predictions present
head results.csv | cut -d',' -f1,7,10,13
# Should see numbers, not blanks

# 4. Check plot generation
ls -lh bcell_plots/
# Should see 6 PNG files
```

---

## Output Organization

```
project/
├── mm10.gtf                          # Reference
├── bcell_regions.bed                 # Input regions
├── bcell_predictions.csv             # Analysis results
│
└── bcell_plots/                      # Visualizations
    ├── rna_seq_distribution.png
    ├── chromatin_accessibility.png
    ├── chip_seq_tf_binding.png
    ├── all_predictions_heatmap.png
    ├── annotation_summary.png
    └── biotype_comparison.png
```

---

## Customization Options

### Different Tissues
```bash
# Brain tissue
python alphagenome_gtf_pipeline_nearest.py regions.bed mm10.gtf \
  --api-key KEY --tissues brain

# Lung tissue
python alphagenome_gtf_pipeline_nearest.py regions.bed mm10.gtf \
  --api-key KEY --tissues lung

# Multiple tissues (compare)
for tissue in brain lung intestine; do
  python alphagenome_gtf_pipeline_nearest.py regions.bed mm10.gtf \
    --api-key KEY --tissues $tissue --output results_$tissue.csv
  python plot_alphagenome_results.py results_$tissue.csv \
    --output plots_$tissue/
done
```

### Sequence Context
```bash
# Smaller context (faster)
--sequence-length SEQUENCE_LENGTH_100KB

# Larger context (more information)
--sequence-length SEQUENCE_LENGTH_500KB
```

### Plotting Options
```bash
# Only specific plots
python plot_alphagenome_results.py results.csv \
  --plots rna_seq dnase \
  --output plots/

# High resolution for publication
--dpi 600

# Exclude uninteresting regions
--exclude-intergenic

# Organize by biotype
--by-biotype

# Organize by annotation type
--by-annotation
```

---

## Troubleshooting

### No Genes Found (no_gene in results)
```bash
# Check GTF
python debug_gtf.py mm10.gtf --bed-file regions.bed

# Use buffer to find nearby genes
python alphagenome_gtf_pipeline_nearest.py regions.bed mm10.gtf \
  --api-key KEY --tissues brain --buffer 50000
```

### Plots Don't Show Data
```bash
# Verify CSV has prediction columns
head results.csv

# Check for data
grep -v "^chrom" results.csv | head | cut -d',' -f7,10,13

# If empty, re-run analysis with correct API key
```

### Too Many Intergenic Regions
```bash
# Likely reason: regions genuinely intergenic
# Use --buffer to find nearby genes
# Or exclude with --exclude-intergenic in plots

python plot_alphagenome_results.py results.csv \
  --output plots/ \
  --exclude-intergenic
```

---

## Next Steps for Publication

1. **Generate high-resolution plots**
   ```bash
   python plot_alphagenome_results.py results.csv \
     --output publication_figs/ --dpi 600 --by-biotype
   ```

2. **Create summary statistics**
   - Copy means/medians from plot output
   - Create supplementary table

3. **Validate findings**
   - Compare with experimental data
   - Check known gene patterns
   - Verify tissue appropriateness

4. **Prepare figures**
   - Insert PNG plots into presentation/paper
   - Add figure legends
   - Include supplementary data

5. **Document methodology**
   - Record commands used
   - Note GTF version/source
   - Document tissue selection rationale

---

## Summary Commands

```bash
# Complete analysis in one go
TISSUE="intestine"  # or brain, lung, colon, kidney

python alphagenome_gtf_pipeline_nearest.py regions.bed mm10.gtf \
  --api-key YOUR_API_KEY \
  --tissues $TISSUE \
  --buffer 10000 \
  --output results_$TISSUE.csv && \
python plot_alphagenome_results.py results_$TISSUE.csv \
  --output plots_$TISSUE/ \
  --by-biotype \
  --exclude-intergenic

# View results
head results_$TISSUE.csv
ls -lh plots_$TISSUE/
```

---

## File Reference

| File | Purpose |
|------|---------|
| alphagenome_gtf_pipeline_nearest.py | Main analysis with nearest genes |
| plot_alphagenome_results.py | Generate plots from results |
| debug_gtf.py | Troubleshoot GTF/BED issues |
| VISUALIZATION_GUIDE.md | How to interpret plots |
| TROUBLESHOOT_NO_GENE.md | Fix "no_gene" annotation |

---

## Final Checklist

- [ ] GTF file downloaded and decompressed
- [ ] BED file prepared (3 columns: chrom, start, end)
- [ ] GTF compatible with BED (chromosomes match)
- [ ] API key ready
- [ ] Correct tissue selected for your cell type
- [ ] Analysis complete (results.csv created)
- [ ] Plots generated successfully
- [ ] Summary statistics printed
- [ ] Results match expected biology
- [ ] Files ready for presentation/publication

---

**Ready to go! Your complete AlphaGenome analysis workflow is set up.** 🎉

For questions on specific plots, see **VISUALIZATION_GUIDE.md**
For issues with annotations, see **TROUBLESHOOT_NO_GENE.md**
For quick troubleshooting, see **QUICK_FIX_NO_GENE.md**
