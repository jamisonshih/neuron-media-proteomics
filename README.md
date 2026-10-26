# Neuronal Culture Media Shape Protein Abundance and Turnover in Human iPSC-Derived Neurons

R scripts for label-free and dSILAC protein analysis.

**Input:** Tables as `.xlsx`, `.csv` or `.tsv`. Sample columns must be numeric and end in a replicate number (`Ctrl1`, `Ctrl2`, `KO_1`, ...); the group is the column name without the number. The gene column must be `PG.Genes`.

| Script | Analysis |
| --- | --- |
| `01_pca.R` | PCA, log10 + per-protein z-score, complete cases |
| `02_sample_clustering.R` | Group Spearman correlation, replicate clustering |
| `03_limma.R` | All pairwise limma comparisons on log2(x + 1), BH; absolute log2FC > 1, padj < 0.05 |
| `04_go_ora.R` | GO over-representation (hypergeometric, BH) |
| `05_heatmap.R` | Z-scored heatmap of the listed genes |
| `06_half_life.R` | Peptide H/L ratios to protein half-lives |
| `07_2d_go.R` | 2D GO enrichment: rank MANOVA (Wilks), BH across BP/CC/MF |
