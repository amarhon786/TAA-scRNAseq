# TAA-scRNAseq

Sex differences in human ascending thoracic aortic aneurysm (TAA), using the public single-cell dataset GSE155468 ([Li et al. 2020, Circulation](https://pubmed.ncbi.nlm.nih.gov/33017217/)). The question is whether the disease looks different in men and women at the gene level, so the DE model has a sex:disease interaction term instead of only testing TAA vs control.

The R scripts go from the raw 10x matrices to figures: Seurat v5 QC, Harmony integration, marker-score annotation into 8 cell types, pseudobulk DESeq2 per cell type, then fgsea / clusterProfiler enrichment.

## Data

GSE155468, 10x scRNA-seq of human ascending aorta, ages 56-78.

| Group | Male | Female |
|---|---|---|
| TAA | 4 | 4 |
| Control | 1 | 2 |

## Method

Counts are summed per sample and cell type (pseudobulk), so the unit of replication is the patient (n = 11) rather than the ~50k cells. Per cell type DESeq2 fits

```
~ sex + disease + sex:disease
```

with Female and Control as reference levels. Five contrasts come out of it:

1. TAA vs Control averaged over sex (`disease + 0.5 * sex:disease`)
2. Male TAA vs Male Control (4 vs 1, low power)
3. Female TAA vs Female Control (4 vs 2)
4. Male TAA vs Female TAA (4 vs 4)
5. The sex:disease interaction, i.e. genes whose disease response differs by sex

Enrichment is pre-ranked GSEA on the Wald statistic against MSigDB Hallmark, KEGG and GO:BP, plus ORA on the significant genes. FDR is Benjamini-Hochberg.

## Findings

1. In contractile SMCs the top TAA vs control genes are immune / interferon genes, not SMC genes: CD69 (log2FC +5.0, padj 3.0e-6), IFNG (+6.4, padj 7e-6), IFITM1, CCL5 and CD52. Modulated SMCs add CTHRC1 (an aortic fibrosis marker, padj 3.5e-6) and the chemokines CCL4, CXCR4, CCL5 and CCL2.
2. Female TAA macrophages have a much stronger interferon / oxidative response than male ones. Hallmark GSEA for male TAA vs female TAA (negative NES = higher in females):

   | Pathway | NES | padj |
   |---|---:|---:|
   | INTERFERON_GAMMA_RESPONSE | -2.82 | 3.6e-24 |
   | OXIDATIVE_PHOSPHORYLATION | -2.67 | 4.8e-20 |
   | ALLOGRAFT_REJECTION | -2.45 | 1.5e-14 |
   | MYC_TARGETS_V1 | -2.43 | 1.5e-14 |
   | MTORC1_SIGNALING | -2.34 | 2.8e-12 |

3. The same female-skewed pattern shows up in SMCs (contractile: IFN-gamma NES -2.83, padj 2e-23 and IFN-alpha -2.93, padj 2e-18; modulated: IFN-gamma -2.94, padj 2e-25, IL6-JAK-STAT3 -2.44, EMT -2.16) and as higher OXPHOS in female NK, T and plasma cells. Women generally mount stronger IFN responses, but I haven't found this pattern described in TAA across several aortic lineages before.
4. Endothelial cells go the other way, with males showing the higher inflammatory signal: TNFA_SIGNALING_VIA_NFKB (NES +2.23, padj 9.7e-10), P53_PATHWAY (+1.87, 4.4e-5) and ANDROGEN_RESPONSE (+1.94, 2.4e-4). The androgen hit is also a decent check that the sex labels are right.

Less certain:

5. Interaction genes in plasma cells include IGHG2 (log2FC +11.1, padj 0.006), HSPA6 (-12.1, 0.006), PTGDS (-15.8, 0.011) and NFKB1 (+4.2, 0.013), and macrophages have IGLC3 (-21.7) and IGHG1 (-9.6) near the FDR cutoff. Those effect sizes are extreme for a cohort this small, so it could be one patient driving it.
6. MT1X in NK cells is the best-powered interaction gene in the dataset (log2FC -3.22, padj 3.8e-6).

## Caveats

- The cohort is small. Male TAA vs Male Control is 4 vs 1 and gives 21 significant genes across all 8 cell types combined, so treat it as exploratory. The interaction term has the same single-male-control problem.
- In modulated SMCs the interaction didn't reach FDR (GRIA2 / CDKN2A / KRT14 look divergent but padj >= 0.44), so those are suggestive only.
- It's one cohort. It would need replication in an independent TAA scRNA-seq dataset (e.g. Pedroza et al. 2022/2023).

## Running it

Needs R >= 4.4. From the project root:

```r
source("R/00_setup.R")           # renv + packages, ~20 min the first time
source("R/01_download.R")        # GEO download, ~5 min, ~3 GB
source("R/02_qc_integrate.R")    # QC, Harmony, clustering, annotation, ~20-40 min
source("R/03_celltype_abundance.R")
source("R/04_pseudobulk_de.R")
source("R/05_enrichment.R")
source("R/06_figures.R")
```

Each script reads what the one before it wrote, so any step can be rerun on its own. `results/figures/` has the two SMC figures (5 panels, PDF) with their legends in the matching `.md` files and `results/tables/` has the DEG workbooks for the two SMC types. The per-contrast CSVs (~60 files) and the raw data are gitignored and regenerate from the scripts.

## Citation

Li Y, Ren P, Dawson A, Vasquez HG, Ageedi W, Zhang C, Luo W, Chen R, Li Y, Kim S, Lu HS, Cassis LA, Coselli JS, Daugherty A, Shen YH, LeMaire SA. Single-Cell Transcriptome Analysis Reveals Dynamic Cell Populations and Differential Gene Expression Patterns in Control and Aneurysmal Human Aortic Tissue. Circulation 2020. PMID: 33017217.
