# CVID_Paper

Analysis code for the RNA-seq study:

> **Pousas Navarro, A.\***, Pour Akaber, S.\*, Waldenström, J.\*, Friman, V.\*, Zanchetta, A., Basílio, J., Barresi, S., Hambrecht, B., Tartaglia, M., Lundqvist, A., Malmström, N., Karlsson, J., Tilevik, A., Bemark, M., Bohle, B., & Grimsholm, O. (2026). Impaired activation of STAT5 pathways in the naïve B cell compartment in patients with common variable immunodeficiency harbouring CD27bright memory B cells. *Frontiers in Immunology*. Accepted for publication. [DOI to be added]
>
> \*These authors contributed equally to this work.

## Overview

This repository contains the R scripts used for the transcriptomic analyses in the study above. We profiled sorted naïve B cells and CD27bright memory B cells (MBCs) from patients with common variable immunodeficiency (CVID) and healthy donors (HDs), at steady state and after in vitro activation. Patients were selected for harbouring a clearly detectable CD27bright MBC population.

## Study design

| Factor | Levels |
|---|---|
| Condition | CVID, HD |
| Status | Activated, Non-activated (resting) |
| Cell type | Naïve B cells (CD24+CD27-CD38+), CD27bright MBCs (CD24+CD27++) |

- **Samples:** 17 patients with CVID and 9 HDs, plus resting naïve and CD27bright MBCs from 4 additional HDs from a previously published dataset (65 RNA-seq samples in total).
- **Activation of naïve B cells:** 1 μg/mL rhCD40L + 10 μg/mL anti-IgM/IgG/IgA + 20 ng/mL rhIL-21, 48 hours.
- **Activation of CD27bright MBCs:** 0.35 μM CpG-B ODN2006 (TLR9 agonist), 48 hours.

Abbreviations used in file and sample names:

| Code | Meaning |
|---|---|
| IN | Resting (non-activated) naïve B cells |
| AN | Activated naïve B cells |
| IB | Resting (non-activated) CD27bright MBCs |
| AB | Activated CD27bright MBCs |
| MBC | Memory B cells |
| HD | Healthy donor |

## Upstream processing (not included as scripts)

Reads were processed before the analyses in this repository as follows:

- Quality control: FastQC (v0.12.1); adapter trimming: Trimmomatic (v0.39); summary: MultiQC (v1.17)
- Alignment to GRCh38 with STAR (v2.6.1a) in 2-pass mode
- Gene quantification with featureCounts (Subread v2.0.6), Ensembl 110 annotation
- Ensembl IDs mapped to gene symbols, keeping protein-coding genes only; low-expression genes (mean log₂ CPM < 1) and immunoglobulin and MHC genes were excluded

## Scripts in this repository

| Script | Purpose | Figure |
|---|---|---|
| `library.R` | Installs and loads all required packages | |
| `PCA_VSTdata.R` | Variance-stabilising transformation and PCA of all samples | 2A |
| `DESeq_ActxNonAct.R` | DESeq2: activated vs. resting cells within each group (CVID and HD; naïve and CD27bright MBCs) | 2B, 2E |
| `Venn_UniqueGenes.R` | Overlap of activation-induced DEGs between CVID and HD; extraction of group-specific genes | 2B, 2E |
| `ORA_HD_Naive.R`, `ORA_CVID_Naive.R` | Over-representation analysis (GO, KEGG, Reactome, Hallmark) of group-specific DEGs, naïve B cells | 2C, 2D |
| `ORA_HD_Memory.R`, `ORA_CVID_Memory.R` | Over-representation analysis of group-specific DEGs, CD27bright MBCs | 2F, 2G |
| `heatmap_RestingNaive.R`, `heatmap_ActivatedNaive.R` | Heatmaps of the top 50 up- and 50 downregulated genes, CVID vs. HD (z-scored VST values) | 3A, 3B |
| `heatmap_RestingMemory.R`, `heatmap_ActivatedMemory.R` | As above, CD27bright MBCs | 4A, 4B |
| `GSEA_Volcano_RestingNaive.R`, `GSEA_Volcano_ActivatedNaive.R` | Volcano plots and GSEA (GO, Reactome, Hallmark, C3 motifs including STAT5), CVID vs. HD, naïve B cells | 3C–3H |
| `GSEA_Volcano_RestingMemory.R`, `GSEA_Volcano_ActivatedMemory.R` | As above, CD27bright MBCs (including TLR9 cascade and STAT3 motif) | 4C–4H |

Differentially expressed genes were defined as Benjamini-Hochberg adjusted p < 0.05, without a log₂ fold-change threshold. Volcano plots highlight genes with |log₂FC| > 1.5.

## Suggested order of execution

1. `library.R`
2. `PCA_VSTdata.R`
3. `DESeq_ActxNonAct.R`
4. `Venn_UniqueGenes.R`
5. `ORA_*.R`
6. `GSEA_Volcano_*.R`
7. `heatmap_*.R`

## Software

R (v4.3.2) with DESeq2 (v1.42.1), clusterProfiler (v4.10.1), ReactomePA (v1.46.0), msigdbr (v10.0.1), GSVA (v2.2.0), org.Hs.eg.db (v3.18.0), pheatmap (v1.0.12), ggplot2 (v3.5.1) and ggvenn (v0.1.10). The full list of packages is in `library.R`.

## Input data

The scripts expect the following inputs:

- `data/countData.xlsx`: gene-level count matrix
- `SampleTable/table2_dds.xlsx`: sample metadata (condition, status, cell type, batch)

Raw sequencing data are publicly available under accession numbers **GSE307095** and **GSE227377**.

## Contact

Anna Pousas Navarro, MD, MSc
[LinkedIn](https://www.linkedin.com/in/anna-carolina-pousas-navarro-73949339/)

Corresponding author of the study: Ola Grimsholm, Medical University of Vienna.

## License

[Choose a license, e.g. MIT, or state "All rights reserved".]
