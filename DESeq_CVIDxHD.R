###############################################################################
# DESeq_CVIDxHD.R
#
# Differential expression analysis: CVID vs. healthy donors (HD), run
# separately for each cell type and activation state.
#
#   IN = resting (non-activated) naïve B cells
#   AN = activated naïve B cells
#   IB = resting (non-activated) CD27bright memory B cells
#   AB = activated CD27bright memory B cells
#
# Outputs (written to data/ and results/):
#   data/df.<code>_CVIDxHD.xlsx      full DESeq2 results (input to GSEA_Volcano_*.R)
#   data/DEG_<code>_padj0.05.xlsx    DEGs with padj < 0.05 (input to heatmap_*.R)
#   results/<code>_top50.xlsx        top 25 up- and 25 downregulated DEGs
#                                    (|log2FC| > 1.5, ranked by |log2FC|)
#
# DEGs are defined as Benjamini-Hochberg padj < 0.05, without a log2FC
# threshold, as described in the manuscript.
###############################################################################

library(DESeq2)
library(dplyr)
library(readxl)
library(openxlsx)

dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Input data
# ---------------------------------------------------------------------------
countData <- read.xlsx("data/countData.xlsx", rowNames = TRUE)

sampleTable <- read_excel("SampleTable/table2_dds.xlsx")
sampleTable <- sampleTable %>%
  mutate(across(c(Condition, CellType, Status, Batch), as.factor),
         Group = factor(paste(Condition, Status, CellType, sep = "_")))

# ---------------------------------------------------------------------------
# 2. DESeq2 model (fitted once and reused for all contrasts)
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData   = sampleTable,
                              design    = ~ Group)
dds <- DESeq(dds)

# ---------------------------------------------------------------------------
# 3. Contrasts: CVID vs. HD within each cell type and state
# ---------------------------------------------------------------------------
contrasts <- list(
  IN = c("Group", "CVID_Nonactivated_NaïveBCs",   "HD_Nonactivated_NaïveBCs"),
  AN = c("Group", "CVID_Activated_NaïveBCs",      "HD_Activated_NaïveBCs"),
  IB = c("Group", "CVID_Nonactivated_BrightMBCs", "HD_Nonactivated_BrightMBCs"),
  AB = c("Group", "CVID_Activated_BrightMBCs",    "HD_Activated_BrightMBCs")
)

padj_cutoff  <- 0.05
lfc_cutoff   <- 1.5   # used only to select top genes for display
n_top_each   <- 25    # top up- and downregulated genes in the top-50 table

# Top up- and downregulated DEGs ranked by absolute log2FC
# (base R subsetting is used so that gene names, stored as row names, are kept)
top_degs <- function(degs, n) {
  up   <- degs[degs$log2FoldChange >  lfc_cutoff, ]
  down <- degs[degs$log2FoldChange < -lfc_cutoff, ]
  rbind(head(up[order(-up$log2FoldChange), ], n),
        head(down[order(down$log2FoldChange), ], n))
}

# ---------------------------------------------------------------------------
# 4. Run all contrasts and export results
# ---------------------------------------------------------------------------
summary_table <- lapply(names(contrasts), function(code) {

  res <- as.data.frame(results(dds, contrast = contrasts[[code]]))

  # Full results table
  write.xlsx(res, file.path("data", paste0("df.", code, "_CVIDxHD.xlsx")),
             asTable = TRUE, rowNames = TRUE)

  # Significant DEGs (genes with NA padj, e.g. filtered by independent
  # filtering, are excluded)
  degs <- res[!is.na(res$padj) & res$padj < padj_cutoff, ]
  write.xlsx(degs, file.path("data", paste0("DEG_", code, "_padj0.05.xlsx")),
             asTable = TRUE, rowNames = TRUE)

  # Top 50 DEGs for display
  write.xlsx(top_degs(degs, n_top_each),
             file.path("results", paste0(code, "_top50.xlsx")),
             asTable = TRUE, rowNames = TRUE)

  data.frame(contrast = code,
             DEGs     = nrow(degs),
             up       = sum(degs$log2FoldChange > 0),
             down     = sum(degs$log2FoldChange < 0))
})

summary_table <- bind_rows(summary_table)
print(summary_table)
write.xlsx(summary_table, "results/DEG_summary_CVIDxHD.xlsx", asTable = TRUE)

# ---------------------------------------------------------------------------
# 5. Session information for reproducibility
# ---------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "results/sessionInfo_DESeq_CVIDxHD.txt")
