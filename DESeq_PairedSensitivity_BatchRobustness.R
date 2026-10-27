# =====================================================================
# DESeq_PairedSensitivity_BatchRobustness.R
#
# Paired sensitivity analysis and batch-robustness analysis of the
# CVID RNA-seq dataset (robustness checks of the primary DGE results).
#
# Outputs (written to ./sensitivity_robustness/):
#
# 1) Paired sensitivity analysis
#    DESeq2 with design ~ Patient + Status for the 5 CVID patients with
#    matched resting/activated samples per cell type, compared against
#    the primary unpaired activation contrasts (~ Group, full cohort).
#
# 2) Batch-robustness analysis
#    Disease contrasts (CVID vs HD within each activation state and
#    cell type) re-analysed with ~ Batch + Condition and compared with
#    the primary ~ Group model (published DEGs).
#
# Requires: library.R (loads readxl, DESeq2)
# Input:    data/countData.xlsx, SampleTable/table2_dds.xlsx
#
# Note: robustness analyses were run with DESeq2 v1.46.0 (R 4.4.3).
# The primary analysis used DESeq2 v1.42.1 (R 4.3.2); results are
# consistent across versions.
# =====================================================================

library(readxl)
library(DESeq2)

dir.create("sensitivity_robustness", showWarnings = FALSE)

##### count_table #####
countData <- read_excel("data/countData.xlsx")
countData <- as.data.frame(countData)
rownames(countData) <- countData[[1]]
countData[[1]] <- NULL

##### sample_table #####
sampleTable <- as.data.frame(read_excel("SampleTable/table2_dds.xlsx"))
sampleTable$Condition <- as.factor(sampleTable$Condition)
sampleTable$CellType  <- as.factor(sampleTable$CellType)
sampleTable$Status    <- as.factor(sampleTable$Status)
sampleTable$Batch     <- as.factor(sampleTable$Batch)

# Patient ID is embedded in the sample name (e.g. "AB.CVID01" -> "CVID01")
sampleTable$Patient <- as.factor(trimws(sub("^[A-Z]+[.]", "", sampleTable$Sample)))

# Harmonise sample names (remove stray whitespace, if any) and match
# the column order of the count matrix
sampleTable$Sample <- trimws(sampleTable$Sample)
rownames(sampleTable) <- sampleTable$Sample
rownames(countData) <- trimws(rownames(countData))
countData <- countData[, rownames(sampleTable)]
stopifnot(all(colnames(countData) == rownames(sampleTable)))

sampleTable$Group <- as.factor(paste(sampleTable$Condition, sampleTable$Status,
                                     sampleTable$CellType, sep = "_"))

##### Primary model (full cohort, ~ Group) #####
# This is the model used for the published DEGs (Figures 2-4).
dds_full <- DESeq(DESeqDataSetFromMatrix(countData = countData,
                                         colData   = sampleTable,
                                         design    = ~ Group))

# All primary contrasts (activation within each group + disease within
# each state x cell type)
prim <- list(
  # Activation contrasts (CVID)
  CVID_Naive_act  = results(dds_full, contrast = c("Group", "CVID_Activated_Na\u00efveBCs",
                                                   "CVID_Nonactivated_Na\u00efveBCs")),
  CVID_MBC_act    = results(dds_full, contrast = c("Group", "CVID_Activated_BrightMBCs",
                                                   "CVID_Nonactivated_BrightMBCs")),
  # Disease contrasts (CVID vs HD)
  RestingNaive    = results(dds_full, contrast = c("Group", "CVID_Nonactivated_Na\u00efveBCs",
                                                   "HD_Nonactivated_Na\u00efveBCs")),
  ActivatedNaive  = results(dds_full, contrast = c("Group", "CVID_Activated_Na\u00efveBCs",
                                                   "HD_Activated_Na\u00efveBCs")),
  RestingMBC      = results(dds_full, contrast = c("Group", "CVID_Nonactivated_BrightMBCs",
                                                   "HD_Nonactivated_BrightMBCs")),
  ActivatedMBC    = results(dds_full, contrast = c("Group", "CVID_Activated_BrightMBCs",
                                                   "HD_Activated_BrightMBCs"))
)

# =====================================================================
# PART 1 - PAIRED SENSITIVITY ANALYSIS
# =====================================================================
# 5 CVID patients have both resting and activated samples per cell type:
#   Naive: CVID08, CVID10, CVID13, CVID14, CVID21
#   MBC:   CVID01, CVID10, CVID13, CVID14, CVID21
# Design: ~ Patient + Status (activation effect estimated within patient)

paired_patients <- list(
  Naive = c("CVID08", "CVID10", "CVID13", "CVID14", "CVID21"),
  MBC   = c("CVID01", "CVID10", "CVID13", "CVID14", "CVID21")
)
celltype_of <- c(Naive = "Na\u00efveBCs", MBC = "BrightMBCs")
prim_act_of <- c(Naive = "CVID_Naive_act", MBC = "CVID_MBC_act")

for (ct in names(paired_patients)) {

  pats <- paired_patients[[ct]]
  sel  <- sampleTable$Patient %in% pats &
          sampleTable$CellType == celltype_of[[ct]] &
          sampleTable$Condition == "CVID"
  st_paired <- sampleTable[sel, ]
  ct_paired <- countData[, rownames(st_paired)]

  dds_paired <- DESeq(DESeqDataSetFromMatrix(countData = ct_paired,
                                             colData   = st_paired,
                                             design    = ~ Patient + Status))
  res_paired <- results(dds_paired, contrast = c("Status", "Activated", "Nonactivated"))

  # ---- Compare paired vs primary unpaired ----
  df_p <- as.data.frame(res_paired)
  df_u <- as.data.frame(prim[[prim_act_of[[ct]]]])
  genes <- intersect(rownames(df_p), rownames(df_u))
  df_p <- df_p[genes, ]; df_u <- df_u[genes, ]

  cmp <- data.frame(
    gene            = genes,
    log2FC_paired   = df_p$log2FoldChange,
    padj_paired     = df_p$padj,
    log2FC_unpaired = df_u$log2FoldChange,
    padj_unpaired   = df_u$padj,
    stringsAsFactors = FALSE
  )
  cmp$Sig_paired   <- !is.na(cmp$padj_paired)   & cmp$padj_paired   < 0.05
  cmp$Sig_unpaired <- !is.na(cmp$padj_unpaired) & cmp$padj_unpaired < 0.05
  cmp$Direction_same <- sign(cmp$log2FC_paired) == sign(cmp$log2FC_unpaired)

  write.csv(cmp,
            file.path("sensitivity_robustness",
                      paste0("Paired_vs_Unpaired_", ct, ".csv")),
            row.names = FALSE)

  # Key genes (STAT5 isoforms) reported in the manuscript
  key <- cmp[cmp$gene %in% c("STAT5A", "STAT5B"), ]
  print(paste0("Paired sensitivity (", ct, "): ",
               sum(cmp$Sig_paired), " DEGs paired vs ",
               sum(cmp$Sig_unpaired), " unpaired; direction concordance ",
               round(100 * mean(cmp$Direction_same, na.rm = TRUE), 1), "%"))
  print(key)
}

# =====================================================================
# PART 2 - BATCH-ROBUSTNESS ANALYSIS
# =====================================================================
# Disease contrasts (CVID vs HD) within each activation state x cell
# type. Batch is NOT confounded with Condition (HD and CVID samples
# occur in both the resting and the activated batches), so
# ~ Batch + Condition is estimable for disease contrasts. Batch IS
# completely confounded with Status (batches 1-3 = all resting,
# batches 4-5 = all activated), so activation contrasts cannot be
# batch-adjusted.

contrasts4 <- list(
  RestingNaive   = list(celltype = "Na\u00efveBCs",  status = "Nonactivated"),
  ActivatedNaive = list(celltype = "Na\u00efveBCs",  status = "Activated"),
  RestingMBC     = list(celltype = "BrightMBCs",     status = "Nonactivated"),
  ActivatedMBC   = list(celltype = "BrightMBCs",     status = "Activated")
)

batch_summary <- data.frame()
for (nm in names(contrasts4)) {

  cc <- contrasts4[[nm]]
  sel <- sampleTable$CellType == cc$celltype & sampleTable$Status == cc$status
  st_sub <- sampleTable[sel, ]
  ct_sub <- countData[, rownames(st_sub)]

  # Plain model (no batch term), subsetted data
  dds_plain <- DESeq(DESeqDataSetFromMatrix(countData = ct_sub,
                                            colData   = st_sub,
                                            design    = ~ Condition))
  res_plain <- results(dds_plain, contrast = c("Condition", "CVID", "HD"))

  # Batch-adjusted model
  dds_batch <- DESeq(DESeqDataSetFromMatrix(countData = ct_sub,
                                            colData   = st_sub,
                                            design    = ~ Batch + Condition))
  res_batch <- results(dds_batch, contrast = c("Condition", "CVID", "HD"))

  # ---- Compare with the primary (~ Group) model ----
  df_pub <- as.data.frame(prim[[nm]])
  df_bt  <- as.data.frame(res_batch)
  df_pl  <- as.data.frame(res_plain)
  genes  <- Reduce(intersect, list(rownames(df_pub), rownames(df_bt), rownames(df_pl)))
  df_pub <- df_pub[genes, ]; df_bt <- df_bt[genes, ]; df_pl <- df_pl[genes, ]

  out <- data.frame(
    gene          = genes,
    log2FC_pub    = df_pub$log2FoldChange,
    padj_pub      = df_pub$padj,
    log2FC_plain  = df_pl$log2FoldChange,
    padj_plain    = df_pl$padj,
    log2FC_batch  = df_bt$log2FoldChange,
    padj_batch    = df_bt$padj,
    stringsAsFactors = FALSE
  )
  out$sig_pub   <- !is.na(out$padj_pub)   & out$padj_pub   < 0.05
  out$sig_plain <- !is.na(out$padj_plain) & out$padj_plain < 0.05
  out$sig_batch <- !is.na(out$padj_batch) & out$padj_batch < 0.05
  out$pub_still_sig <- out$sig_pub & out$sig_batch
  out$dir_concordant <- sign(out$log2FC_pub) == sign(out$log2FC_batch)

  write.csv(out,
            file.path("sensitivity_robustness",
                      paste0("BatchRobustness_", nm, ".csv")),
            row.names = FALSE)

  batch_summary <- rbind(batch_summary, data.frame(
    contrast                 = nm,
    n_samples                = ncol(ct_sub),
    n_batches                = nlevels(droplevels(st_sub$Batch)),
    published_DEGs           = sum(out$sig_pub),
    sig_plain                = sum(out$sig_plain),
    sig_batch                = sum(out$sig_batch),
    pub_still_sig            = sum(out$pub_still_sig),
    pub_still_sig_pct        = round(100 * sum(out$pub_still_sig) / max(sum(out$sig_pub), 1), 1),
    direction_concordance_pct = round(100 * mean(out$dir_concordant[out$sig_pub], na.rm = TRUE), 1)
  ))
}

print(batch_summary)
write.csv(batch_summary,
          file.path("sensitivity_robustness", "BatchRobustness_Summary.csv"),
          row.names = FALSE)

# =====================================================================
# SESSION INFO
# =====================================================================
sink(file.path("sensitivity_robustness", "sessionInfo.txt"))
print(sessionInfo())
sink()
