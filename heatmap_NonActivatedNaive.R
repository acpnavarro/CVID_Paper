library(readxl)
library(pheatmap)
library(ggplot2)
library(gridExtra) 

#Read
IN <- read_excel("data/DEG_IN_padj0.05.xlsx")
IN <- as.data.frame(IN)
#####sample_table####
sampleTable=read_excel("SampleTable/table2_dds.xlsx")
sampleTable$Condition=as.factor(sampleTable$Condition)
sampleTable$CellType=as.factor(sampleTable$CellType)
sampleTable$Status=as.factor(sampleTable$Status)
sampleTable$Batch=as.factor(sampleTable$Batch)
sampleTable$Group <- as.factor(paste(sampleTable$Condition, sampleTable$Status, sampleTable$CellType, sep = "_"))
#Selecting the 100 most up and down regulated
library(dplyr)
# Top 50 most downregulated genes (log2FoldChange < 0)
top50_down <- IN %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>%  # most negative first
  dplyr::slice(1:50)

# Top 50 most upregulated genes (log2FoldChange > 0)
top50_up <- IN %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%  # most positive first
  dplyr::slice(1:50)

top100<- bind_rows(top50_up, top50_down)
# VST data
vst_dds <- readRDS("C:/Users/Anna/OneDrive/Documents/CVID_project/norm_dds.RDS")

##prepare the matrix
# Extract gene names
top_genes <- top100$row.names
# Subset VST matrix
mat <- assay(vst_dds)[top_genes, ]
#Z score
s_mat=t(scale(t(mat)))
# select the samples
s_mat=as.data.frame(s_mat)
s_mat2<-s_mat[, grep("^IN.CVID|^IN.HD", names(s_mat))]

##prepare annotation sample table
ann_sample<-data.frame(sampleTable$Sample, sampleTable$Group)
colnames(ann_sample)<-c("Sample","Group")
df_ann=ann_sample %>% 
  filter(str_detect(Sample, "^IN.CVID|^IN.HD"))
df_ann <- data.frame(df_ann, row.names = 1)
df_ann$Group<-gsub("*_Nonactivated_NaïveBCs", "", df_ann$Group)

# Match annotation to heatmap matrix
df_ann <- df_ann[colnames(s_mat2), , drop = FALSE]

# Annotation colors
color_ann <- list(Group = c("CVID" = "#2a2a2a", "HD" = "#bdbdbd"))

# Heatmap
library(RColorBrewer)

# Custom color palette for smoother transitions
my_palette <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100)

# Enhanced heatmap
s_mat2 <- as.matrix(s_mat2)
attr(s_mat2, "name") <- NULL
pheatmap(s_mat2,
         name = "Z-score",  
         color = my_palette,
         scale = "none",
         clustering_method = "average",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         cutree_rows = 3,
         cutree_cols = 2,
         show_rownames = FALSE,
         annotation_col = df_ann,
         annotation_colors = color_ann,
         legend_breaks = c(-2, 0, 2),
         legend_labels = c("-2", "0", "2"),
         angle_col = "45",
         show_colnames = FALSE,
         fontsize = 10,
         fontsize_col = 8,
         border_color = NA,
         main = "Resting Naïve B Cells: CVID vs. HD")

dev.off()
