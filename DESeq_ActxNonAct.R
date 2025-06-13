library(readxl)
#####count_table#####
countData <- read_excel("data/countData.xlsx")
# Set the first column as row names
countData <- as.data.frame(countData)
rownames(countData) <- countData[[1]]   
countData[[1]] <- NULL
####sample_table####
sampleTable=read_excel("SampleTable/table2_dds.xlsx")
sampleTable$Condition=as.factor(sampleTable$Condition)
sampleTable$CellType=as.factor(sampleTable$CellType)
sampleTable$Status=as.factor(sampleTable$Status)
sampleTable$Batch=as.factor(sampleTable$Batch)
sampleTable$Group <- as.factor(paste(sampleTable$Condition, sampleTable$Status, sampleTable$CellType, sep = "_"))
###DESeq2###
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = sampleTable,
                              design = ~ Group)

##dds2
dds2 <-DESeq(dds)
resultsNames(dds2)
res_1 <- results(dds2, contrast = c("Group", "HD_Activated_BrightMBCs", "HD_Nonactivated_BrightMBCs"))
res_2 <- results(dds2, contrast = c("Group", "HD_Activated_NaïveBCs", "HD_Nonactivated_NaïveBCs"))
res_3 <- results(dds2, contrast = c("Group", "CVID_Activated_BrightMBCs", "CVID_Nonactivated_BrightMBCs"))
res_4 <- results(dds2, contrast = c("Group", "CVID_Activated_NaïveBCs", "CVID_Nonactivated_NaïveBCs"))
###pheatmap package###
#transform in data frame
DEGs_1=data.frame(res_1[(res_1$padj)<0.05,]) #DEG_HDMBC.xlsx
DEGs_2=data.frame(res_2[(res_2$padj)<0.05,]) #DEG_HDNaive.xlsx
DEGs_3=data.frame(res_3[(res_3$padj)<0.05,]) #DEG_CVIDMBC.xlsx
DEGs_4=data.frame(res_4[(res_4$padj)<0.05,]) #DEG_CVIDNaive.xlsx

