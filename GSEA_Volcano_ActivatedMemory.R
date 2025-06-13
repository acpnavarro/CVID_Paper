#Annotation
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#create a gene list from a df resulted from DESeq2 
df <-read_excel("data/df.AB_CVIDxHD.xlsx")
        
# log2 fold change 
original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$row.names
##gene_list_entrez
#Convert gene IDs 
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
ids <-as.data.frame(ids)
# remove duplicate IDS 
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$row.names %in% dedup_ids$SYMBOL,]
# Create a new column in df2 with the corresponding ENTREZ IDs
df2$ENTREZ = dedup_ids$ENTREZID
# Create a vector of the gene unuiverse
entrez_geneList <- df2$log2FoldChange
# Name vector with ENTREZ ids
names(entrez_geneList) <- df2$ENTREZ
# omit any NA values 
entrez_geneList<-na.omit(entrez_geneList)
# sort the list in decreasing order (required for clusterProfiler)
entrez_geneList = sort(entrez_geneList, decreasing = TRUE)
head(entrez_geneList)
###############################################################################
##gseGO
set.seed(123)
gse <- gseGO(geneList=entrez_geneList, 
             ont ="BP", 
             keyType = "ENTREZID",
             minGSSize = 10, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH",
             eps =0,
             seed=TRUE,
             exponent =1,
             by="fgsea")

p <- as.data.frame(gse) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("GO: Activated Memory CVID vs. HD")

##REACTOME
set.seed(123)
react <- gsePathway(geneList = entrez_geneList,
                    organism = "human",
                    minGSSize = 10,
                    maxGSSize  = 500,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    verbose=TRUE,
                    seed=TRUE,
                    by="fgsea")

p <- as.data.frame(react) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:15) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("Reactome: Activated Memory CVID vs. HD")

##hallmarks
hall <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(hall)
em4<- GSEA(entrez_geneList, TERM2GENE = hall)
head(em4)
data_em4 <- as.data.frame(em4)

p_hall <- as.data.frame(em4) %>%
  dplyr::filter(!is.na(NES)) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::group_by(sign(NES)) %>%
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES), fill = p.adjust)) + 
  geom_bar(stat = "identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low = 'red', high = 'blue', guide = guide_colorbar(reverse = TRUE)) + 
  ylab(NULL) +
  ggtitle("Hallmark Pathways: Activated Memory CVID vs. HD")


#Motifs
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame
#GSEA motifs
C3_t2g <- msigdbr(species = "Homo sapiens", category = "C3") %>% 
  dplyr::select(gs_name, entrez_gene)
head(C3_t2g)
em2 <- GSEA(entrez_geneList, TERM2GENE = C3_t2g)
head(em2)
data_em2 <- as.data.frame(em2)

# Extract NES and FDR (adjusted p-value) for STAT3
stat3_data <- data_em2 %>% filter(ID == "STAT3_02")
NES_stat3 <- round(stat3_data$NES, 2)
FDR_stat3 <- signif(stat3_data$p.adjust, 2)

# Create the gseaplot2 with NES and FDR in the title
stat3 <- gseaplot2(em2,
                   geneSetID = "STAT3_02",
                   title = paste0("STAT3 Motif (Activated Memory CVID vs. HD)\nNES = ", NES_stat3, ", FDR = ", FDR_stat3))

#####Volcano_DEGs
###volcano plot###
df <-read_excel("data/df.AB_CVIDxHD.xlsx")
# Set the first column as row names
df <- as.data.frame(df)
rownames(df) <- df[[1]]   
df[[1]] <- NULL 

# Theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(0.5), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(0.5), color = 'black'),
              plot.title = element_text(hjust = 0.5, size = 15),
              axis.title = element_text(size = 30, color = "black", face = "bold"),
              legend.title = element_text(color = "black", size = 15),
              legend.text = element_text(color = "black", size=10)
            ))

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
df$diffexpressed <- "NO"
# if log2Foldchange > 1.5 and padjvalue < 0.05, set as "UP"
df$diffexpressed[df$log2FoldChange > 1.5 & df$padj < 0.05] <- "UP"
# if log2Foldchange < -1.5 and padjvalue < 0.05, set as "DOWN"
df$diffexpressed[df$log2FoldChange < -1.5 & df$padj < 0.05] <- "DOWN"
# Explore a bit
head(subset(df, diffexpressed == "DOWN")[order(subset(df, diffexpressed == "DOWN")$padj), ])
# create a column for symbol in df 
df$symbol <- rownames(df)
# create a column delabel in df and select the top 10 symbols
# Top 5 UP-regulated genes by log2FC
top_up <- head(df[df$diffexpressed == "UP", ][order(-df$log2FoldChange[df$diffexpressed == "UP"]), "symbol"], 5)

# Top 5 DOWN-regulated genes by log2FC
top_down <- head(df[df$diffexpressed == "DOWN", ][order(df$log2FoldChange[df$diffexpressed == "DOWN"]), "symbol"], 5)

df$delabel <- ifelse(df$symbol %in% c(top_up, top_down), df$symbol, NA)

ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(aes(color=df$diffexpressed), alpha=0.8, size = 2, shape=20 ) +
  scale_color_manual(values = c("#00AFBB", "grey", "#e2062c"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  coord_cartesian(ylim = c(0, 10), xlim = c(-5, 10)) +
  labs(color = 'Gene expression',
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  theme(axis.title = element_text(size = 30,
                                  color = "black",
                                  face = "bold"))+
  theme(legend.title = element_text(color = "black", size = 15), legend.text = element_text(color = "black", size=10))+
  scale_x_continuous(breaks = seq(-10, 10, 5))+
  labs(title="Activated Memory B cells: CVID vs. HD")+
  theme(plot.title = element_text(size = 15))+
  geom_text_repel(aes(label = delabel), size = 4, max.overlaps = Inf)

ggsave("volcanoAB_CVIDxHD.png", dpi = 300, width = 10, height = 8)

