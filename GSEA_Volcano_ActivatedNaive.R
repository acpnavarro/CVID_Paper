#Annotation
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#create a gene list from a df resulted from DESeq2 
df <-read_excel("data/df.AN_CVIDxHD.xlsx")
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

data_go <- as.data.frame(gse)

##plots
p <- as.data.frame(gse) %>%
  mutate(
    num_genes = sapply(strsplit(core_enrichment, "/"), length),
    GeneRatio = num_genes / setSize) %>%
  dplyr::mutate(type = dplyr::if_else(NES > 0, "upregulated", "downregulated")) %>%
  dplyr::arrange(NES) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:20) %>%
  ggplot(., aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  facet_grid(~ type) +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("GO: Activated Naïve CVID vs. HD")


p <- as.data.frame(gse) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:10) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("GO: Activated Naïve CVID vs. HD")

##Single GO plot -  somatic recombination
# Extract NES and adjusted p-value
nes_value <- round(gse@result[gse@result$ID == "GO:0002204", "NES"], 2)
fdr_value <- signif(gse@result[gse@result$ID == "GO:0002204", "p.adjust"], 2)

# Create a plot title with both values
plot_title <- paste0("Somatic Recombination of Immunoglobulin Genes Involved in Immune Response\n(Activated Naïve CVID vs. HD)\nNES = ", 
                     nes_value, ", FDR = ", fdr_value)

# Generate the GSEA plot
somatic <- gseaplot2(gse,
                     geneSetID = "GO:0002204",
                     title = plot_title)

##Single GO plot - isotype switching
# Extract NES and adjusted p-value
nes_value <- round(gse@result[gse@result$ID == "GO:0045190", "NES"], 2)
fdr_value <- signif(gse@result[gse@result$ID == "GO:0045190", "p.adjust"], 2)

# Create a plot title with both values
plot_title <- paste0("Isotype Switching\n(Activated Naïve CVID vs. HD)\nNES = ", 
                     nes_value, ", FDR = ", fdr_value)

# Generate the GSEA plot
isotype<- gseaplot2(gse,
                    geneSetID = "GO:0045190",
                    title = plot_title)

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

data_react <- as.data.frame(react)

p <- as.data.frame(react) %>%
  dplyr::arrange(p.adjust) %>% 
  dplyr::group_by(sign(NES)) %>% 
  dplyr::slice(1:15) %>%
  ggplot(., aes(x = NES, y = fct_reorder(Description, NES),fill=p.adjust)) + 
  geom_bar(stat="identity") +
  theme_bw(base_size = 10) +
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  ylab(NULL) +
  ggtitle("Reactome: Activated Naïve CVID vs. HD")

##Single Reactome plot - Senescence
# Extract NES and adjusted p-value
nes_value <- round(react@result[react@result$ID == "R-HSA-2559582", "NES"], 2)
fdr_value <- signif(react@result[react@result$ID == "R-HSA-2559582", "p.adjust"], 2)

# Create a plot title with both values
plot_title <- paste0("Senescence-Associated Secretory Phenotype (SASP)\n(Activated Naïve CVID vs. HD)\nNES = ", 
                     nes_value, ", FDR = ", fdr_value)

# Generate the GSEA plot
sene <- gseaplot2(react,
                  geneSetID = "R-HSA-2559582",
                  title = plot_title)

##Hallmarks
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
  ggtitle("Hallmark Pathways: Activated Naïve CVID vs. HD")


stat5<- gseaplot2(em4, geneSetID = "HALLMARK_IL2_STAT5_SIGNALING", 
                  title = "IL2 STAT5 Signaling")

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

#####Volcano_DEGs
###volcano plot###
df <-read_excel("data/df.AN_CVIDxHD.xlsx")
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
df$diffexpressed[df$log2FoldChange > 1.5 & df$padj< 0.05] <- "UP"
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
  coord_cartesian(ylim = c(0, 15), xlim = c(-5, 10)) +
  labs(color = 'Gene expression',
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  theme(axis.title = element_text(size = 30,
                                  color = "black",
                                  face = "bold"))+
  theme(legend.title = element_text(color = "black", size = 15), legend.text = element_text(color = "black", size=10))+
  scale_x_continuous(breaks = seq(-10, 10, 5))+
  labs(title="Activated Naïve B cells: CVID vs. HD")+
  theme(plot.title = element_text(size = 15))+
  geom_text_repel(aes(label = delabel), size = 4, max.overlaps = Inf)

ggsave("volcanoAN_CVIDxHD.png", dpi = 300, width = 10, height = 8)

