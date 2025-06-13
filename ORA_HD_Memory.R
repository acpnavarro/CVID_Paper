# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(msigdbr)
library(openxlsx) 
library(enrichplot)
library(ggplot2)
library(patchwork)

# Unique genes list
unique_genes <- read_excel("C:/Users/Anna/OneDrive/Documents/CVID_project/Paper/Venn/unique_genes_HD_MBC.xlsx")
colnames(unique_genes) <- "SYMBOL"

#Convert unique gene list to ENTREZ IDs
unique_entrez <- bitr(unique_genes$SYMBOL, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db)

table(duplicated(unique_entrez$SYMBOL))
duplicated_symbols <- unique_entrez$SYMBOL[duplicated(unique_entrez$SYMBOL)]
unique_entrez[unique_entrez$SYMBOL %in% duplicated_symbols, ]
unique_entrez_filtered <- unique_entrez[!duplicated(unique_entrez$SYMBOL), ]
# Run ORA
ego <- enrichGO(gene         = unique_entrez_filtered$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                minGSSize    = 5,
                readable     = TRUE)



#KEGG
gene_list <- unique_entrez_filtered$ENTREZID
kegg_ora <- enrichKEGG(gene         = gene_list,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)

# Convert to readable gene symbols
kegg_ora <- setReadable(kegg_ora, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

#Reactome
reactome_ora <- enrichPathway(gene         = gene_list,
                              organism     = "human",
                              pvalueCutoff = 0.05,
                              readable     = TRUE)
#Hallmark
gene_list <- unique(as.character(unique_entrez_filtered$ENTREZID))
# Get Hallmark gene sets
library(msigdbr)
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
term2gene <- hallmark_genesets[, c("gs_name", "entrez_gene")]

# Run enrichment
hallmark_ora <- enricher(
  gene          = gene_list,
  TERM2GENE     = term2gene,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Save results to Excel files
write.xlsx(as.data.frame(ego), file = "ORA_GO_HD_M.xlsx", rowNames = FALSE)
write.xlsx(as.data.frame(kegg_ora), file = "ORA_KEGG_HD_M.xlsx", rowNames = FALSE)
write.xlsx(as.data.frame(reactome_ora), file = "ORA_Reactome_HD_M.xlsx", rowNames = FALSE)
write.xlsx(as.data.frame(hallmark_ora), file = "ORA_Hallmark_HD_M.xlsx", rowNames = FALSE)

# Create individual plots
go_barplot <- barplot(ego, showCategory = 10, title = "GO Biological Process HD_MBC")
hallmark_dotplot <- dotplot(hallmark_ora, showCategory = 10, title = "Hallmark Pathways HD_MBC")

# Combine plots side by side
combined_plot <- go_barplot + hallmark_dotplot + plot_layout(ncol = 2)

# Show combined plot
print(combined_plot)
ggsave("combined_ORA_GO_Hallmark_HD_MBC.png", plot = combined_plot, width = 14, height = 6, dpi = 300)

# Create individual plots
go_barplot <- barplot(ego, showCategory = 10, title = "GO Biological Process HD_MBC")
reactome_dotplot <- dotplot(reactome_ora, showCategory = 10, title = "Reactome Enrichment HD_MBC")

# Combine plots side by side
combined_plot <- go_barplot + reactome_dotplot + plot_layout(ncol = 2)

# Show combined plot
print(combined_plot)
ggsave("combined_ORA_GO_Reactome_HD_MBC.png", plot = combined_plot, width = 14, height = 6, dpi = 300)
