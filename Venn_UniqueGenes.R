library(ggvenn)
library(ggplot2)
library(openxlsx)

#Read DEGs Activated x Non-Activated
CVID_MBC= read_excel("data/DEG_CVIDMBC.xlsx")
HD_MBC= read_excel("data/DEG_HDMBC.xlsx")
CVID_Naive= read_excel("data/DEG_CVIDNaive.xlsx")
HD_Naive= read_excel("data/DEG_HDNaive.xlsx")

# Convert gene columns to character vectors
genes <- list(
  CVID_MBC   = as.character(CVID_MBC[[1]]),
  HD_MBC     = as.character(HD_MBC[[1]]),
  CVID_Naive = as.character(CVID_Naive[[1]]),
  HD_Naive   = as.character(HD_Naive[[1]])
)

#Venn
#Naive BCs
# Create the list of genes
genes_two <- list(
  CVID_Naive = genes$CVID_Naive,
  HD_Naive = genes$HD_Naive
)

# Create the Venn diagram with ggvenn
venn_plot <- ggvenn(
  genes_two,
  fill_color = c("#E69F00", "#56B4E9"),  
  stroke_size = 0.6,
  set_name_size = 5,
  text_size = 5
)

# Add title and annotation
final_plot <- venn_plot + 
  ggtitle("Differentially expressed genes: Naïve B Cells (Activated vs. Non-Activated)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic")
  ) +
  labs(caption = "CVID_Naïve = CVID Naïve B cells (Activated vs. Non-Activated)\nHD_Naïve = Healthy Donor Naïve B cells (Activated vs. Non-Activated)")

# Save the plot (change path if needed)
ggsave("venn_naive_bcells.png", plot = final_plot, width = 8, height = 6, dpi = 300)

#VENN NAIVE WITHOUT LEGENDS
# Create the Venn diagram with ggvenn
venn_plot <- ggvenn(
  genes_two,
  fill_color = c("#E69F00", "#56B4E9"), 
  stroke_size = 0.6,
  set_name_size = 5,
  text_size = 5
)

# Add title and annotation
final_plot <- venn_plot + 
  ggtitle("Differentially expressed genes: Naïve B Cells (Activated vs. Non-Activated)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic")
  ) 
# Save the plot 
ggsave("venn_naive_bcells_2.png", plot = final_plot, width = 8, height = 6, dpi = 300)

# MBC
# Create the list of genes
genes_3 <- list(
  CVID_MBC = genes$CVID_MBC,
  HD_MBC = genes$HD_MBC
)

# Create the Venn diagram with ggvenn
venn_plot <- ggvenn(
  genes_3,
  fill_color = c("#00cdc0", "#c3b6cc"),    stroke_size = 0.6,
  set_name_size = 5,
  text_size = 5
)

# Add a title and annotation using ggplot2

final_plot2 <- venn_plot + 
  ggtitle("Differentially expressed genes: Memory B Cells (Activated vs. Non-Activated)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic")
  ) +
  labs(caption = "CVID_MBC = CVID Bright Memory B cells (Activated vs. Non-Activated)\nHD_MBC = Healthy Donor Bright Memory B cells (Activated vs. Non-Activated)")
# Save the plot (change path if needed)
ggsave("venn_memory_bcells.png", plot = final_plot2, width = 9, height = 6, dpi = 300)

#VENN MBC WITHOUT LEGEND
# Create the Venn diagram with ggvenn
venn_plot <- ggvenn(
  genes_3,
  fill_color = c("#00cdc0", "#c3b6cc"),    stroke_size = 0.6,
  set_name_size = 5,
  text_size = 5
)

# Add a title and annotation using ggplot2

final_plot2 <- venn_plot + 
  ggtitle("Differentially expressed genes: Memory B Cells (Activated vs. Non-Activated)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.caption = element_text(hjust = 0, size = 10, face = "italic")
  ) 
# Save the plot (change path if needed)
ggsave("venn_memory_bcells_2.png", plot = final_plot2, width = 9, height = 6, dpi = 300)

#####################################################

# Extract unique genes for each group
unique_CVID_Naive <- setdiff(genes$CVID_Naive, genes$HD_Naive)
unique_HD_Naive <- setdiff(genes$HD_Naive, genes$CVID_Naive)
unique_CVID_MBC <- setdiff(genes$CVID_MBC,genes$HD_MBC)
unique_HD_MBC <- setdiff(genes$HD_MBC,genes$CVID_MBC)

# Convert vectors to data frames (with column name "Gene")
df_unique_CVID_Naive <- data.frame(Gene = unique_CVID_Naive)
df_unique_HD_Naive   <- data.frame(Gene = unique_HD_Naive)
df_unique_CVID_MBC   <- data.frame(Gene = unique_CVID_MBC)
df_unique_HD_MBC     <- data.frame(Gene = unique_HD_MBC)

# Save each as a separate Excel file
write.xlsx(df_unique_CVID_Naive, "unique_genes_CVID_Naive.xlsx", rowNames = FALSE)
write.xlsx(df_unique_HD_Naive,   "unique_genes_HD_Naive.xlsx", rowNames = FALSE)
write.xlsx(df_unique_CVID_MBC,   "unique_genes_CVID_MBC.xlsx", rowNames = FALSE)
write.xlsx(df_unique_HD_MBC,     "unique_genes_HD_MBC.xlsx", rowNames = FALSE)

