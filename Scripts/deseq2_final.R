# Check if 'BiocManager' is installed and install 'DESeq2'
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)  # For generating heatmaps

# Set working directory
setwd("C:/Users/singh/Desktop/DH307/E-MTAB-9950")

# Load count data
countData <- read.csv('gene_count_matrix.csv', header = TRUE, sep = ",")
head(countData)

# Load metadata
metaData <- read.csv('metadata.csv', header = TRUE, sep = ",")
head(metaData)

# Construct DESeqDataSet Object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              design = ~ dex, tidy = TRUE)

# Run DESeq function
dds <- DESeq(dds)

# Check result table
res <- results(dds)
head(res)

# Sort summary list with adjusted p-value
res <- res[order(res$padj),]
head(res)

# Set the output folder
output_folder <- "output_folder"
dir.create(output_folder, showWarnings = FALSE)  # Create output folder

# Save plots for all genes in a PDF file
pdf_file <- file.path(output_folder, "all_plots.pdf")
pdf(pdf_file)

# Original boxplots
par(mfrow=c(2, 3))
plotCounts(dds, gene = "ENSG00000224864.4|ENSG00000224864", intgroup = "dex")
plotCounts(dds, gene = "ENSG00000125740.15|FOSB", intgroup = "dex")
plotCounts(dds, gene = "ENSG00000239839.7|DEFA3", intgroup = "dex")
plotCounts(dds, gene = "ENSG00000123080.12|CDKN2C", intgroup = "dex")
plotCounts(dds, gene = "ENSG00000110203.10|FOLR3", intgroup = "dex")
plotCounts(dds, gene = "ENSG00000144668.12|ITGA9", intgroup = "dex")

# MA plot
plotMA(dds, ylim = c(-2, 2), main = "MA Plot")

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
# Extracting data for the 6 genes of interest
genes_of_interest <- c("ENSG00000224864.4|ENSG00000224864", "ENSG00000125740.15|FOSB",
                       "ENSG00000239839.7|DEFA3", "ENSG00000123080.12|CDKN2C",
                       "ENSG00000110203.10|FOLR3", "ENSG00000144668.12|ITGA9")

selected_genes_data <- counts(dds)[genes_of_interest,]

# Scale the data
heatmap_data_scaled <- scale(selected_genes_data)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Create the heatmap
pheatmap(heatmap_data_scaled, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, col = colors)

# PCA plot
plotPCA(vsd, intgroup = "dex")

dev.off()  # Close the PDF device

# Save results table
write.csv(res, file.path(output_folder, "deseq_results.csv"), row.names = FALSE)