# Load DESeq2 library
library(DESeq2)

# Read count matrix from StringTie output
countMatrix <- read.table("output_transcripts.gtf", header=FALSE, comment.char = "#")

# Sample information (modify this according to your sample conditions)
sampleInfo <- data.frame(condition = c("A", "A", "B", "B"))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countMatrix[, c(7, 8, 9, ...)], colData = sampleInfo, design = ~condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)

# Save results
write.csv(res, file = "differential_expression_results_DESeq2.csv")
