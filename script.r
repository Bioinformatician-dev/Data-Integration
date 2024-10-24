# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(reshape2)

# Step 1: Load Data
# Load RNA-Seq data
rna_data <- read.csv("rna_seq_data.csv", row.names = 1)  # Rows are genes, columns are samples
# Load Proteomics data
protein_data <- read.csv("proteomics_data.csv", row.names = 1)  # Rows are proteins, columns are samples

# Optionally load Genomic data
# genomic_data <- read.csv("genomic_data.csv", row.names = 1)

# Step 2: Normalize RNA-Seq Data
# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = rna_data,
                              colData = data.frame(condition = rep(c("Control", "Treatment"), each=3)), 
                              design = ~ condition)

# Normalize the data
dds <- DESeq(dds)
rlog_data <- rlog(dds, blind = TRUE)

# Extract normalized counts
norm_counts <- assay(rlog_data)

# Step 3: Ensure Common Samples
# Identify common samples between RNA and protein data
common_samples <- intersect(colnames(norm_counts), colnames(protein_data))
norm_counts <- norm_counts[, common_samples]
protein_data <- protein_data[, common_samples]

# Step 4: Correlation Analysis
# Calculate correlation
cor_results <- cor(t(norm_counts), t(protein_data), method = "pearson")

# View correlation matrix
print(cor_results)

# Step 5: Visualize Correlation Heatmap
pheatmap(cor_results, 
         main = "Correlation between Gene Expression and Protein Levels",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 10,
         display_numbers = TRUE)

# Step 6: Save Results
write.csv(cor_results, file = "correlation_results.csv")

# Optional: Further analysis with genomic data (if available)
# Example: Gene set enrichment analysis based on significant genes from RNA-Seq
# Assume you have a list of significant genes based on some threshold
# significant_genes <- rownames(rna_data)[which(rowMeans(norm_counts) > threshold)]

# If genomic data is available, you could integrate it similarly:
# genomic_data <- genomic_data[, common_samples]
# Perform analysis to identify relationships

# Optional: Save normalized RNA-Seq and protein data
w
