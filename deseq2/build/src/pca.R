# Load necessary libraries
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))

# Input files
deseq2_results_file <- "/w5home/bmoore/sequencing/MSpurgeon/RNAseq/03.DEseq/006_deseq_output_filtered/normalized_counts_filtered.csv" # Replace with your DESeq2 results file path
output_png <- "/w5home/bmoore/sequencing/MSpurgeon/RNAseq/03.DEseq/008_plots/PCA_norm_count"

# Load DESeq2 results
# Assuming DESeq2 results have genes as row names and normalized expression values
deseq2_results <- read.csv(deseq2_results_file, row.names = 1)

# Extract the expression matrix
expression_matrix <- as.matrix(deseq2_results)

# remove all samples that don't contain "Pan24"
#expression_matrix <- expression_matrix[, grepl("Pan24", colnames(expression_matrix))]

# remove rows that have fewer than 3 samples with expression > 10
expression_matrix <- expression_matrix[rowSums(expression_matrix > 10) >= 3, ]

# print head of expression matrix
# print(head(expression_matrix))
# stop()

# log2-transform the expression matrix (log2(x + 10))
expression_matrix <- log2(expression_matrix + 10)

# Perform PCA on the expression matrix
cat("Performing PCA on expression matrix...\n")
pca <- prcomp(t(expression_matrix), scale. = TRUE)

# Create a data frame for PCA results
pca_data <- as.data.frame(pca$x)
pca_data$Sample <- rownames(pca_data)

# Extract cell line information from sample names
pca_data$Sampletype <- case_when(
    grepl("4_", pca_data$Sample) ~ "NIKS_pLXSB_sandwich",
    grepl("5_", pca_data$Sample) ~ "NIKS_ER_sandwich",
    grepl("6_", pca_data$Sample) ~ "MCC_Raft",
    grepl("7_", pca_data$Sample) ~ "NHDF_pLXSB_sandwich",
    grepl("8_", pca_data$Sample) ~ "NHDF_ER_sandwich",
    TRUE ~ "Unknown"  # Fallback for unmatched names
)
print(pca_data$Sample)
pca_data$Rep <- case_when(
    grepl("A", pca_data$Sample) ~ "A",
    grepl("B", pca_data$Sample) ~ "B",
    grepl("C", pca_data$Sample) ~ "C",
    grepl("D", pca_data$Sample) ~ "D",
    grepl("E", pca_data$Sample) ~ "E",
    grepl("F", pca_data$Sample) ~ "F",
    grepl("G", pca_data$Sample) ~ "G",
    grepl("H", pca_data$Sample) ~ "H",
    grepl("I", pca_data$Sample) ~ "I",
    grepl("J", pca_data$Sample) ~ "J",
    TRUE ~ "Unknown"  # Fallback for unmatched names
)
head(pca_data$Rep)
# Plot the PCA results
cat("Generating PCA plot...\n")
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Sampletype, shape=Sampletype)) +
    geom_point() +
    stat_ellipse(geom = "polygon", aes(fill = Sampletype), alpha=0.2, level=0.6) +
    geom_text(aes(label = Rep), vjust = -1, size = 3) +
    labs(title = "PCA Plot", x = "PC1", y = "PC2") +
    theme_classic()

# Save PCA plot to PNG
ggsave(paste0(output_png,"_sampletype.png"), plot = pca_plot, width = 8, height = 6, dpi = 300)

cat("PCA plot saved to", output_png, "\n")