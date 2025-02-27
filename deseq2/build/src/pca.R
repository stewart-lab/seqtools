# Load necessary libraries
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))

# Input files
deseq2_results_file <- "/Users/rmillikin/Desktop/2025-01-10-RNASeq/006_deseq_output_filtered/normalized_counts_filtered.csv" # Replace with your DESeq2 results file path
output_png <- "/Users/rmillikin/Desktop/2025-01-10-RNASeq/009_plots/PCA_Pan24.png"

# Load DESeq2 results
# Assuming DESeq2 results have genes as row names and normalized expression values
deseq2_results <- read.csv(deseq2_results_file, row.names = 1)

# Extract the expression matrix
expression_matrix <- as.matrix(deseq2_results)

# remove all samples that don't contain "Pan24"
expression_matrix <- expression_matrix[, grepl("Pan24", colnames(expression_matrix))]

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
pca_data$CellLine <- case_when(
    grepl("Pan24", pca_data$Sample) ~ "Pan24",
    grepl("Pan33", pca_data$Sample) ~ "Pan33",
    grepl("Pan38", pca_data$Sample) ~ "Pan38",
    TRUE ~ "Unknown"  # Fallback for unmatched names
)

pca_data$Treatment <- case_when(
    grepl("Control", pca_data$Sample) ~ "Control",
    grepl("VIP", pca_data$Sample) ~ "VIP",
    grepl("Zen", pca_data$Sample) ~ "Zen",
    grepl("Combo", pca_data$Sample) ~ "Combo",
    TRUE ~ "Unknown"  # Fallback for unmatched names
)

pca_data$Rep <- case_when(
    grepl("Rep1", pca_data$Sample) ~ "Rep1",
    grepl("Rep2", pca_data$Sample) ~ "Rep2",
    grepl("Rep3", pca_data$Sample) ~ "Rep3",
    grepl("Rep4", pca_data$Sample) ~ "Rep4",
    TRUE ~ "Unknown"  # Fallback for unmatched names
)

# Plot the PCA results
cat("Generating PCA plot...\n")
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment, shape=Treatment)) +
    geom_point() +
    stat_ellipse(geom = "polygon", aes(fill = Treatment), alpha=0.2, level=0.6) +
    geom_text(aes(label = Rep), vjust = -1, size = 1.2) +
    labs(title = "PCA Plot", x = "PC1", y = "PC2") +
    theme_classic()

# Save PCA plot to PNG
ggsave(output_png, plot = pca_plot, width = 8, height = 6, dpi = 300)

cat("PCA plot saved to", output_png, "\n")