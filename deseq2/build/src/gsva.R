# Load necessary libraries
suppressPackageStartupMessages(library("GSVA"))
suppressPackageStartupMessages(library("GSEABase"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("limma"))

# Define the command-line options
option_list <- list(
    make_option(c("--data"), type = "character", default = NULL,
                help = "Path to the normalized count matrix .csv file", metavar = "character"),
    make_option(c("--gene_set"), type = "character", default = NULL,
                help = "Path to a .gmt file with a list of gene sets to perform GSVA with", metavar = "character"),
    make_option(c("--filter"), type = "character", default = NULL,
                help = "comma-separated list of gene set names to filter the final result by", metavar = "character"),
    make_option(c("--min"), type = "numeric", default = -1.0,
                help = "Minimum value for the heatmap color scale", metavar = "numeric"),
    make_option(c("--mid"), type = "numeric", default = 0.0,
                help = "Middle value for the heatmap color scale", metavar = "numeric"),
    make_option(c("--max"), type = "numeric", default = 1.0,
                help = "Maximum value for the heatmap color scale", metavar = "numeric"),
    make_option(c("--output"), type = "character", default = NULL,
                help = "Output file path for the heatmap plot. Should end in .png", metavar = "character"),
    make_option(c("--method"), type = "character", default = "gsva",
                help = "Algorithm to use to calculate pathway activity scores. Options: gsva, ssgsea", metavar = "character")
    
)

# Parse the command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Input files
gmt_file <- opt$gene_set
deseq2_results_file <- opt$data
output_png <- opt$output
gsva_activity_min <- opt$min
gsva_activity_mid <- opt$mid
gsva_activity_max <- opt$max
method <- opt$method

# check that the input files are provided and have the correct extensions
if (is.null(deseq2_results_file) || !grepl(".csv$", deseq2_results_file)) {
  stop("The DESeq2 result file must be provided using the --data flag and must end with .csv")
}

if (is.null(output_png) || !grepl(".png$", output_png)) {
  stop("The output file path must be provided using the --output flag and must end with .png")
}

# Load .gmt file if it's provided
gene_sets <- NULL
if (!is.null(gmt_file) && grepl(".gmt$", gmt_file)) {
    gene_sets <- getGmt(gmt_file)
} else {
    # If no gene set file is provided, stop
    stop("The gene set file must be provided using the --gene_set flag and must end with .gmt")
}

# Load DESeq2 results (genes as row names, normalized expression values)
deseq2_results <- read.csv(deseq2_results_file, header = TRUE)

# Set row names to the "gene" column and remove it from data frame
if ("gene" %in% colnames(deseq2_results)) {
  rownames(deseq2_results) <- deseq2_results$gene
  deseq2_results <- deseq2_results[, !colnames(deseq2_results) %in% "gene"]
} else {
  stop("ERROR: The 'gene' column was not found in the header.")
}

# Remove non-expression columns
columns_to_remove <- c("baseMean", "lfcSE", "log2FoldChange", "stat", "pvalue", "padj")
deseq2_results <- deseq2_results[, !colnames(deseq2_results) %in% columns_to_remove]

# Convert to matrix
expression_matrix <- as.matrix(deseq2_results)

# remove all columns that don't contain 'Pan24'
expression_matrix <- expression_matrix[, grepl("Pan24", colnames(expression_matrix))]

# remove rows that have fewer than 3 samples with expression > 10
expression_matrix <- expression_matrix[rowSums(expression_matrix > 10) >= 3, ]

gsvaPar <- NULL
if (method == "gsva") {
    # round expression matrix to integers
    expression_matrix <- round(expression_matrix)

    # Create a GSVA parameter object
    gsvaPar <- gsvaParam(expression_matrix, gene_sets, minSize=10, maxSize=500, tau=1.0, kcdf="Poisson")
} else if (method == "ssgsea") {
    # log transform
    expression_matrix <- log(expression_matrix + 1, base = 2)

    # subtract row means
    expression_matrix <- expression_matrix - rowMeans(expression_matrix)

    # Create a ssGSEA parameter object
    gsvaPar <- ssgseaParam(expression_matrix, gene_sets, minSize=10, maxSize=500, alpha=0.25)
} else {
    stop("ERROR: The method must be either 'gsva' or 'ssgsea'")
}

# write to .csv
# file_name <- gsub(".png", "_expression_matrix.csv", output_png)
# write.csv(expression_matrix, file_name, row.names = TRUE)

# Calculate pathway activity scores
gsva_results <- gsva(gsvaPar)

# Filter the results geneset name if specified
filtered_gsva_results <- NULL
if (!is.null(opt$filter)) {
    filter_genesets <- unlist(strsplit(opt$filter, ","))
    filtered_gsva_results <- gsva_results[rownames(gsva_results) %in% filter_genesets, , drop = FALSE]
} else {
    filtered_gsva_results <- gsva_results
}

# create .csv file with GSVA results
# csv_name <- gsub(".png", ".csv", output_png)
# write.csv(gsva_results, file = csv_name, row.names = TRUE)

# print(method, 'done. generating heatmap...')

# Create a color function for the heatmap
col_fun <- colorRamp2(c(gsva_activity_min, gsva_activity_mid, gsva_activity_max), c("blue", "white", "red"), space = "LUV")

# if not filtering, cluster by rows
cluster_rows <- is.null(opt$filter)

# Generate a heatmap
heatmap <- Heatmap(
    filtered_gsva_results,
    name = "Pathway Activity",
    col = col_fun,
    cluster_rows = cluster_rows,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(
        title = "Activity",
        at = c(gsva_activity_min, gsva_activity_mid, gsva_activity_max),    # Set the legend ticks
        labels = c(as.character(gsva_activity_min), as.character(gsva_activity_mid), as.character(gsva_activity_max))  # Define custom labels
    )
)

# save as .png
png(output_png, width = 2400, height = 1600, res = 300)
draw(heatmap, heatmap_legend_side = "right")
dev.off()

# cat("Pathway activity heatmap saved to", output_png, "\n")

# save as .svg
svg_filename <- gsub(".png", ".svg", output_png)
ggsave(svg_filename, plot = grid.grabExpr(draw(heatmap)), width = 8, height = 5.33)








# make a heatmap with the significant pathways annotated


# define contrasts for limma analysis
contrasts <- list(
  c("Pan24Control", "Pan24VIP"),
  c("Pan24Control", "Pan24Zen"),
  c("Pan24Control", "Pan24Combo")
#   c("Pan38Control", "Pan38VIP"),
#   c("Pan38Control", "Pan38Zen"),
#   c("Pan38Control", "Pan38Combo")
)

# Get sample names from GSVA results
sample_names <- colnames(gsva_results)

# Store extracted GSVA data for each contrast
contrast_results <- list()

for (contrast in contrasts) {
    # get control and treatment sample names
    control_samples <- grep(contrast[1], sample_names, value = TRUE)
    treatment_samples <- grep(contrast[2], sample_names, value = TRUE)

    if (length(control_samples) == 0 || length(treatment_samples) == 0) {
        cat("No samples found for contrast", contrast[1], "vs", contrast[2], "\n")
        next
    }

    # get GSVA results for control and treatment samples
    gsva_subset <- gsva_results[, c(control_samples, treatment_samples), drop = FALSE]

    # make design matrix for limma
    condition_factor <- factor(c(
        rep("control", length(control_samples)),
        rep("treatment", length(treatment_samples))
    ))
    design <- model.matrix(~condition_factor)

    # perform statistical testing on the GSVA activity scores for this contrast using limma
    fit <- lmFit(gsva_subset, design)
    fit <- eBayes(fit)
    contrast_table <- topTable(fit, coef = 2, adjust.method = "fdr", number = Inf)  # coef = 2 corresponds to "treatment" vs "control"

    # Store results
    contrast_results[[paste(contrast[1], "vs", contrast[2], sep = "_")]] <- contrast_table

    # Write results to a CSV
    contrast_table <- data.frame(GeneSet = rownames(contrast_table), contrast_table, row.names = NULL)
    file_name <- gsub(".png", paste(contrast[1], "vs", contrast[2], ".csv", sep = "_"), output_png)
    write.csv(contrast_table, file_name, row.names = FALSE)
}


# Extract all unique conditions from the contrasts file
# condition_names <- unique(unlist(contrasts))

# reorder as Pan24Control, Pan38Control, Pan24VIP, Pan38VIP, Pan24Zen, Pan38Zen, Pan24Combo, Pan38Combo
# condition_names <- c("Pan24Control", "Pan38Control", "Pan24VIP", "Pan38VIP", "Pan24Zen", "Pan38Zen", "Pan24Combo", "Pan38Combo")
condition_names <- c("Pan24Control", "Pan24VIP", "Pan24Zen", "Pan24Combo")
# condition_names <- c("Pan38Control", "Pan38VIP", "Pan38Zen", "Pan38Combo")

# Initialize a matrix to store averaged GSVA scores
avg_gsva_matrix <- matrix(NA, nrow = nrow(gsva_results), ncol = length(condition_names))
rownames(avg_gsva_matrix) <- rownames(gsva_results)
colnames(avg_gsva_matrix) <- condition_names

# Compute the mean GSVA score for each condition
for (condition in condition_names) {
  matching_samples <- grep(condition, sample_names, value = TRUE)  # Find sample names containing condition substring
  avg_gsva_matrix[, condition] <- rowMeans(gsva_results[, matching_samples, drop = FALSE], na.rm = TRUE)
}

# subtract the mean of each row
# avg_gsva_matrix <- avg_gsva_matrix - rowMeans(avg_gsva_matrix)

# write avg_gsva_matrix to csv
# file_name <- gsub(".png", "_avg_gsva_matrix.csv", output_png)
# write.csv(avg_gsva_matrix, file_name, row.names = TRUE)


# Identify significant gene sets from all contrasts
# significant_gene_sets <- unique(unlist(lapply(contrast_results, function(contrast_table) {
#     rownames(contrast_table[contrast_table$adj.P.Val <= 0.05, ])
# })))

# top_n <- 30

# take top 30 downregulated gene sets for Pan24Combo vs Pan24Control
# top30_pan24_combo <- contrast_results[["Pan24Control_vs_Pan24Combo"]]
# top30_pan24_combo <- top30_pan24_combo[top30_pan24_combo$logFC < 0, ]
# top30_pan24_combo <- top30_pan24_combo[order(top30_pan24_combo$adj.P.Val), ]
# top30_pan24_combo <- rownames(top30_pan24_combo[1:top_n, ])

# print('top 30 gene sets for Pan24Combo vs Pan24Control:')
# print(top30_pan24_combo)

# take top 30 downregulated gene sets for Pan38Combo vs Pan38Control
# top30_pan38_combo <- contrast_results[["Pan38Control_vs_Pan38Combo"]]
# top30_pan38_combo <- top30_pan38_combo[top30_pan38_combo$logFC < 0, ]
# top30_pan38_combo <- top30_pan38_combo[order(top30_pan38_combo$adj.P.Val), ]
# top30_pan38_combo <- rownames(top30_pan38_combo[1:top_n, ])

# print('top 30 gene sets for Pan38Combo vs Pan38Control:')
# print(top30_pan38_combo)


# get intersection
# significant_gene_sets <- intersect(top30_pan24_combo, top30_pan38_combo)


# Filter the results geneset name if specified
filtered_gsva_matrix <- NULL
if (!is.null(opt$filter)) {
    significant_gene_sets <- unlist(strsplit(opt$filter, ","))
    filtered_gsva_matrix <- avg_gsva_matrix[significant_gene_sets, , drop = FALSE]
} else {
    filtered_gsva_matrix <- avg_gsva_matrix
}



# print('significant gene sets:')
# print(significant_gene_sets)

# Subset the averaged GSVA matrix to only include significant pathways
# filtered_gsva_matrix <- avg_gsva_matrix[significant_gene_sets, , drop = FALSE]

# write to csv
# file_name <- gsub(".png", "avg_gsva_matrix_filtered.csv", output_png)
# write.csv(filtered_gsva_matrix, file_name, row.names = TRUE)

# Check the dimensions after filtering
cat("Number of significant gene sets:", nrow(filtered_gsva_matrix), "\n")




# heatmap
# Create an annotation matrix for significance
default_pvalue <- 1.0
significance_matrix <- matrix(default_pvalue, nrow = nrow(filtered_gsva_matrix), ncol = ncol(filtered_gsva_matrix))
rownames(significance_matrix) <- rownames(filtered_gsva_matrix)
colnames(significance_matrix) <- colnames(filtered_gsva_matrix)

# Loop through contrasts and mark significant conditions with an asterisk
for (contrast in contrasts) {
    # get the gsva result from contrast_results
    contrast_name <- paste(contrast[1], "vs", contrast[2], sep = "_")
    contrast_table <- contrast_results[[contrast_name]]
    control <- contrast[1]
    treatment <- contrast[2]
    
    # extract adjusted p-values for this contrast
    gene_sets <- rownames(contrast_table)
    df <- data.frame(GeneSet = rownames(contrast_table), contrast_table, row.names = NULL)

    # save padj in the significance matrix
    for (gene_set in gene_sets) {
        if (gene_set %in% rownames(significance_matrix)) {
            padj <- df[df$GeneSet == gene_set, "adj.P.Val"]
            significance_matrix[gene_set, treatment] <- padj
        }
    }
}

# Define function to add asterisks for significant values
asterisk_fun <- function(j, i, x, y, width, height, fill) {
    # padj <- significance_matrix[i, j]
    # annot <- ""
    # if (padj < 0.01) {
    #     annot <- "**"
    # } else if (padj < 0.05) {
    #     annot <- "*"
    # }

    # if (annot != "") {   
    #     # y-center the asterisk
    #     # y <- y - height / 3

    #     grid.text(annot, x, y, gp = gpar(fontsize = 10))
    # }
}

# Generate the clustered heatmap with significance annotations
heatmap <- Heatmap(
    filtered_gsva_matrix,
    name = "Pathway Activity",
    col = col_fun,
    cluster_rows = cluster_rows,  # Cluster pathways
    cluster_columns = FALSE,  # Keep conditions in the provided order
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(
        title = "Activity",
        at = c(gsva_activity_min, gsva_activity_mid, gsva_activity_max),  # Set legend ticks
        labels = c(as.character(gsva_activity_min), as.character(gsva_activity_mid), as.character(gsva_activity_max))
    ),
    cell_fun = asterisk_fun  # Add asterisks for significance
)

# write to .png
output_png_name <- gsub(".png", "_gsva_heatmap.png", output_png)
png(output_png_name, width = 2400, height = 1600, res = 300)
draw(heatmap, heatmap_legend_side = "right")
dev.off()