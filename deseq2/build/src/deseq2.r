suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("optparse"))

# Define the command-line options
option_list <- list(
    make_option(c("--counts"), type = "character", default = NULL,
                help = "Path to the counts CSV file", metavar = "character"),
    make_option(c("--design"), type = "character", default = NULL,
                help = "Path to the design CSV file", metavar = "character"),
    make_option(c("--contrasts"), type = "character", default = NULL,
                help = "Path to the contrasts CSV file", metavar = "character"),
    make_option(c("--output_dir"), type = "character", default = NULL,
                help = "Path to the output directory", metavar = "character")
)


# Parse the command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that required arguments are provided
if (is.null(opt$counts)) {
  stop("The counts file must be provided using the --counts flag")
}

if (is.null(opt$design)) {
  stop("The design file must be provided using the --design flag")
}

if (is.null(opt$contrasts)) {
  stop("The contrasts file must be provided using the --contrasts flag")
}

if (is.null(opt$output_dir)) {
  stop("The output directory must be provided using the --output-dir flag")
}

# Function to check if a file exists and has the .csv extension
check_csv_file <- function(file_path, file_type) {
  if (!file.exists(file_path)) {
    stop(paste("The", file_type, "file does not exist:", file_path))
  }
  if (tools::file_ext(file_path) != "csv") {
    stop(paste("The", file_type, "file must have a .csv extension:", file_path))
  }
}

# Check that the counts and design files exist and are CSV files
check_csv_file(opt$counts, "counts")
check_csv_file(opt$design, "design")
check_csv_file(opt$contrasts, "contrasts")

# remove trailing slash from output_dir if it exists
opt$output_dir <- gsub("/$", "", opt$output_dir)

# make the output directory
dir.create(opt$output_dir, showWarnings = FALSE)


# The design file is in CSV format and must have the headers "sample" and "condition"
# The contrasts file is in CSV format and must have the headers "treatment" and "control"


# Read the design file.
design <- read.csv2(opt$design, header=TRUE, stringsAsFactors=F, sep=",")
design$condition = factor(design$condition)

# Read the contrasts file.
contrasts <- read.csv2(opt$contrasts, header=TRUE, sep=",")

# Isolate the sample names.
sample_names <- design$sample

# read the counts data
df = read.csv2(opt$counts, header=TRUE, row.names=1, sep=",")

# Check and clean rownames (gene names)
message("Original rownames (first 5): ", paste(head(rownames(df)), collapse=", "))
rownames(df) <- trimws(rownames(df))  # Remove leading/trailing whitespace
rownames(df) <- gsub("['\"]", "", rownames(df))  # Remove quotes
message("Cleaned rownames (first 5): ", paste(head(rownames(df)), collapse=", "))

# convert to numeric while preserving row names
df2 <- df
df2[] <- lapply(df, function(x) as.numeric(as.character(x)))  # Using df2[] preserves row names

# round the count data to integers
countData = round(df2[, sample_names])

# create DESEq2 dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~condition)

# run DESeq2
dds = DESeq(dds)

# get normalized counts
norm_counts <- counts(dds, normalized=TRUE)

# Map Ensembl IDs to gene symbols for the normalized count matrix
ensembl_genes <- rownames(norm_counts)

# Add a check to ensure we have valid Ensembl IDs
message("First few row names: ", paste(head(ensembl_genes), collapse=", "))

# Filter for only valid Ensembl IDs (they typically start with ENSG for human genes)
valid_ensembl <- grep("^ENSG", ensembl_genes, value = TRUE)
print(length(valid_ensembl))
if(length(valid_ensembl) == 0) {
    warning("No valid Ensembl IDs found. Check your input data format.")
    # If no valid IDs, just use original names
    unique_gene_symbols <- ensembl_genes
} else {
    print(paste0("org.Hs.eg.db columns ",columns(org.Hs.eg.db)))
    # Map only valid Ensembl IDs
    gene_symbols <- mapIds(
        org.Hs.eg.db,
        keys = valid_ensembl,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
    )
    print(paste0("gene symbols: ",head(gene_symbols), " ensembl IDs: ", head(names(gene_symbols))))
    gene_symbols.df = as.data.frame(as.matrix(gene_symbols))
    print(head(gene_symbols.df))
    # Create a named vector for all genes, including non-Ensembl IDs
    ensembl_genes <- unlist(ensembl_genes)
    ens_vec <- unlist(rownames(gene_symbols.df))
    colnames(gene_symbols.df) <- "gene"
    print(head(gene_symbols.df))
    sym_vec <- unlist(gene_symbols.df$gene)
    all_symbols <- ensembl_genes
    
    all_symbols[ens_vec] <- sym_vec
    unique_gene_symbols <- make.unique(as.character(all_symbols))
    print(paste0("unique gene symbols: ",head(unique_gene_symbols)))
}

# Add gene symbols as the first column and write the normalized counts
norm.counts <- as.data.frame(norm_counts)
head(rownames(norm_counts))
norm.counts <- merge(gene_symbols.df, norm.counts, by="row.names", all=TRUE)
#print(head(norm_counts))
#norm_counts <- data.frame(gene = unique_gene_symbols, as.data.frame(norm_counts),row.names = NULL)
write.csv(norm_counts, file = file.path(opt$output_dir, "normalized_counts.csv"), row.names = FALSE, quote = FALSE)
message("# Normalized count matrix written to:", file.path(opt$output_dir, "normalized_counts.csv"))

# Save the individual contrast results
for (i in 1:nrow(contrasts)) {
    contrast <- contrasts[i,]
    print(contrast)
    result_name <- paste(contrast$treatment, "vs", contrast$control, sep = "_")
    output_file <- paste0(opt$output_dir, "/", result_name, ".csv")
    
    # Extract DESeq2 results for this contrast
    res <- results(dds, contrast = c("condition", contrast$treatment, contrast$control))
    result_df <- as.data.frame(res)
    print(paste0("contrast result dim: ", dim(result_df)))
    # Add gene symbols to the contrast results
    #result_df$gene <- unique_gene_symbols[1:nrow(result_df)]
    result_df <- merge(gene_symbols.df, result_df, by="row.names", all=TRUE)  # merge by row names (by=0 or by="row.names")
    print(head(result_df))
    # Subset normalized counts for treatment and control groups
    control_samples <- design[design$condition == contrast$control, "sample"]
    treatment_samples <- design[design$condition == contrast$treatment, "sample"]
    
    print("Debug: Structure of norm_counts")
    print(str(norm_counts))
    print("Debug: beginning of norm_counts")
    print(head(norm_counts))
    
    # Get the normalized counts matrix without the gene column
    #norm_counts_matrix <- as.matrix(norm_counts[, !colnames(norm_counts) %in% "gene"])
    #rownames(norm_counts_matrix) <- norm_counts[["gene"]]  # Use [[ ]] instead of $
    
    # Subset the normalized counts for control and treatment samples
    #normalized_control_counts <- norm_counts[rownames(result_df), control_samples, drop = FALSE]
    #normalized_treatment_counts <- norm_counts[rownames(result_df), treatment_samples, drop = FALSE]
    normalized_control_counts <- norm_counts[result_df$Row.names, control_samples, drop = FALSE]
    normalized_treatment_counts <- norm_counts[result_df$Row.names, treatment_samples, drop = FALSE]

    print(paste0("norm ctrl dim: ", dim(normalized_control_counts)))
    print(paste0("norm treatment dim: ", dim(normalized_treatment_counts)))
    
    # Add normalized counts for treatment and control groups to the results
    result_df <- cbind(result_df, normalized_control_counts, normalized_treatment_counts)
    print(head(result_df))
    
    # Reorder the columns to have the gene symbol as the first column
    result_df <- result_df[, c("gene", setdiff(names(result_df), "gene"))]
    
    # Write the results to a file
    write.csv(result_df, file = output_file, row.names = FALSE, quote = FALSE)
    message(paste("# Output for contrast", result_name, "written to:", output_file))
}