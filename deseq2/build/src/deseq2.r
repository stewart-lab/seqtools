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
design <- read.csv(opt$design, stringsAsFactors=F)
design$condition = factor(design$condition)

# Read the contrasts file.
contrasts <- read.csv(opt$contrasts)

# Isolate the sample names.
sample_names <- design$sample

# read the counts data
df = read.csv(opt$counts, header=TRUE, row.names=1)

# round the count data to integers
countData = round(df[, sample_names])

# create DESEq2 dataset
dds = DESeqDataSetFromMatrix(countData=countData, colData=design, design=~condition)

# run DESeq2
dds = DESeq(dds)

# get normalized counts
norm_counts <- counts(dds, normalized=TRUE)

# save the individual contrast results
for (i in 1:nrow(contrasts)) {
    contrast <- contrasts[i,]

    # make the output file name for this contrast
    result_name <- paste(contrast$treatment, "vs", contrast$control, sep="_")
    output_file <- paste0(opt$output_dir, "/", result_name, ".csv")

    # extract DESeq2 results for this contrast and convert to data frame
    res <- results(dds, contrast=c("condition", contrast$treatment, contrast$control))
    result_df <- as.data.frame(res)

    # extract the sample names for the contrast's treatment and control groups
    treatment_design_rows <- design[design$condition == contrast$treatment, ]
    control_design_rows <- design[design$condition == contrast$control, ]

    treatment_samples <- treatment_design_rows$sample
    control_samples <- control_design_rows$sample

    # subset normalized counts for treatment and control groups
    normalized_treatment_counts <- norm_counts[, treatment_samples]
    normalized_control_counts <- norm_counts[, control_samples]

    # Add normalized counts for each sample in treatment and control groups to result_df
    result_df <- cbind(result_df, normalized_treatment_counts, normalized_control_counts)

    ensembl_genes <- rownames(result_df)

    gene_symbols <- mapIds(org.Hs.eg.db,
                       keys=ensembl_genes,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")

    # Replace NAs with original Ensembl IDs
    gene_symbols[is.na(gene_symbols)] <- names(gene_symbols)[is.na(gene_symbols)]

    # Make unique row names if there are duplicates
    unique_gene_symbols <- make.unique(as.character(gene_symbols))

    # Add gene symbols to the result_df
    row.names(result_df) <- unique_gene_symbols

    # add gene names as the first column
    result_df <- cbind(gene=rownames(result_df), result_df)

    # write individual result to file
    write.csv(result_df, file=output_file, row.names=FALSE, quote=FALSE)
    
    message(paste("# Output for contrast", result_name, "written to:", output_file))
}
