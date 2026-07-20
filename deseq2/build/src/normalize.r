suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("optparse"))

# Define the command-line options
option_list <- list(
    make_option(c("--counts"),
        type = "character", default = NULL,
        help = "Path to the counts CSV file", metavar = "character"
    ),
    make_option(c("--output_dir"),
        type = "character", default = NULL,
        help = "Path to the output directory", metavar = "character"
    )
)

# Parse the command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that required arguments are provided
if (is.null(opt$counts)) {
    stop("The counts file must be provided using the --counts flag")
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
# remove trailing slash from output_dir if it exists
opt$output_dir <- gsub("/$", "", opt$output_dir)

# make the output directory
dir.create(opt$output_dir, showWarnings = FALSE)

# read the counts data
df <- read.csv2(opt$counts, header = TRUE, row.names = 1, sep = ",")

# Check and clean rownames (gene names)
message("Original rownames (first 5): ", paste(head(rownames(df)), collapse = ", "))
rownames(df) <- trimws(rownames(df)) # Remove leading/trailing whitespace
rownames(df) <- gsub("['\"]", "", rownames(df)) # Remove quotes
message("Cleaned rownames (first 5): ", paste(head(rownames(df)), collapse = ", "))

# convert to numeric while preserving row names
df2 <- df
df2[] <- lapply(df, function(x) as.numeric(as.character(x))) # Using df2[] preserves row names

# round the count data to integers
countData <- round(df2)

# get sample names
sample_names <- colnames(countData)

# make deseq2 object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = data.frame(sample = sample_names), design = ~1)

# run DESeq2
dds <- DESeq(dds)

# get normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# save normalized counts
write.csv(norm_counts, file.path(opt$output_dir, "normalized_counts.csv"), row.names = TRUE)

# write packages
writeLines(capture.output(sessionInfo()), file.path(opt$output_dir, "sessionInfo.txt"))
