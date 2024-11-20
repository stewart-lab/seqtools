suppressPackageStartupMessages(library("gprofiler2"))
suppressPackageStartupMessages(library("optparse"))


# Define the command-line options
option_list <- list(
    make_option(c("--genes"), type = "character", default = NULL,
                help = "Path to the gene list file", metavar = "character"),
    make_option(c("--organism"), type = "character", default = NULL,
                help = "Organism to use for enrichment analysis", metavar = "character"),
    make_option(c("--mthreshold"), type = "numeric", default = 500,
                help = "Max members per term", metavar = "numeric"),
    make_option(c("--output_dir"), type = "character", default = NULL,
                help = "Path to the output directory", metavar = "character")
)


# Parse the command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# read the gene list (header col = 'genes')
data <- read.csv(opt$genes, header = TRUE, stringsAsFactors = FALSE)

# perform the enrichment analysis
gostres <- gost(
    data$gene, 
    organism = opt$organism, 
    significant = TRUE,
    evcodes = TRUE
)

# Flatten list columns in the results
flattened_results <- gostres$result
list_columns <- sapply(flattened_results, is.list)

# Convert list columns to strings (comma-separated)
flattened_results[list_columns] <- lapply(flattened_results[list_columns], function(column) {
  sapply(column, function(entry) {
    if (is.null(entry)) {
      return(NA)  # Handle NULL entries
    }
    paste(entry, collapse = ", ")  # Concatenate list elements into a single string
  })
})

# Sort the results by recall (descending order)
flattened_results <- flattened_results[order(flattened_results$recall, decreasing = TRUE), ]

# remove any terms with term_size >= mthreshold
flattened_results <- flattened_results[flattened_results$term_size < opt$mthreshold, ]


# get the gene list for each significant term - call gConvert
gconvertres <- gconvert(
    query = flattened_results$term_id,
    organism = opt$organism,
    target = "HGNC",
    mthreshold = Inf,
    filter_na = TRUE
)

# Create a mapping from term_id to genes
gene_map <- split(gconvertres$target, gconvertres$input)

# Add a new column 'term_members' to flattened_results
flattened_results$term_members <- sapply(flattened_results$term_id, function(id) {
  genes <- gene_map[[id]]
  if (is.null(genes)) {
    return(NA)  # If no genes are found for the term_id
  } else {
    paste(genes, collapse = ";")  # Concatenate gene names into a single string
  }
})


# write the results to a file. file name is same as 'genes' input file
file_name <- tools::file_path_sans_ext(basename(opt$genes))
write.csv(flattened_results, file.path(opt$output_dir, paste0(file_name, ".csv")), row.names = FALSE)