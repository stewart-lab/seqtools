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

# Extract the base filename without the extension
output_file_name <- tools::file_path_sans_ext(basename(opt$genes))

# Replace '_significant' with '_pathways' in the file name, if present to get the output file name
output_file_name <- sub("_significant$", "_pathways", output_file_name)


# read the gene list (header col = 'genes')
data <- read.csv(opt$genes, header = TRUE, stringsAsFactors = FALSE)

############### enrichment analysis via gProfiler ################
# perform the enrichment analysis
gostres <- gost(
    data$gene, 
    organism = opt$organism, 
    significant = TRUE
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

# Sort the results by p_value (descending order)
flattened_results <- flattened_results[order(flattened_results$p_value, decreasing = FALSE), ]

# remove any terms with term_size >= mthreshold
flattened_results <- flattened_results[flattened_results$term_size < opt$mthreshold, ]

# if there are no results, write 'no significant terms found' to the output file
if (nrow(flattened_results) == 0) {
  # write.csv(data.frame(term_id = NA, term_name = "No significant terms found", p_value = NA, term_size = NA), file.path(opt$output_dir, paste0(output_file_name, ".csv")), row.names = FALSE)
  stop("No significant terms found")
}

# if any term_name contains 'https://', swap the term_id and term_name columns.
# this isn't really a 'bug' in gProfiler per se, but it uses the second column
# of a .gmt file (the description) as the term_name. but the first column (the term ID)
# is probably more representative of the term name.
first_term_name <- flattened_results$term_name[1]
if (grepl("https://", first_term_name)) {
  temp <- flattened_results$term_name
  flattened_results$term_name <- flattened_results$term_id
  flattened_results$term_id <- temp
}

# Write the results to the output file
write.csv(flattened_results, file.path(opt$output_dir, paste0(output_file_name, ".csv")), row.names = FALSE)


# # if the organism contains an underscore, assume it is a custom organism and stop the script here
# if (grepl("_", opt$organism)) {
#   stop(".gmt file (custom organism) detected. Skipping gConvert.")
# }

# ############### get gene members for each significant term ################
# gconvertres <- gconvert(
#     query = flattened_results$term_id,
#     organism = opt$organism,
#     target = "HGNC", # convert to HGNC gene symbols
#     mthreshold = Inf,
#     filter_na = TRUE
# )

# # Create a mapping from term_id to genes
# gene_map <- split(gconvertres$target, gconvertres$input)

# # Add a new column 'term_members' to flattened_results
# flattened_results$term_members <- sapply(flattened_results$term_id, function(id) {
#   genes <- gene_map[[id]]
#   if (is.null(genes)) {
#     return(NA)  # If no genes are found for the term_id
#   } else {
#     paste(genes, collapse = ";")  # Concatenate gene names into a single string
#   }
# })

# # Write the results to the output file
# write.csv(flattened_results, file.path(opt$output_dir, paste0(output_file_name, ".csv")), row.names = FALSE)