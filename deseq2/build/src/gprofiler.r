suppressPackageStartupMessages(library("gprofiler2"))
suppressPackageStartupMessages(library("optparse"))
library(dplyr)
library("ggplot2")
library(reshape2)
library(devEMF)
library(ggrepel)
library(tidyr)

# Define the command-line options
option_list <- list(
    make_option(c("--genes"), type = "character", default = NULL,
                help = "Path to the gene list file", metavar = "character"),
    make_option(c("--organism"), type = "character", default = NULL,
                help = "Organism to use for enrichment analysis, hsapiens, 
                or gmt_ids: GO_BP: gp__httq_So1p_ckU, KEGG legacy: gp__KbhP_AZWh_sxU, 
                KEGG medicus: gp__hHRY_iYNB_uoQ"
                , metavar = "character"), 
    make_option(c("--mthreshold"), type = "numeric", default = 500,
                help = "Max members per term", metavar = "numeric"),
    make_option(c("--padj"), type = "numeric", default = 0.05,
                help = "Minimum adj p-value for genes to use", metavar = "numeric"),
    make_option(c("--lfc"), type = "numeric", default = 1,
                help = "Minimum log2 fold change value for genes to use", metavar = "numeric"),
    make_option(c("--output_dir"), type = "character", default = NULL,
                help = "Path to the output directory", metavar = "character"),
    make_option(c("--output_name"), type = "character", default = "_GO_enrich",
                help = "name appended to end of output file", metavar = "character"),
    make_option(c("--title"), type = "character", default = "GO enrichment",
                help = "Title of GO term plot", metavar = "character")
                )



# Parse the command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Extract the base filename without the extension
output_file_name <- tools::file_path_sans_ext(basename(opt$genes))

# get the output file name
output_file_name <- paste0(output_file_name, opt$output_name)


# read the gene list (header col = 'genes')
data <- read.csv(opt$genes, header = TRUE, stringsAsFactors = FALSE)
print(colnames(data))
# filter 
if(sign(opt$lfc)==1){
     data <- subset(data, log2FoldChange >= opt$lfc & padj < opt$padj)
}else if (sign(opt$lfc)==-1) {
     data <- subset(data, log2FoldChange <= opt$lfc & padj < opt$padj)
}

# check for empty dataframe
if(dim(data)[1] == 0){
  print("No DE genes to analyze")
} else {
############### enrichment analysis via gProfiler ################
# perform the enrichment analysis
gostres <- gost(
    data$gene, 
    organism = opt$organism,
    measure_underrepresentation = FALSE, 
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

# get adjusted p
flattened_results$p.adj <- p.adjust(flattened_results$p_value, method = "fdr", n = length(flattened_results$p_value))

# remove any terms with term_size >= mthreshold
flattened_results <- flattened_results[flattened_results$term_size < opt$mthreshold, ]

print(flattened_results)
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
write.csv(flattened_results, file.path(opt$output_dir, 
          paste0(output_file_name, ".csv")), row.names = FALSE)
}

# make plot

bar_data <- flattened_results[,c("term_name","p.adj")]
print(bar_data)
# add neg.log.adj.p
bar_data$neg.log.adj.pvalue <- -log(bar_data$p.adj)
# replace NAs from log(0)
bar_data$neg.log.adj.pvalue %>% replace_na(300)
# subset significant ones
bar_data <- subset(bar_data, neg.log.adj.pvalue>=1.3)

# Remove "GOBP" prefix and split term names by "_" and take the last 3 parts
bar_data$term_name <- gsub("^GOBP_", "", as.character(bar_data$term_name))
bar_data$term_name <- sapply(strsplit(as.character(bar_data$term_name), "_"), function(x) {
    last_parts <- tail(x, 3)
    paste(last_parts, collapse="_")
})

# Keep only the term with lowest p-value (highest neg.log.adj.pvalue) for duplicates
bar_data <- bar_data %>%
  group_by(term_name) %>%
  slice(which.max(neg.log.adj.pvalue)) %>%
  ungroup()

# plot
p1<-ggplot(bar_data, aes(x=reorder(term_name, neg.log.adj.pvalue), y=neg.log.adj.pvalue)) + 
  geom_col(aes(fill=neg.log.adj.pvalue), position="identity") + 
  theme_minimal()+ 
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "red"))(100))+
  labs(title= opt$title, x ="Term", y = "Neg Log adjusted P-value")+
  coord_flip()

# Print data for debugging
print("Checking for any duplicate terms or unusual values:")
print(bar_data[, c("term_name", "neg.log.adj.pvalue")])

# make pdf
nd=file.path(opt$output_dir, paste0(output_file_name,"_barplot.pdf"))
pdf(file = nd, height=11, width=8.5)
print(p1)
dev.off()