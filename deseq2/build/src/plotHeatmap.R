# see documentation for ComplexHeatmap to see annotation options, etc.
# https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))

# Define the command-line options
option_list <- list(
    make_option(c("--data"), type = "character", default = NULL,
                help = "Path to the data CSV file", metavar = "character"),
    make_option(c("--gene_set"), type = "character", default = NULL,
                help = "Path to a CSV file with a list of genes to filter the data", metavar = "character"),
    
    make_option(c("--min"), type = "numeric", default = -2.0,
                help = "Minimum value for the heatmap color scale", metavar = "numeric"),
    make_option(c("--mid"), type = "numeric", default = 0.0,
                help = "Middle value for the heatmap color scale", metavar = "numeric"),
    make_option(c("--max"), type = "numeric", default = 2.0,
                help = "Maximum value for the heatmap color scale", metavar = "numeric"),
    make_option(c("--min_color"), type = "character", default = "#185295",
                help = "Hex code for the minimum color in the heatmap", metavar = "character"),
    make_option(c("--mid_color"), type = "character", default = "white",
                help = "Hex code for the middle color in the heatmap", metavar = "character"),
    make_option(c("--max_color"), type = "character", default = "#c41f1f",
                help = "Hex code for the maximum color in the heatmap", metavar = "character"),
    make_option(c("--cluster_rows"), type = "logical", default = FALSE,
                help = "Whether to cluster the rows of the heatmap", metavar = "logical"),
    make_option(c("--cluster_columns"), type = "logical", default = FALSE,
                help = "Whether to cluster the columns of the heatmap", metavar = "logical"),
    make_option(c("--aspect_ratio"), type = "numeric", default = NULL,
                help = "Aspect ratio for the heatmap plot cells (w/h)", metavar = "numeric"),
    make_option(c("--log2"), type = "logical", default = FALSE,
                help = "Whether to log2 transform the data before plotting", metavar = "logical"),
    make_option(c("--output"), type = "character", default = "/mnt/data/heatmap.png",
                help = "Output file path for the heatmap plot", metavar = "character"),
    make_option(c("--hide_row_names"), type = "logical", default = TRUE,
                help = "Whether to hide row names in the heatmap", metavar = "logical"),
    make_option(c("--title"), type = "character", default = "",
                help = "Title for the heatmap plot", metavar = "character"),
    make_option(c("--legend_title"), type = "character", default = "Normalized Log2 Counts",
                help = "Title for the heatmap legend", metavar = "character"),
    make_option(c("--subtract_row_means"), type = "logical", default = FALSE,
                help = "Whether to subtract row means from the data before plotting", metavar = "logical")

)

# Parse the command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that required arguments are provided
if (is.null(opt$data)) {
  stop("The data file must be provided using the --data flag")
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

# Check that the data file exists and is a CSV file
check_csv_file(opt$data, "data")

# Set the random seed for reproducibility
set.seed(42)

# Read the data from the CSV file
data <- read.csv(opt$data, header = TRUE)

# Set row names to the "gene" column and remove it from data frame
if ("gene" %in% colnames(data)) {
  rownames(data) <- data$gene
  data <- data[, !colnames(data) %in% "gene"]
} else {
  stop("ERROR: The 'gene' column was not found in the header.")
}

# Remove specified columns
columns_to_remove <- c("baseMean", "lfcSE", "log2FoldChange", "stat", "pvalue", "padj")
data <- data[, !colnames(data) %in% columns_to_remove]


# Filter rows based on the gene_set if provided
if (!is.null(opt$gene_set)) {
  check_csv_file(opt$gene_set, "gene_set")
  gene_set <- read.csv(opt$gene_set, header = FALSE, stringsAsFactors = FALSE)[, 1]
  data <- data[rownames(data) %in% gene_set, ]
}


# Log2 transform the data if specified
if (opt$log2) {
  data <- log2(data + 1)
}

# Subtract row means if specified
if (opt$subtract_row_means) {
  data <- data - rowMeans(data)
}

# Convert data frame to matrix
data_matrix <- as.matrix(data)

# Define custom color function
col_fun <- colorRamp2(c(opt$min, opt$mid, opt$max), c(opt$min_color, opt$mid_color, opt$max_color))

h = 10 # inches
w = 10 # inches

column_count = ncol(data_matrix)
row_count = nrow(data_matrix)
data_aspect_ratio = column_count / row_count


# if aspect_ratio is provided, use it to calculate the width and height
if (!is.null(opt$aspect_ratio)) {
  desired_aspect_ratio = opt$aspect_ratio
} else {
  desired_aspect_ratio = data_aspect_ratio
}

scaler = desired_aspect_ratio / data_aspect_ratio
data_w = w * 0.7
data_h = h * 0.7

# if scaler < 1 (i.e., there are more columns than rows), scale the height
if (scaler < 1) {
    data_h = data_h * scaler
} else {
 # else scale the width
    data_w = data_w / scaler
}

# Define the heatmap object
my_heatmap <- Heatmap(data_matrix,
                      name = opt$legend_title,
                      col = col_fun,
                      cluster_rows = opt$cluster_rows,
                      cluster_columns = opt$cluster_columns,
                      clustering_method_rows = "ward.D2",
                      clustering_method_columns = "ward.D2",
                      clustering_distance_rows = "euclidean",
                      clustering_distance_columns = "euclidean",
                      row_names_side = "right",
                      column_names_side = "bottom",
                      width = unit(data_w, "inches"),
                      height = unit(data_h, "inches"),
                      row_names_gp = gpar(fontsize = 6),
                      column_names_gp = gpar(fontsize = 8),
                      show_row_names = !opt$hide_row_names)

# save as .png
png_filename <- ifelse(grepl(".png$", opt$output), opt$output, paste0(opt$output, ".png"))
png(filename = png_filename, width = w * 300, height = h * 300, res = 300)
draw(my_heatmap, column_title = opt$title, column_title_gp = gpar(fontsize = 12))
dev.off()

# save as .svg
svg_filename <- ifelse(grepl(".svg$", opt$output), opt$output, paste0(opt$output, ".svg"))
# svg(filename = svg_filename, width = w, height = h)
# draw(my_heatmap, column_title = opt$title, column_title_gp = gpar(fontsize = 12))
# dev.off()
ggsave(svg_filename, plot = grid.grabExpr(draw(my_heatmap)), width = w, height = h)






# complexHeatmapVersion <- packageVersion("ComplexHeatmap")
# current_datetime <- Sys.time()

# output_text <- paste("Heatmaps were generated with the ComplexHeatmap R package version", complexHeatmapVersion, "with the clustering 
#     method set to 'ward.D2' and the distance metric set to 'euclidean'. Heatmap bounds were set to [", opt$min, ",", opt$max, "].")

# write output_text to a file
# write(output_text, file = "/mnt/data/README_Heatmap.txt")