suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("circlize"))

# get the .csv from optparse
# Define the command-line options
option_list <- list(
    make_option(c("--data"), type = "character", default = NULL,
                help = "Path to the data CSV file", metavar = "character")
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

# Load the data for the heatmap. header = sample names, row.names = gene names
data <- read.csv(opt$data, header = TRUE, row.names = 1)

# Convert data frame to matrix
data_matrix <- as.matrix(data)

# Define custom color function
col_fun <- colorRamp2(c(-2.0, -0.3, 0.3, 2.0), c("blue", "white", "white", "red"))

# TODO: get SVG, PNG, PDF, etc.
png(filename = "/mnt/data/heatmap.png", width = 3.25 * 300, height = 3.25 * 300, res = 300)

# Plot heatmap without clustering
Heatmap(data_matrix,
        name = "Expression",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        column_names_side = "top")

# you can use sigclust2 to get statistical significance of clusters in a clustermap. TODO.

# Close the device to save the file
dev.off()