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
    
    make_option(c("--min"), type = "numeric", default = -1.0,
                help = "Minimum value for the heatmap color scale", metavar = "numeric"),
    make_option(c("--mid"), type = "numeric", default = 0.0,
                help = "Middle value for the heatmap color scale", metavar = "numeric"),
    make_option(c("--max"), type = "numeric", default = 1.0,
                help = "Maximum value for the heatmap color scale", metavar = "numeric"),
    make_option(c("--min_color"), type = "character", default = "#185295",
                help = "Hex code for the minimum color in the heatmap", metavar = "character"),
    make_option(c("--mid_color"), type = "character", default = "white",
                help = "Hex code for the middle color in the heatmap", metavar = "character"),
    make_option(c("--max_color"), type = "character", default = "#c41f1f",
                help = "Hex code for the maximum color in the heatmap", metavar = "character"),
    make_option(c("--cluster"), type = "logical", default = TRUE,
                help = "Whether to cluster the rows and columns of the heatmap", metavar = "logical"),
    make_option(c("--aspect_ratio"), type = "numeric", default = NULL,
                help = "Aspect ratio for the heatmap plot cells (w/h)", metavar = "numeric")

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
data <- read.csv(opt$data, header = TRUE, row.names = 1)

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

# TODO: get SVG, PNG, PDF, etc.
png(filename = "/mnt/data/heatmap.png", width = w * 300, height = h * 300, res = 300)

# Plot heatmap without clustering
Heatmap(data_matrix,
        name = "L2FC",
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 1), # white border around each cell
        cluster_rows = opt$cluster,
        cluster_columns = opt$cluster,
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        row_names_side = "right",
        column_names_side = "bottom",
        width = unit(data_w, "inches"),
        height = unit(data_h, "inches"),
        # font size for row names
        row_names_gp = gpar(fontsize = 6),
        # font size for column names
        column_names_gp = gpar(fontsize = 8),
        # show_row_names = FALSE
    )

# you can use sigclust2 to get statistical significance of clusters in a clustermap. TODO.

# Close the device to save the file
dev.off()

complexHeatmapVersion <- packageVersion("ComplexHeatmap")
# current_datetime <- Sys.time()

output_text <- paste("Heatmaps were generated with the ComplexHeatmap R package version", complexHeatmapVersion, "with the clustering 
    method set to 'ward.D2' and the distance metric set to 'euclidean'. Heatmap bounds were set to [", opt$min, ",", opt$max, "].")

# write output_text to a file
write(output_text, file = "/mnt/data/README_Heatmap.txt")

# # create output folder from the current datetime
# output_folder <- format(current_datetime, "%Y-%m-%d_%H-%M-%S")
# dir.create(file.path("/mnt/data", output_folder))
