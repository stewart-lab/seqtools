# Load required libraries
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("optparse"))

# Define command-line options
option_list <- list(
  make_option(c("--data"), type = "character", default = NULL,
              help = "Path to the data CSV file", metavar = "character"),
  make_option(c("--output"), type = "character", default = "scatter_plot.png",
              help = "Output file path for the scatter plot", metavar = "character"),
  make_option(c("--title"), type = "character", default = "Scatter Plot with Trendline",
              help = "Title for the scatter plot", metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that required arguments are provided
if (is.null(opt$data)) {
  stop("The data file must be provided using the --data flag")
}

# Check if the file exists and is a CSV file
check_csv_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste("The file does not exist:", file_path))
  }
  if (tools::file_ext(file_path) != "csv") {
    stop(paste("The file must have a .csv extension:", file_path))
  }
}

check_csv_file(opt$data)

# Load the data
data <- read.csv(opt$data, header = TRUE)

# Ensure there are at least three columns
if (ncol(data) < 3) {
  stop("The input CSV file must have at least three columns.")
}

# Extract the second and third columns for x and y
x_col <- colnames(data)[2]
y_col <- colnames(data)[3]
data$x <- data[[x_col]]
data$y <- data[[y_col]]

# Calculate Pearson correlation coefficient
pearson_corr <- cor(data$x, data$y, use = "complete.obs", method = "pearson")
pearson_label <- paste0("Pearson r = ", round(pearson_corr, 2))

# Determine axis limits
axis_limit <- max(2, max(abs(data$x), abs(data$y), na.rm = TRUE))

# Create scatter plot
scatter_plot <- ggplot(data, aes(x = x, y = y)) +
  geom_point(size = 2, alpha = 0.6, color = "#185295") +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text", x = -axis_limit + 0.5, y = axis_limit - 0.5, label = pearson_label, 
           hjust = 0, vjust = 1, size = 5, color = "black") +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed", size = 0.3) +  # Add dashed y-axis
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed", size = 0.3) +  # Add dashed x-axis
  scale_x_continuous(limits = c(-axis_limit, axis_limit)) +
  scale_y_continuous(limits = c(-axis_limit, axis_limit)) +
  labs(title = opt$title,
       x = x_col,
       y = y_col) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save the plot to a file
ggsave(filename = opt$output, plot = scatter_plot, width = 8, height = 8, dpi = 300)

