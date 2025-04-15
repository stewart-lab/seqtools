suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))

# Define command-line options
option_list <- list(
    make_option(c("--data"), type = "character", default = NULL,
                help = "Path to the data CSV file", metavar = "character"),
    make_option(c("--output"), type = "character", default = "volcano_plot.png",
                help = "Output file name for the plot", metavar = "character"),
    make_option(c("--padj_cutoff"), type = "numeric", default = NULL,
                help = "padj cutoff for significance", metavar = "numeric"),
    make_option(c("--l2fc_cutoff"), type = "numeric", default = 1.0,
                help = "Log2 fold change cutoff for significance", metavar = "numeric"),
    make_option(c("--title"), type = "character", default = "Volcano Plot",
                help = "Title for the plot", metavar = "character")
)

# Parse command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that required arguments are provided
if (is.null(opt$data)) {
    stop("The data file must be provided using the --data flag")
}

# Load data
data <- read.csv(opt$data, header = TRUE)
rownames(data) <- data$gene

# Check for required columns
if (!all(c("padj", "log2FoldChange") %in% colnames(data))) {
    stop("The input data must contain columns named 'padj' and 'log2FoldChange'")
}

# check that at least 1 data point is present
if (nrow(data) == 0) {
    stop("Cannot create volcano plot; the input data file is empty. (data path: ", opt$data, ")")
}

# Calculate axis limits
x_min <- min(-5, min(data$log2FoldChange, na.rm = TRUE))
x_max <- max(5, max(data$log2FoldChange, na.rm = TRUE))
y_max <- max(35, max(-log10(data$padj), na.rm = TRUE))

# add a 'significant' column to the data frame
data$significant <- (abs(data$log2FoldChange) >= opt$l2fc_cutoff) & (data$padj <= opt$padj_cutoff)

# Create a color palette. if the data is significant, use red, otherwise use gray
color_palette <- ifelse(data$significant, "#c41f1f", "gray")

# Create volcano plot
volcano_plot <- ggplot(data, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("TRUE" = "#c41f1f", "FALSE" = "gray")) +
    geom_hline(yintercept = -log10(opt$padj_cutoff), linetype = "dashed", color = "#185295") +
    geom_vline(xintercept = c(-opt$l2fc_cutoff, opt$l2fc_cutoff), linetype = "dashed", color = "#185295") +
    scale_x_continuous(limits = c(x_min, x_max)) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(title = opt$title,
         x = "Log2 Fold Change",
         y = "-Log10 Adj. P-value",
         color = "Significance") +
    theme_minimal() +
    theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA)
    )

# save as .png
png_filename <- ifelse(grepl(".png$", opt$output), opt$output, paste0(opt$output, ".png"))
ggsave(filename = png_filename, plot = volcano_plot, width = 10, height = 8, dpi = 300)
