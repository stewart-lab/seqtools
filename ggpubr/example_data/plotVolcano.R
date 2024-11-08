suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))

# Define command-line options
option_list <- list(
    make_option(c("--data"), type = "character", default = NULL,
                help = "Path to the data CSV file", metavar = "character"),
    make_option(c("--output"), type = "character", default = "volcano_plot.png",
                help = "Output file name for the plot", metavar = "character"),
    make_option(c("--pval_cutoff"), type = "numeric", default = NULL,
                help = "p-value cutoff for significance", metavar = "numeric"),
    make_option(c("--l2fc_cutoff"), type = "numeric", default = 1.0,
                help = "Log2 fold change cutoff for significance", metavar = "numeric")
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

# print the first few lines
print("Read in this data matrix:")
print(head(data))

# Check for required columns
if (!all(c("pvalue", "log2FoldChange") %in% colnames(data))) {
    stop("The input data must contain columns named 'pvalue' and 'log2FoldChange'")
}

# Create volcano plot
volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = (abs(log2FoldChange) >= opt$l2fc_cutoff) & (pvalue <= opt$pval_cutoff)),
               size = 2, alpha = 0.6) +
    scale_color_manual(values = c("gray", "#c41f1f")) +
    geom_hline(yintercept = -log10(opt$pval_cutoff), linetype = "dashed", color = "#185295") +
    geom_vline(xintercept = c(-opt$l2fc_cutoff, opt$l2fc_cutoff), linetype = "dashed", color = "#185295") +
    labs(title = "Volcano Plot",
         x = "Log2 Fold Change",
         y = "-Log10 P-value",
         color = "Significance") +
    theme_minimal() + 
    theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA)
    )

# Save plot
ggsave(filename = opt$output, plot = volcano_plot, width = 10, height = 8, dpi = 300)
