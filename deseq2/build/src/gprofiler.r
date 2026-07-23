suppressPackageStartupMessages(library("gprofiler2"))
suppressPackageStartupMessages(library("optparse"))
library(dplyr)
library("ggplot2")
library(reshape2)
library(devEMF)
library(ggrepel)
library(tidyr)
library(rrvgo)
library(org.Dr.eg.db)
library("org.Hs.eg.db")
library(tibble)

# Define the command-line options
option_list <- list(
  make_option(c("--genes"),
    type = "character", default = NULL,
    help = "Path to the gene list file", metavar = "character"
  ),
  make_option(c("--organism"),
    type = "character", default = NULL,
    help = "Organism to use for enrichment analysis, hsapiens,
                or gmt_ids: GO_BP: gp__httq_So1p_ckU, KEGG legacy: gp__KbhP_AZWh_sxU,
                KEGG medicus: gp__hHRY_iYNB_uoQ, zebrafish: gp__HunH_2Tcy_ESg, gp__ClfJ_YDuX_eRo, drerio",
    metavar = "character"
  ),
  make_option(c("--genes_only"),
    type = "character", default = FALSE,
    help = "Input is only genes: T/F", metavar = "character"
  ),
  make_option(c("--mthreshold"),
    type = "numeric", default = 500,
    help = "Max members per term", metavar = "numeric"
  ),
  make_option(c("--padj"),
    type = "numeric", default = NULL,
    help = "Minimum adj p-value for genes to use", metavar = "numeric"
  ),
  make_option(c("--lfc"),
    type = "numeric", default = NULL,
    help = "Minimum log2 fold change value for genes to use", metavar = "numeric"
  ),
  make_option(c("--output_dir"),
    type = "character", default = "./",
    help = "Path to the output directory", metavar = "character"
  ),
  make_option(c("--output_name"),
    type = "character", default = "_GO_enrich",
    help = "name appended to end of output file", metavar = "character"
  ),
  make_option(c("--title"),
    type = "character", default = "GO enrichment",
    help = "Title of GO term plot", metavar = "character"
  ),
  make_option(c("--lower"),
    type = "character", default = FALSE,
    help = "convert gene names to lower case: T/F",
    metavar = "character"
  ),
  make_option(c("--file2"),
    type = "character", default = NULL,
    help = "Path to a second gene list file. If given, it is run through the same
                enrichment pipeline and its terms are added to the bar plots on the
                negative side (log p.adj instead of -log p.adj), opposite file1's terms.",
    metavar = "character"
  ),
  make_option(c("--color1"),
    type = "character", default = "blue",
    help = "Bar color for file1 in the bar plots (file2, when given, is always red).",
    metavar = "character"
  ),
  make_option(c("--label1"),
    type = "character", default = "R",
    help = "Legend label for file1 in the bar plots (file2, when given, is always labeled NR).",
    metavar = "character"
  )
)

# Parse the command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# run gProfiler enrichment + rrvgo term reduction for one gene-list file.
# returns list(bar_data, reduced) with un-signed (positive) scores;
# sign for the diverging plot is applied later, when combining file1/file2.
run_gprofiler_analysis <- function(gene_file, output_file_name) {
  # read the gene list (header col = 'genes')
  if (opt$genes_only == TRUE) {
    data <- read.table(gene_file, header = FALSE, stringsAsFactors = FALSE)
    colnames(data) <- "gene"
  } else {
    data <- read.table(gene_file, header = TRUE, stringsAsFactors = FALSE)
    data <- tibble::rownames_to_column(data, "gene")
  }
  print(colnames(data))
  # convert to lower
  if (opt$lower == TRUE) {
    data$gene <- tolower(data$gene)
  }
  # filter
  if (is.null(opt$lfc) && is.null(opt$padj)) {
    print("using all genes")
  } else {
    if (sign(opt$lfc) == 1) {
      data <- subset(data, log2FoldChange >= opt$lfc & padj < opt$padj)
    } else if (sign(opt$lfc) == -1) {
      data <- subset(data, log2FoldChange <= opt$lfc & padj < opt$padj)
    }
  }

  # check for empty dataframe
  if (dim(data)[1] == 0) {
    stop("No DE genes to analyze in ", gene_file)
  }

  ############### enrichment analysis via gProfiler ################
  gostres <- gost(
    data$gene,
    organism = opt$organism,
    measure_underrepresentation = FALSE,
    significant = TRUE
  )

  if (is.null(gostres) || is.null(gostres$result) || nrow(gostres$result) == 0) {
    stop("No significant GO terms found for ", gene_file)
  }

  # Flatten list columns in the results
  flattened_results <- gostres$result
  list_columns <- sapply(flattened_results, is.list)

  # Convert list columns to strings (comma-separated)
  flattened_results[list_columns] <- lapply(flattened_results[list_columns], function(column) {
    sapply(column, function(entry) {
      if (is.null(entry)) {
        return(NA) # Handle NULL entries
      }
      paste(entry, collapse = ", ") # Concatenate list elements into a single string
    })
  })

  # Sort the results by p_value (descending order)
  flattened_results <- flattened_results[order(flattened_results$p_value, decreasing = FALSE), ]

  # get adjusted p
  flattened_results$p.adj <- p.adjust(flattened_results$p_value, method = "fdr", n = length(flattened_results$p_value))

  # remove any terms with term_size >= mthreshold
  flattened_results <- flattened_results[flattened_results$term_size < opt$mthreshold, ]

  print(flattened_results)
  if (nrow(flattened_results) == 0) {
    stop("No significant terms found for ", gene_file)
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
  write.csv(flattened_results, file.path(
    opt$output_dir,
    paste0(output_file_name, ".csv")
  ), row.names = FALSE)

  ############### bar plot data ################
  bar_data <- flattened_results[, c("term_id", "term_name", "p.adj")]
  # add neg.log.adj.p
  bar_data$neg.log.adj.pvalue <- -log(bar_data$p.adj)
  # subset significant ones
  bar_data <- subset(bar_data, neg.log.adj.pvalue >= 1.3)

  # Remove "GOBP" prefix and split term names by "_" and take the last 3 parts
  bardata1 <- bar_data
  bardata1$term_name <- gsub("^GOBP_", "", as.character(bardata1$term_name))
  # bardata1$term_name <- sapply(strsplit(as.character(bardata1$term_name), "_"), function(x) {
  #   last_parts <- tail(x, 3)
  #   paste(last_parts, collapse = "_")
  # })

  # Keep only the term with lowest p-value (highest neg.log.adj.pvalue) for duplicates
  bardata1 <- bardata1 %>%
    dplyr::group_by(term_name) %>%
    dplyr::slice(which.max(neg.log.adj.pvalue)) %>%
    dplyr::ungroup()

  ############### reduce similar go terms via rrvgo ################
  go_results <- bar_data
  sim_matrix <- calculateSimMatrix(
    go_results$term_id,
    orgdb    = "org.Dr.eg.db",
    ont      = "BP", # or "MF", "CC"
    method   = "Rel" # Relevance similarity — generally recommended
  )
  scores <- setNames(-log10(go_results$p.adj), go_results$term_id)

  # rrvgo's reduceSimMatrix clusters terms via hclust, which requires at
  # least 2 terms. With too few significant terms (or terms calculateSimMatrix
  # couldn't map in orgdb), there's nothing to reduce -- pass the term(s)
  # through as-is instead of clustering.
  # calculateSimMatrix returns a scalar NA (not NULL, not a 0-row matrix)
  # when none of the terms could be mapped in orgdb, so nrow()/rownames()
  # aren't safe to call on it directly.
  has_sim_matrix <- !is.null(sim_matrix) && !is.null(dim(sim_matrix))
  n_mapped_terms <- if (has_sim_matrix) nrow(sim_matrix) else 0
  if (n_mapped_terms < 2) {
    message(
      "Only ", n_mapped_terms, " term(s) available for rrvgo reduction -- ",
      "skipping clustering for ", gene_file
    )
    # rownames() of a 0-row matrix can itself be NULL rather than
    # character(0); as.character(NULL) normalizes it so the "go" column
    # below doesn't get silently dropped from the data.frame.
    mapped_ids <- if (has_sim_matrix) as.character(rownames(sim_matrix)) else character(0)
    mapped_names <- go_results$term_name[match(mapped_ids, go_results$term_id)]
    reduced <- data.frame(
      go = mapped_ids,
      term = mapped_names,
      parentTerm = mapped_names,
      score = scores[mapped_ids],
      stringsAsFactors = FALSE
    )
  } else {
    reduced <- reduceSimMatrix(
      sim_matrix,
      scores,
      threshold = 0.7, # similarity cutoff; 0.7 = moderately aggressive collapsing
      orgdb     = "org.Dr.eg.db"
    )
  }

  print(head(reduced[, c("go", "term", "parentTerm", "score")]))

  # write out reduced terms
  write.csv(reduced, file.path(
    opt$output_dir,
    paste0(output_file_name, "_reduced.csv")
  ), row.names = FALSE)

  # Keep only the term with lowest p-value (highest score) for duplicates
  reduced <- reduced %>%
    dplyr::group_by(term) %>%
    dplyr::slice(which.max(score)) %>%
    dplyr::ungroup()

  list(bar_data = as.data.frame(bardata1), reduced = as.data.frame(reduced))
}

# apply a sign to a file's scores, keeping only the columns needed for plotting
signed_bar_data <- function(bar_data, sign) {
  data.frame(
    term_name = bar_data$term_name,
    score = sign * bar_data$neg.log.adj.pvalue
  )
}
signed_reduced <- function(reduced, sign) {
  data.frame(
    term = reduced$term,
    score = sign * reduced$score
  )
}

# order term_col's factor levels by file1's score first, then by file2's score
# for any terms not already placed by file1
order_terms_by_file <- function(df, term_col) {
  is_file1 <- df$file == "file1"
  order1 <- df[[term_col]][is_file1][order(df$score[is_file1])]

  is_file2 <- df$file == "file2"
  order2 <- df[[term_col]][is_file2][order(df$score[is_file2])]
  order2 <- order2[!(order2 %in% order1)]

  # coord_flip() puts the first factor level at the bottom of the plot and the
  # last level at the top, so file2 goes first here to land file1 on top.
  df[[term_col]] <- factor(df[[term_col]], levels = c(order2, order1))
  df
}

# run the enrichment + reduction pipeline for file1 (required) and file2 (optional)
output_file_name <- paste0(tools::file_path_sans_ext(basename(opt$genes)), opt$output_name)
results1 <- run_gprofiler_analysis(opt$genes, output_file_name)

results2 <- NULL
if (!is.null(opt$file2)) {
  output_file_name2 <- paste0(tools::file_path_sans_ext(basename(opt$file2)), opt$output_name)
  results2 <- run_gprofiler_analysis(opt$file2, output_file_name2)
}

# combine scores for plotting: file1 stays positive (-log p.adj), file2 is
# negated (log p.adj), so it lands on the opposite side of the bar plot.
combined_bar <- signed_bar_data(results1$bar_data, 1)
combined_bar$file <- rep("file1", nrow(combined_bar))
combined_reduced <- signed_reduced(results1$reduced, 1)
combined_reduced$file <- rep("file1", nrow(combined_reduced))

if (!is.null(results2)) {
  bar2 <- signed_bar_data(results2$bar_data, -1)
  bar2$file <- rep("file2", nrow(bar2))
  combined_bar <- rbind(combined_bar, bar2)

  reduced2 <- signed_reduced(results2$reduced, -1)
  reduced2$file <- rep("file2", nrow(reduced2))
  combined_reduced <- rbind(combined_reduced, reduced2)
}

# order file1/file2 rows together now that both are present, so file1 lands on top
combined_bar <- order_terms_by_file(combined_bar, "term_name")
combined_reduced <- order_terms_by_file(combined_reduced, "term")

print(combined_reduced)
# plot: full GO term bar plot
p1 <- ggplot(combined_bar, aes(x = term_name, y = score, fill = file)) +
  geom_col(position = "identity") +
  theme_minimal() +
  # scale_fill_gradientn(colours = colorRampPalette(c("blue", "red"))(100)) +
  scale_fill_manual(
    values = c("file1" = opt$color1, "file2" = "red"),
    labels = c("file1" = opt$label1, "file2" = "NR")
  ) +
  labs(
    title = opt$title, x = "Term",
    y = "log10(p.adj) NR / -log10(p.adj) R", fill = NULL
  ) +
  coord_flip()

# make pdf
nd <- file.path(opt$output_dir, paste0(output_file_name, "_combined_barplot.pdf"))
pdf(file = nd, height = 11, width = 8.5)
print(p1)
dev.off()

# plot: rrvgo-reduced GO term bar plot
p2 <- ggplot(combined_reduced, aes(x = term, y = score, fill = file)) +
  geom_col(position = "identity") +
  theme_minimal() +
  scale_fill_manual(
    values = c("file1" = opt$color1, "file2" = "red"),
    labels = c("file1" = opt$label1, "file2" = "NR")
  ) +
  labs(
    title = opt$title, x = "Term",
    y = "log10(p.adj) NR / -log10(p.adj) R", fill = NULL
  ) +
  coord_flip()

# make pdf
nd <- file.path(opt$output_dir, paste0(output_file_name, "_combined_barplot_reduced.pdf"))
pdf(file = nd, height = 11, width = 8.5)
print(p2)
dev.off()

# write out options to output, so the run can be repeated
opt_log_file <- file.path(opt$output_dir, paste0(output_file_name, "_options.txt"))
opt_lines <- sapply(names(opt), function(n) {
  if (n == "help") {
    return(NULL)
  }
  paste0("--", n, " ", opt[[n]])
})
opt_lines <- unlist(opt_lines)
writeLines(
  c(paste("# gprofiler.r options used on", Sys.time()), opt_lines),
  opt_log_file
)

# write packages
writeLines(capture.output(sessionInfo()), file.path(opt$output_dir, "sessionInfo.txt"))
