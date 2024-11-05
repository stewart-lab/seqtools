suppressPackageStartupMessages(library("tximport"))
suppressPackageStartupMessages(library("optparse"))

# Define the command-line options
option_list <- list(
    make_option(c("--rsem_dir"), type = "character", default = NULL,
                help = "Path to the RSEM results directory", metavar = "character"),
    make_option(c("--output_dir"), type = "character", default = NULL,
                help = "Directory to write the counts CSV file to", metavar = "character")
)

# Parse the command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that required arguments are provided
if (is.null(opt$rsem_dir)) {
  stop("The RSEM results directory must be provided using the --rsem_dir flag")
}

if (is.null(opt$output_dir)) {
  stop("The output directory must be provided using the --output_dir flag")
}

# read in the RSEM results files
files <- list.files(opt$rsem_dir, pattern="*.genes.results", full.names=TRUE)
names(files) <- gsub(".genes.results", "", basename(files))
txi.rsem <- tximport(files, type="rsem", txIn=FALSE, txOut=FALSE)

# write the counts to a CSV file
write.csv(txi.rsem$counts, file=file.path(opt$output_dir, "counts.csv"))
