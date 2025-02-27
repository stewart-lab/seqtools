suppressPackageStartupMessages(library("gprofiler2"))
suppressPackageStartupMessages(library("optparse"))


# Define the command-line options
option_list <- list(
    make_option(c("--gmt"), type = "character", default = NULL,
                help = "Path to gmt file to upload", metavar = "character")
)


# Parse the command-line options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!file.exists(opt$gmt)) {
    stop("The provided GMT file does not exist.")
}

gmt_id <- upload_GMT_file(gmtfile = opt$gmt)[[1]]

# print(gmt_id)

cat(paste(gmt_id, "\n"))