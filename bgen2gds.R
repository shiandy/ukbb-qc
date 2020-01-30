# ======================================================================
# Convert bgen file to gds to read into R.
# ======================================================================

# -------------------------------------------
# Parse command line arguments
# -------------------------------------------
library("optparse")

parser <- OptionParser(description = "Convert bgen file to gds.")
parser <- add_option(parser, c("-b", "--bgen"), type = "character",
                     help = "Input bgen file.")
parser <- add_option(parser, c("-o", "--output"), type = "character",
                     help = "Output gds file.")
parser <- add_option(parser, c("-t", "--threads"), type = "integer",
                     default = 1,
                     help = "Number of threads to run [default: %default]")
opt <- parse_args(parser)

if (!file.exists(opt$bgen)) {
  stop(sprintf("Input bgen file %s does not exist.", opt$bgen))
}

if (opt$threads < 1) {
  stop("Number of threads must be >= 1.")
}

# -------------------------------------------
# Do the conversion
# -------------------------------------------
library("gds2bgen")

seqBGEN_Info(opt$bgen)

seqBGEN2GDS(opt$bgen, opt$output, storage.option = "LZMA_RA",
            float.type = "packed16", geno = FALSE, dosage = TRUE,
            prob = FALSE, parallel = opt$threads)
