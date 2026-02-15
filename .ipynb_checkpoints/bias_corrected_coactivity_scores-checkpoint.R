suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--uncorrected_file",
                    help = "path to file containing uncorrected ArchR correlations")
parser$add_argument("--blacklist_peaks_file",
                   help = "Single-column text file with list of peaks to exclude (these peaks will not be considered in peak-gene correlations)")
parser$add_argument("--background_coacc_file",
                    help = "path to file containing uncorrected ArchR correlations")
parser$add_argument("--gene_universe_file",
                   help = "file with first column 'gene'",
                   default = "gene_TSS.txt")
parser$add_argument("--focal_feature",
                   help = "level at which coactivity score is defined ('peak' or 'gene')",
                   default = "peak")
parser$add_argument("--outfile",
                   help = "file to send output")

args <- parser$parse_args()

uncorrected_file = args$uncorrected_file
background_coacc_file = args$background_coacc_file
blacklist_peaks_file = args$blacklist_peaks_file
gene_universe_file = args$gene_universe_file
focal_feature = args$focal_feature
outfile = args$outfile

# Get filename
if (!dir.exists(dirname(outfile))) {
    dir.create(dirname(outfile))
}

coacc = read.table(uncorrected_file, sep = "\t", header = TRUE)
head(coacc)

# Filter to non-promoter peaks
if (!is.null(blacklist_peaks_file)) {
    blacklist_peaks = read.table(blacklist_peaks_file, sep = "\t")$V1
    blacklist_peaks = paste0("chr", blacklist_peaks)
    print(sprintf("Restricting to non-blacklisted peaks only! From %s peaks...", length(unique(coacc$peak))))
    coacc = coacc[!(coacc$peak %in% blacklist_peaks),]
    print(sprintf("...to %s peaks!", length(unique(coacc$peak))))
}

# Filter to genes in gene universe
print(sprintf("Restricting to gene universe! From %s genes...", length(unique(coacc$gene))))
gene_universe = read.table(gene_universe_file, sep = "\t", header = TRUE)
coacc = coacc[coacc$gene %in% gene_universe$gene,]
print(sprintf("...to %s genes!", length(unique(coacc$gene))))

# Square correlations
coacc["corr_sq"] = coacc[corr_colname] ** 2

# Correct each r2 by subtracting bias
bias = read.table(background_coacc_file, header = TRUE)$noise[1]
coacc["corrected_corr_sq"] = coacc["corr_sq"] - bias
coacc["cis_genes"] = 1
coacc['focal_feature'] <- coacc[, focal_feature]

# Compute scores
coacc_scores <- aggregate(coacc[c("corrected_corr_sq", "cis_genes")], by=list(focal_feature=coacc$focal_feature), FUN=sum)
coacc_scores$coactivity_score <- coacc_scores$corrected_corr_sq
coacc_scores$corrected_corr_sq <- NULL

# Save scores
write.table(coacc_scores, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
print("Done!")