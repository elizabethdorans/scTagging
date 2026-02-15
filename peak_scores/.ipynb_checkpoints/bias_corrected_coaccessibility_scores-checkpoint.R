suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--uncorrected_file",
                    help = "path to file containing uncorrected ArchR correlations")
parser$add_argument("--blacklist_peaks_file",
                   help = "Single-column text file with list of peaks to exclude (these peaks will not be considered in peak-peak coaccessibilities)")
parser$add_argument("--background_coacc_file",
                    help = "path to file containing background coasccessibility")
parser$add_argument("--outfile",
                   help = "file to send output")

args <- parser$parse_args()

uncorrected_file = args$uncorrected_file
background_coacc_file = args$background_coacc_file
blacklist_peaks_file = args$blacklist_peaks_file
outfile = args$outfile

# Read in peak-peak correlations
coacc = read.table(uncorrected_file, sep = "\t", header = TRUE)
head(coacc)

if (!is.null(blacklist_peaks_file)) {
    blacklist_peaks = read.table(blacklist_peaks_file, sep = "\t")$V1
    blacklist_peaks = paste0("chr", blacklist_peaks)
    print(sprintf("Restricting to non-blacklisted peaks only! From %s peaks...", length(unique(coacc$peak1))))
    coacc = coacc[!(coacc$peak1 %in% blacklist_peaks) & !(coacc$peak2 %in% blacklist_peaks),]
    print(sprintf("...to %s peaks!", length(unique(coacc$peak1))))
}

# Get squared correlations
coacc["corr_sq"] = coacc["correlation"] ** 2

# Correct each r2 by subtracting (1 - r2) * bias
bias = read.table(background_coacc_file, header = TRUE)$noise[1]
coacc["corrected_corr_sq"] = coacc["corr_sq"] - bias
coacc["cis_peaks"] = 1

# Compute scores by summing corrected r2
coacc_scores <- aggregate(coacc[c("corrected_corr_sq", "cis_peaks")], by=list(peak=coacc$peak1), FUN=sum)

# Add 1 to incude the focal peak
coacc_scores$coaccessibility_score <- coacc_scores$corrected_corr_sq
coacc_scores$corrected_corr_sq <- NULL

coacc_scores$coaccessibility_score <- coacc_scores$coaccessibility_score + 1
coacc_scores$cis_peaks <- coacc_scores$cis_peaks + 1

# Save coaccessibility scores
write.table(coacc_scores, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
print("Done!")