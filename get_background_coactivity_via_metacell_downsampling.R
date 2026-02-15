library(dplyr)
library(tidyr)
library(rhdf5)
library(ArchR)
suppressPackageStartupMessages(library(data.table))
library(argparse)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--norm_atac_file",
                   help = "matrix of normalized ATAC, rows = peaks, columns = cells/metacells")
parser$add_argument("--norm_rna_file",
                   help = "matrix of normalized ATAC, rows = peaks, columns = cells/metacells")
parser$add_argument("--blacklist_peaks_file",
                   help = "Single-column text file with list of peaks to exclude (these peaks will not be considered in peak-peak coaccessibilities)")
parser$add_argument("--gene_universe_file",
                   help = "file with first two columns 'gene', 'chr'",
                   default = "../gene_TSS.txt")
parser$add_argument("--num_peak_gene_pairs", default = 100000,
                   help = "number of peak-gene pairs")
parser$add_argument("--seed", default = 511,
                   help = "seed value for random number generation")
parser$add_argument("--outfile",
                   help = "file to send output")

args <- parser$parse_args()

norm_atac_file = args$norm_atac_file
norm_rna_file = args$norm_rna_file
blacklist_peaks_file = args$blacklist_peaks_file
gene_universe_file = args$gene_universe_file
num_peak_gene_pairs = as.numeric(args$num_peak_gene_pairs)
seed = as.numeric(args$seed)
outfile = args$outfile

print(outfile)
out_dir = dirname(outfile)
print(out_dir)
if (!dir.exists(out_dir)) {
    print(sprintf("Creating directory"))
    dir.create(out_dir, recursive = TRUE)
}

# Get normal ATAC and RNA matrices
norm_atac = readRDS(norm_atac_file)
norm_rna = readRDS(norm_rna_file)

norm_atac = norm_atac[Matrix::rowSums(norm_atac) != 0,]
norm_rna = norm_rna[Matrix::rowSums(norm_rna) != 0,]

if (!is.null(blacklist_peaks_file)) {
    blacklist_peaks = read.table(blacklist_peaks_file, sep = "\t")$V1
    blacklist_peaks = paste0("chr", blacklist_peaks)
    print(sprintf("Restricting to non-blacklisted peaks only! From %s peaks...", dim(norm_atac)[1]))
    norm_atac = norm_atac[!(rownames(norm_atac) %in% blacklist_peaks),]
    print(sprintf("...to %s peaks!", dim(norm_atac)[1]))
}

# Peak set with chromosomes
peaks = rownames(norm_atac)
peak_chrom = sapply(strsplit(peaks, "-"), '[[', 1)
names(peak_chrom) = peaks

# Gene set with chromosomes
gene_universe = read.table(gene_universe_file, sep = "\t", header = TRUE)
gene_universe = gene_universe[,1:2]
genes = intersect(rownames(norm_rna), gene_universe$gene)

norm_rna = norm_rna[rownames(norm_rna) %in% genes,]
gene_universe = gene_universe[gene_universe$gene %in% genes,]

genes = gene_universe$gene
gene_chrom = paste("chr", gene_universe$chr, sep = "")
names(gene_chrom) = genes

# Get metacells for 50% downsampling
set.seed(seed)
n_metacell_orig = dim(norm_atac)[2]
n_metacell_downsampled = floor(0.5 * n_metacell_orig)
sample_col_idx = sample(seq(1:n_metacell_orig), n_metacell_downsampled, replace = FALSE)

# Get 50% downsampled ATAC and RNA matrices
downsampled_atac <- norm_atac[,sample_col_idx]
print(dim(downsampled_atac))
downsampled_rna <- norm_rna[,sample_col_idx]
print(dim(downsampled_rna))

# Select random off-chromosome pairs
set.seed(seed)
peak = sample(peaks, num_peak_gene_pairs, replace = TRUE)
gene = rep(NA, num_peak_gene_pairs)
for (i in 1:length(peak)) {
    focal_peak = peak[i]
    focal_chrom = peak_chrom[[focal_peak]]
    gene[i] = sample(names(gene_chrom[gene_chrom != focal_chrom]), 1)
}

# Compute squared correlations with original and downsampled matrices
coaccessibility_sq_norm = rep(NA, num_peak_gene_pairs)
coaccessibility_sq_downsampled = rep(NA, num_peak_gene_pairs)
for (i in 1:num_peak_gene_pairs) {
    coaccessibility_sq_norm[i] = cor(norm_atac[peak[i],], norm_rna[gene[i],]) ** 2
    coaccessibility_sq_downsampled[i] = cor(downsampled_atac[peak[i],], downsampled_rna[gene[i],]) ** 2
    print(i)
}

coaccessibility_sq_norm = replace_na(coaccessibility_sq_norm, 0)
coaccessibility_sq_downsampled = replace_na(coaccessibility_sq_downsampled, 0)

# Take averages and average difference
avg_coacc_norm = mean(coaccessibility_sq_norm)
avg_coacc_downsampled = mean(coaccessibility_sq_downsampled)
avg_coacc_difference = avg_coacc_downsampled - avg_coacc_norm
out = data.frame(noise = avg_coacc_difference,
                 avg_coacc_norm = avg_coacc_norm,
                 avg_coacc_downsampled = avg_coacc_downsampled)

write.table(out, outfile, row.names = FALSE, quote = FALSE, sep = "\t")