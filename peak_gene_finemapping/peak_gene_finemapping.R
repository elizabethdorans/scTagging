library(data.table)
library(argparse)
library(dplyr)
library('susieR')

set.seed(511)

parser <- ArgumentParser()

parser$add_argument("--link_universe_file",
                   help = "File containing the set of candidate peak-gene links to consider (coactivity_file will be intersected with these links)")
parser$add_argument("--coactivity_file",
                   help = "File containing summary statistics (linking scores) for each candidate peak-gene links")
parser$add_argument("--rna_file",
                   help = "Matrix of RNA measurements")
parser$add_argument("--atac_file",
                   help = "Matrix of ATAC measurements")
parser$add_argument("--num_causal_peaks",
                   default = 1,
                   help = "Number of causal signals to allow in SuSiE algorithm")
parser$add_argument("--gene_set_file",
                   help = "Single-column text file with list of genes to restrict to (only these genes will be fine-mapped)")
parser$add_argument("--blacklist_peaks_file",
                   help = "Single-column text file with list of peaks to exclude (these peaks will not be considered in fine-mapping)")
parser$add_argument("--weights_file",
                   help = "file with columns 'peak', 'weight'")
parser$add_argument("--outfile")

args <- parser$parse_args()

link_universe_file = args$link_universe_file
coactivity_file = args$coactivity_file
rna_file = args$rna_file
atac_file = args$atac_file
num_causal_peaks = as.numeric(args$num_causal_peaks)
gene_set_file = args$gene_set_file
blacklist_peaks_file = args$blacklist_peaks_file
weights_file = args$weights_file
outfile = args$outfile

# Define function for creating empty results dataframe
create_empty_dataframe <- function() {
    res <- data.frame(peak = character(), gene = numeric(), pip = character())
    return(res)
}

print(sprintf("Assuming %s causal peaks", num_causal_peaks))

out_dir = dirname(outfile)
if (!dir.exists(out_dir)) {
    print(sprintf("Creating directory"))
    dir.create(out_dir, recursive = TRUE)
    }

# Read in peak_gene_links (restrict to first 3 columns and name the third column "Correlation")
peak_gene_links = read.table(coactivity_file, header = TRUE)
peak_gene_links = peak_gene_links[,1:3]
colnames(peak_gene_links) = c("peak", "gene", "Correlation")

if (!is.null(link_universe_file)) {
    link_universe = read.table(link_universe_file, header = TRUE)
    print(sprintf("Restricting to candidate link universe! From %s pairs...", nrow(peak_gene_links)))
    peak_gene_links <- merge(peak_gene_links, link_universe)
    print(sprintf("...to %s pairs!", nrow(peak_gene_links)))
}

if (!is.null(gene_set_file)) {
    genes = read.table(gene_set_file)$V1
    print(sprintf("Restricting to supplied genes only! From %s pairs...", nrow(peak_gene_links)))
    peak_gene_links <- peak_gene_links[peak_gene_links$gene %in% genes,]
    print(sprintf("...to %s pairs!", nrow(peak_gene_links)))
}

if (!is.null(blacklist_peaks_file)) {
    blacklist_peaks = read.table(blacklist_peaks_file, sep = "\t")$V1
    blacklist_peaks = paste0("chr", blacklist_peaks)
    print(sprintf("Restricting to non-blacklisted peaks only! From %s pairs...", nrow(peak_gene_links)))
    peak_gene_links = peak_gene_links[!(peak_gene_links$peak %in% blacklist_peaks),]
    print(sprintf("...to %s pairs!", nrow(peak_gene_links)))
}

if (!is.null(weights_file)) {
    print("Executing weighted regression!")
    weights = read.csv(weights_file, sep = '\t')[1:2]
    colnames(weights) = c('peak', 'weight')
    peak_gene_links = merge(peak_gene_links, weights)
}

# Read in RNA matrix
print("Reading in RNA matrix")
rna = readRDS(rna_file)
rna = rna[Matrix::rowSums(rna) > 0,]
sample_size = dim(rna)[2]

# Read in ATAC matrix
print("Reading in ATAC matrix")
atac = readRDS(atac_file)
atac = atac[Matrix::rowSums(atac) > 0,]

# Subset to candidate peak-gene links
print("Defining focal genes")
genes = intersect(unique(peak_gene_links$gene), rownames(rna))
print("Defining focal peaks")
peaks = intersect(unique(peak_gene_links$peak), rownames(atac))

rna = rna[genes, , drop = FALSE]
atac = atac[peaks, , drop = FALSE]

# Scale RNA
rna_scaled = t(as.matrix(rna))
rna_scaled = scale(rna_scaled, center = TRUE, scale = TRUE)

peak_gene_links = peak_gene_links[peak_gene_links$gene %in% genes,]
peak_gene_links = peak_gene_links[peak_gene_links$peak %in% peaks,]
print(sprintf("subsetted to %s links", nrow(peak_gene_links)))

# Write actual candidate links used in fine-mapping
print(sprintf("%s_candidate_links", outfile))
write.table(peak_gene_links, sprintf("%s_candidate_links", outfile), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Fine-mapping
res = data.frame()

for (gene in genes) {
    print(gene)
    
    # Define focal gene and peaks
    focal_peaks = peak_gene_links[peak_gene_links$gene == gene,]$peak
    
    if (length(focal_peaks) > 1) {

        # Set focal RNA and ATAC data
        y = rna_scaled[,gene]
        X = t(as.matrix(atac[focal_peaks,]))
        X = scale(X,center = TRUE,scale = TRUE)
        
        focal_gene_L = min(num_causal_peaks, ncol(X))
        
        # Get prior weights (if applicable)
        if (!is.null(weights_file)) {
            print("executing weighted regression")
            focal_weights = peak_gene_links[peak_gene_links$gene == gene,]$weight
            total_weight = sum(focal_weights)
            prior_weights = focal_weights / total_weight
        } else {
            prior_weights = NULL
        }
        
        # Fit SuSiE model using cell-level data
        fitted <- susie(X, 
                        y, 
                        L = focal_gene_L,
                        prior_weights = prior_weights)
        
        # Extract PIPs
        pips <- data.frame(peak = focal_peaks, gene = gene, pip = fitted$pip)
        
        res <- rbind(res, pips)

    }
    
}

if (nrow(res) == 0) {
    res <- create_empty_dataframe()
}

print(head(res))
write.table(res, outfile, sep = "\t", row.names = FALSE, quote = FALSE)