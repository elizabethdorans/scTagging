library(dplyr)
library(rhdf5)
library(ArchR)
library(parallel)
library(data.table)
library(argparse)

# Read in arguments
parser <- ArgumentParser()

parser$add_argument("--archr_proj_dir",
                    help = "name of folder containing ArchR project")
parser$add_argument("--out_dir",
                    help = "name of folder to send output")
parser$add_argument("--distance_threshold",
                    default = 1000000,
                    help = "distance threshold for peak-gene linking")

args <- parser$parse_args()

archr_proj_dir = args$archr_proj_dir
suffix = args$suffix
distance_threshold = as.numeric(args$distance_threshold)

if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

proj <- loadArchRProject(path = archr_proj_dir)

# Add iterative LSI for ATAC (peak matrix)
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(resolution = c(2), 
                       sampleCells = 10000, 
                       maxClusters = 6, 
                       n.start= 10),
  saveIterations = FALSE,
  useMatrix = "PeakMatrix", 
  depthCol = "nFrags",
  binarize = TRUE,
  name = "LSI_ATAC",
  force = TRUE
)

# Add iterative LSI for gene expression (RNA count matrix)
proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(resolution = c(2), 
                       sampleCells = 10000, 
                       maxClusters = 6, 
                       n.start= 10),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  force = TRUE
)

# Combined dims
proj <- addCombinedDims(proj, 
                        reducedDims = c("LSI_ATAC", "LSI_RNA"), 
                        name =  "LSI_Combined")

# Compute coaccessibility
proj <- addCoAccessibility(
  ArchRProj = proj,
  maxDist = distance_threshold,
  reducedDims = "LSI_Combined"
)

coacc <- getCoAccessibility(
    ArchRProj = proj,
    corCutOff = -1.5,
    resolution = 1,
    returnLoops = FALSE
)

peaks_idx = as.data.frame(metadata(coacc)[[1]])
peaks = paste0(peaks_idx$seqnames, "-", peaks_idx$start - 1, "-", peaks_idx$end)
coacc$peak1 = peaks[coacc$queryHits]
coacc$peak2 = peaks[coacc$subjectHits]

coacc = as.data.frame(coacc[c("peak1", "peak2", "correlation", "FDR")])

# Save coaccessibilities
write.table(coacc, sprintf("%s/coaccessibility_dist%s.tsv", out_dir, distance_threshold), sep = "\t", quote = FALSE)

# Compute and retrieve peak-gene links
proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  maxDist = distance_threshold,
  reducedDims = "LSI_Combined",
  useMatrix = "GeneExpressionMatrix",
  addEmpiricalPval = FALSE
)

p2g <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = -1.5,
    FDRCutOff = 1,
    resolution = 0,
    returnLoops = FALSE
)

peaks_idx = as.data.frame(metadata(p2g)$peakSet)
peaks = paste0(peaks_idx$seqnames, "-", peaks_idx$start - 1, "-", peaks_idx$end)
gene_idx = as.data.frame(metadata(p2g)$geneSet)
genes = gene_idx$name

peak_gene_links = data.frame(p2g)
peak_gene_links$peak = peaks[peak_gene_links$idxATAC]
peak_gene_links$gene = genes[peak_gene_links$idxRNA]

peak_gene_links = peak_gene_links[c("peak", "gene", "Correlation", "FDR")]

# Save peak-gene links
write.table(peak_gene_links, sprintf("%s/peak_gene_links_%s.tsv", out_dir, cell_type, distance_threshold), sep = "\t", quote = FALSE)

## Save metacell matrices

# Define output folders
rna_out_dir = sprintf("%s/RNA_metacell_matrices", out_dir)
if (!dir.exists(rna_out_dir)) {
    print("Creating rna_out_dir")
    dir.create(rna_out_dir, recursive = TRUE)
}

atac_out_dir = sprintf("%s/atac_metacell_matrices", out_dir)
if (!dir.exists(atac_out_dir)) {
    print("Creating atac_out_dir")
    dir.create(atac_out_dir, recursive = TRUE)
}

rna_outfile = sprintf("%s/RNA_metacell_matrix.rds", rna_out_dir)
atac_outfile = sprintf("%s/ATAC_metacell_matrix.rds", atac_out_dir)

# Save metacell x RNA matrix (using genes from the peak-gene matrices)
seRNA <- readRDS(metadata(p2g)$seRNA)
norm_RNA = assays(seRNA)[["RNA"]]
rownames(norm_RNA) = genes

saveRDS(norm_RNA, rna_outfile)

# Save metacell x atac matrix (using peaks from the peak-gene matrices)
seATAC <- readRDS(metadata(p2g)$seATAC)
norm_atac = assays(seATAC)[["ATAC"]]
rownames(norm_atac) = peaks

saveRDS(norm_atac, atac_outfile)

# Save ArchR project
saveArchRProject(ArchRProj = proj, outputDirectory = archr_proj_dir, overwrite = TRUE, load = FALSE)

print("Done!")