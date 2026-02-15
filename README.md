# scTagging
This repo contain code for implementing analyses in Dorans et al. medRxiv, "Distinguishing causal from tagging enhancers using single-cell multiome data." See below (and analysis-specific folders) for steps to compute ATAC peak-level and gene-level scores, implement stratified co-accessibility score regression, and fine-map ATAC peaks to genes.

Clone the repository using the following command: 

`git clone https://github.com/elizabethdorans/scTagging.git`

Submit jobs to your own remote cluster (e.g. using [sbatch](https://slurm.schedmd.com/sbatch.html)):
`cmd="<command>"; sbatch --time=<time> --mem=<mem> ... --wrap="$cmd"`
Approximate memory and time requirements are given for computationally intensive tasks, but these will need to be adjusted for different data sets.

# Step 0: Compute peak-peak co-accessibility and peak-gene links using ArchR.

The script `run_archr_peak_gene_linking_coaccessibility_metacells.R` takes as input an ArchR project for a single-cell RNA+ATAC multiome data set, computes peak-peak correlations for all peak-peak pairs with distance <1Mb, computes peak-gene correlations for all peak-gene pairs with distance <1Mb, and saves RNA and ATAC metacell-level matrices for downstream analyses.

Example command: [~1 hour, ~20G]

`Rscript run_archr_peak_gene_linking_coaccessibility_metacells.R --archr_proj_dir <archr_proj_dir> --out_dir <out_dir>`

<archr_proj_dir>: Path to a folder containing an ArchR project with ATAC and RNA data from a single-cell RNA+ATAC multiome data set (same cell barcodes in both RNA and ATAC matrices).\
<out_dir>: Path to a folder where outputs will be saved (peak-peak co-accessibilities, peak-gene links, and RNA and ATAC metacell-level matrices).\
<distance_threshold>: [OPTIONAL] Maximum peak-peak and peak-gene distance for computing peak-peak and peak-gene correlations (default is 1Mb).
                    
Outputs: 

1) Peak-peak correlations at <out_dir>/coaccessibility_dist<distance_threshold>.tsv
2) Peak-gene correlations at <out_dir>/peak_gene_links_dist<distance_threshold>.tsv
3) RNA metacell matrix at <out_dir>/RNA_metacell_matrices/RNA_metacell_matrix.rds
3) ATAC metacell matrix at <out_dir>/ATAC_metacell_matrices/ATAC_metacell_matrix.rds

# peak_scores/

This folder contains code for computing co-accessibility and co-activity scores for each ATAC peak in a single-cell RNA+ATAC multiome data set.

# gene_scores/ 

This folder contains code for computing gene co-expression and gene co-activity scores for each gene in a single-cell RNA+ATAC multiome data set.

# S-CASC/

This folder contains code to implement stratified co-accessibility score regression.

# peak_gene_finemapping/

This folder contains code to implement peak-gene fine-mapping.

