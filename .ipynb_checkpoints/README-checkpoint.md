# scTagging
This repo contain code for implementing analyses in Dorans et al. medRxiv, "Distinguishing causal from tagging enhancers using single-cell multiome data." See the manuscript at XX for more details. See below (and analysis-specific folders) for steps to compute ATAC peak-level and gene-level scores, implement stratified co-accessibility score regression, and fine-map ATAC peaks to genes.

Clone the repository using the following command: 

`git clone https://github.com/elizabethdorans/scTagging.git`

Download required packages in `software_packages.txt`.

Submit jobs to your own remote cluster (e.g. using [sbatch](https://slurm.schedmd.com/sbatch.html)):
`cmd="\<command>"; sbatch --time=<time> --mem=<mem> ... --wrap="$cmd"`
Approximate memory and time requirements are given for computationally intensive tasks, but these will need to be adjusted for different data sets.

# Step 1: Compute peak-peak co-accessibility and peak-gene links using ArchR.

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

## Running E2G methods

See README.md files in method-specific folders for further steps in running each method!

## Post-processing for IGVF portal:

The script `postprocessing_for_IGVF_portal.R` takes as input peak-gene link predictions, restricts to a given gene universe, and produces a file with format appropriate for the IGVF portal.

**NOTE: In order to use the included gene universe file `CollapsedGeneBounds.hg38.bed`, run `postprocessing_for_IGVF_portal.R ` from this folder (`E2G_Method_Tutorials`) or specify the path to this file from your location as an argument to `--genes_file`.**

Example command: 

`Rscript reformatting_for_IGVF_portal.R --input_file <input_file> --output_file <output_file>  --genes_file IGVF_portal_genes_file.tsv --cell_type <cell_type> --sample_summary_short <sample_summary_short> --sample_term_id <sample_term_id> --method <method> --version <version> --metadata <metadata> --score_column <score_column> --score_type <score_type>`

<input_file>: A .tsv file containing E2G method predictions (includes at minimum columns 'peak', 'gene', and 'Score').\
<output_file>: A .tsv file containing containing E2G method predictions restricted to gene universe and reformatted for IGVF portal.\
<genes_file>: [DEFAULT CollapsedGeneBounds.hg38.bed] Gene universe file (inclues gene symbol in 'name' column and Ensembl ID in 'Ensembl_ID' column). The default gene universe file CollapsedGeneBounds.hg38.bed was obtained from https://github.com/EngreitzLab/CRISPR_comparison/blob/main/resources/genome_annotations/CollapsedGeneBounds.hg38.bed.\
<cell_type>: Cell type used to generate predictions (for 'CellType' column).\
<sample_summary_short>: Short string describing sample information (for header).\
<sample_term_id>: Sample ID (for header).\
\<method>: E2G method used to generate predictions (for header).\
\<version>: Version of E2G method used to generate predictions (for header).\
\<metadata>: IGVF data portal accession.\
\<score_column>: Name of column in input file to be renamed 'Score.'\
\<score_type>: Type of score (from options {positive_score, negative_score, p_value, adj_p_value, divergent, boolean}) (for header).\
                    
Outputs: 

1) Reformatted predictions at <output_file>

# peak_scores/

This folder contains code for computing co-accessibility and co-activity scores for each ATAC peak in a single-cell RNA+ATAC multiome data set.

# gene_scores/ 

This folder contains code for computing gene co-expression and gene co-activity scores for each gene in a single-cell RNA+ATAC multiome data set.

# S-CASC/

This folder contains code to implement stratified co-accessibility score regression.

# peak_gene_finemapping/

