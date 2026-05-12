# Gene scores
This folder contains code for computing gene co-expression and gene co-activity scores for each gene in a single-cell RNA+ATAC multiome data set.

# Gene co-expression scores

## Step 1: Define cis gene-gene pairs

See `generate_proximal_gene_gene_pairs.py` for code to compute gene co-expression scores from gene-gene correlations and expected upward bias.

Example command: [~5 minutes, ~2G]

`python generate_proximal_gene_gene_pairs.py --rna_matrix_file <rna_matrix_file> --outfile <outfile>`

`<rna_matrix_file>`: Path to RNA metacell-level matrix.\
`<outfile>`: Path to file to save output to.\
`<gene_universe_file>`: [OPTIONAL] File with gene TSS coordinates (see gene_TSS.txt in the main folder for format). If running from this folder, then you do not need to supply it. If running outside this folder, then supply a revised path to ../gene_TSS.txt.\
`<max_distance>`: [OPTIONAL] Maximum gene-gene distance for computing gene-gene correlations (default is 1Mb).\
`<split_by_chromosome>`: [OPTIONAL] If supplied, will generate separate outfiles for each chromosome.

Outputs: 

1) Cis gene-gene pairs at `<outfile>`.

## Step 2: Compute gene-gene correlations

The script `gene_gene_coexpression.R` takes as input an RNA matrix for a single-cell RNA+ATAC multiome data set, candidate gene-gene pairs, and computes gene-gene correlations.

Example command: [~20 minutes, ~10G]

`Rscript gene_gene_coexpression.R --rna_matrix <rna_matrix> --gene_gene_pairs_file <gene_gene_pairs_file> --outfile <outfile>`

`<rna_matrix>`: Path to RNA metacell-level matrix.\
`<gene_gene_pairs_file>`: Path to file with candidate gene-gene pairs.\
`<outfile>`: Path to a file where output will be saved.

Outputs: 

1) Gene-gene correlations at `<outfile>`.

## Step 3: Compute expected upward bias in squared gene-gene correlations.

The script `get_background_coaccessibility_via_metacell_downsampling.R` (in the main folder) takes as input an RNA matrix for a single-cell RNA+ATAC multiome data set and computes the upward bias in squared gene-gene correlation expected due to finite sample size; see Methods section of the manuscript for further details).

Example command: [~1 hour, ~10G]

`Rscript get_background_coaccessibility_via_metacell_downsampling.R --norm_atac_file <norm_atac_file> --outfile <outfile>`

`<norm_atac_file>`: Path to RNA metacell-level matrix.\
`<outfile>`: Path to a file where output will be saved.\
`<blacklist_peaks_file>`: [OPTIONAL] Single-column text file with a list of genes to exclude.\
`<num_peak_pairs>`: [OPTIONAL] Number of random gene-gene pairs to use for the computation. Default is 100000 gene-gene pairs.
    
Outputs: 

1) Expected upward bias in squared gene-gene correlations due to finite sample size ('noise' column) at `<outfile>`.

## Step 4: Compute gene co-expression scores

See `gene_coexpression_scores.ipynb` for code to compute gene co-expression scores from gene-gene correlations and expected upward bias.

# Gene co-activity scores

## Step 1: Compute expected upward bias in squared peak-gene correlations.

See Co-activity scores: Step 1 in main folder for code to compute expected upward bias in squared peak-gene correlations. 

## Step 2: Compute gene co-activity scores

The script `bias_corrected_coactivity_scores.R` in the main folder takes as input peak-gene correlations and expected upward bias in squared correlations, and computes peak or gene co-activity scores.

Example command: [~10 min, ~5G]

`Rscript bias_corrected_coactivity_scores.R --uncorrected_file <uncorrected_file> --background_coacc_file <background_coacc_file> --outfile <outfile> --focal_feature gene`

`<uncorrected_file>`: Path to file containing uncorrected ArchR peak-gene correlations.\
`<background_coacc_file>`: Path to file containing expected upward bias in peak-gene squared correlations (output of step 1 above).\
`<outfile>`: Path to a file where output will be saved.\
`<blacklist_peaks_file>`: [OPTIONAL] Single-column text file with a list of peaks to exclude.\
`<gene_universe_file>`: [OPTIONAL] File with genes to include (see gene_TSS.txt in the main folder for format). If running from the peak_scores/ folder, then do not need to supply. If running outside the peak_scores/ folder, then supply revised path to ../gene_TSS.txt.\
`<focal_feature>`: Specify as "gene" to compute gene co-activity scores.

Outputs: 

1) Gene co-activity scores at `<outfile>`.
