# Peak scores
This folder contains code for computing co-accessibility and co-activity scores for each ATAC peak in a single-cell RNA+ATAC multiome data set.

# Co-accessibility scores

## Step 1: Compute expected upward bias in squared peak-peak correlations.

The script `get_background_coaccessibility_via_metacell_downsampling.R` (in the main folder) takes as input an ATAC matrix for a single-cell RNA+ATAC multiome data set and computes the upward bias in peak-peak squared correlation expected due to finite sample size; see Methods section of the manuscript for further details).

Example command: [~1 hour, ~10G]

`Rscript get_background_coaccessibility_via_metacell_downsampling.R --norm_atac_file <norm_atac_file> --outfile <outfile>`

<norm_atac_file>: Path to ATAC metacell-level matrix.\
<\outfile>: Path to a file where output will be saved.\
<blacklist_peaks_file>: [OPTIONAL] Single-column text file with a list of peaks to exclude.\
<num_peak_pairs>: [OPTIONAL] Number of random peak-peak pairs to use for the computation. Default is 100000 peak-peak pairs.\
    
Outputs: 

1) Expected upward bias in squared peak-peak correlations due to finite sample size ('noise' column) at <\outfile>.

## Step 2: Compute co-accessibility scores

The script `bias_corrected_coaccessibility_scores.R` takes as input peak-peak correlations and expected upward bias in squared correlations, and computes co-accessibility scores.

Example command: [~20 min, ~10G]

`Rscript bias_corrected_coaccessibility_scores.R --uncorrected_file <uncorrected_file> --background_coacc_file <background_coacc_file> --outfile <outfile>`

<uncorrected_file>: Path to file containing uncorrected ArchR peak-peak correlations.\
<background_coacc_file>: Path to file containing expected upward bias in peak-peak squared correlations (output of step 1 above).\
<\outfile>: Path to a file where output will be saved.\
<blacklist_peaks_file>: [OPTIONAL] Single-column text file with a list of peaks to exclude.\

Outputs: 

1) Co-accessibility scores at <\outfile>.

# Co-activity scores

## Step 1: Compute expected upward bias in squared peak-gene correlations.

The script `get_background_coactivity_via_metacell_downsampling.R` (in the main folder) takes as input RNA and ATAC matrices for a single-cell RNA+ATAC multiome data set and computes the upward bias in squared peak-gene correlation expected due to finite sample size; see Methods section of the manuscript for further details).

Example command: [~1 hour, ~10G]

`Rscript get_background_coactivity_via_metacell_downsampling.R --norm_atac_file <norm_atac_file> --norm_rna_file <norm_rna_file> --outfile <outfile>`

<norm_atac_file>: Path to ATAC metacell-level matrix.\
<norm_rna_file>: Path to ATAC metacell-level matrix.\
<\outfile>: Path to a file where output will be saved.\
<blacklist_peaks_file>: [OPTIONAL] Single-column text file with a list of peaks to exclude.\
<gene_universe_file>: [OPTIONAL] File with genes to include (see gene_TSS.txt in the main folder for format). If running from the peak_scores/ folder, then do not need to supply. If running outside the peak_scores/ folder, then supply revised path to ../gene_TSS.txt.\
<num_peak_gene_pairs>: [OPTIONAL] Number of random peak-peak pairs to use for the computation. Default is 100000 peak-peak pairs.\
    
Outputs: 

1) Expected upward bias in squared peak-gene correlations due to finite sample size ('noise' column) at <\outfile>.

## Step 2: Compute co-activity scores

The script `bias_corrected_coactivity_scores.R` in the main folder takes as input peak-gene correlations and expected upward bias in squared correlations, and computes co-activity scores.

Example command: [~10 min, ~5G]

`Rscript bias_corrected_coactivity_scores.R --uncorrected_file <uncorrected_file> --background_coacc_file <background_coacc_file> --outfile <outfile> --focal_feature peak`

<uncorrected_file>: Path to file containing uncorrected ArchR peak-peak correlations.\
<background_coacc_file>: Path to file containing expected upward bias in peak-peak squared correlations (output of step 1 above).\
<\outfile>: Path to a file where output will be saved.\
<blacklist_peaks_file>: [OPTIONAL] Single-column text file with a list of peaks to exclude.\
<gene_universe_file>: [OPTIONAL] File with genes to include (see gene_TSS.txt in the main folder for format). If running from the peak_scores/ folder, then do not need to supply. If running outside the peak_scores/ folder, then supply revised path to ../gene_TSS.txt.\
<focal_feature>: Specified as "peak" (by default) to compute peak co-activity scores.\

Outputs: 

1) Co-activity scores at <\outfile>.