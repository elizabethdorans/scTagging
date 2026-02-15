# S-CASC
This folder contains code to implement stratified co-accessibility score regression (see Methods section of the manuscript for more details).

## Step 0: Compute peak co-activity scores

See 'Co-activity scores' in `../peak_scores` folder for code to compute co-activity scores.

## Step 1: Create peak categories/annotations

The script `create_annot.py` takes as input a bedfile of genomic regions and a bedfile of ATAC peaks, and generates a peak category (or "annotation") (denoting overlap between each ATAC peak and the set of genomic regions). This command should be run for each set of genomic regions of interest.

Example command: [~2 minutes, ~2G]

`Rscript create_annot.py --peaks_bedfile <peaks_bedfile> --annot_bedfile <annot_bedfile> --outfile <outfile>`

`<peaks_bedfile>`: Path to a bedfile containing ATAC peak coordinates.\
`<annot_bedfile>`: Path to a bedfile containing coordinates of genomic regions of interest.\
`<outfile>`: File to save output to.\
`<proportion_bp_overlap>`: If supplied, the "value" or category membership of each peak will be computed as the proportion of base pairs in the peak that overlap the set of genomic regions of interest. By default, the value for each peak is equal to 1 (if there is any overlap with the set of genomic regions) or 0 (if no overlap).
                    
Outputs: 

1) Peak category memberships at `<outfile>`.

## Step 2: Define number of _cis_ genes per peak

The script `number_nearby_genes.py` defines the number of nearby genes per peak (a covariate in S-CASC; see manuscript for method details). 

Example command: [~2 minutes, ~2G]

`Rscript number_nearby_genes.py --peaks_bedfile <peaks_bedfile> --coactivity_scores_file <coactivity_scores_file> --outfile <outfile>`

`<peaks_bedfile>`: Path to a bedfile containing ATAC peak coordinates.\
`<coactivity_scores_file>`: Path to file with co-activity scores.\
`<outfile>`: File to save output to.

Outputs: 

1) Number of _cis_ genes per peak at `<outfile>`.

## Step 3: Run S-CASC

See `stratified_coaccessibility_score_regression.ipynb` for code to run S-CASC.
