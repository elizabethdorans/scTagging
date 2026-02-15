# Peak-gene fine-mapping
This folder contains code to fine-map peaks to genes (see manuscript for additional method details).

## Step 0: Compute peak-gene links and save metacel matrices

See Step 0 in main folder for code to compute peak-gene links and save RNA and ATAC metacell matrices.

## Step 1: Fine-map peaks to genes

The script `peak_gene_finemapping.R` takes as input RNA and ATAC matrices and generates fine-mapped PIPs for each candidate peak-gene pair.

Example command: [~1 hour, ~5G]

`Rscript peak_gene_finemapping.R --coactivity_file <coactivity_file> --annot_bedfile <annot_bedfile> --outfile <outfile>`

`<coactivity_file>`: Path to peak-gene links file.\
`<rna_file>`: Path to RNA metacell-level matrix.\
`<atac_file>`: Path to ATAC metacell-level matrix.\
`<outfile>`: File to save output to.\
`<num_causal_peaks>`: [OPTIONAL] The maximum number of causal peaks per gene ('L' parameter in the susie() function). Default is 1. Maximum is 10.
`<weights_file>`: [OPTIONAL] A file with prior weights to be used in fine-mapping (default behavior is flat priors). Generated from stratified co-accessibility score regression (see the .ipynb notebook in the `S-CASC/` folder).\
`<link_universe_file>`: [OPTIONAL] 2-column .tsv file (columns = \['peak', 'gene'\]. If supplied, only the peak-gene pairs in the link universe file will be considered in fine-mapping (e.g. the universe of CRISPR-tested peak-gene pairs).\
`<blacklist_peaks_file>`: [OPTIONAL] Single-column text file with a list of peaks to exclude from fine-mapping.\
`<gene_set_file>`: [OPTIONAL] Single-column text file with a list of genes to include in fine-mapping.
                    
Outputs: 

1) Fine-mapped PIPs for each candidate peak-gene pair at `<outfile>`.
