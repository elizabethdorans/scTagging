import pandas as pd
import argparse
import numpy as np
import pybedtools
import pyreadr
import functions as fn
import os

parser = argparse.ArgumentParser()

parser.add_argument('--gene_list_file')
parser.add_argument('--rna_matrix_file')
parser.add_argument('--outfile')
parser.add_argument('--gene_universe_file',
                    default = "../gene_TSS.txt")
parser.add_argument('--max_distance', default = 1e+06)
parser.add_argument('--split_by_chromosome', action = 'store_true')

args = parser.parse_args()

gene_list_file = args.gene_list_file
rna_matrix_file = args.rna_matrix_file
gene_universe_file = args.gene_universe_file
outfile = args.outfile
max_distance = int(args.max_distance)
split_by_chromosome = args.split_by_chromosome

# Define set of genes of interest (measured genes with nonzero expression)
if gene_list_file != None:
    with open(gene_list_file) as f:
        measured_genes = f.read().splitlines()
elif rna_matrix_file != None:
    rna_mtx = pyreadr.read_r(rna_matrix_file)[None]
    rna_mtx = rna_mtx[rna_mtx.sum(axis = 1) > 0]
    measured_genes = rna_mtx.index.tolist()
else:
    print('Must supply either a list of genes or an RNA matrix!')

# Define distance windows around genes
gene_universe = pd.read_csv(gene_universe_file, sep = '\t')
gene_universe = gene_universe[gene_universe['gene'].isin(measured_genes)]
gene_universe['start'] = np.maximum(gene_universe['tss'] - max_distance, 0)
gene_universe['end'] = gene_universe['tss'] + max_distance
gene_universe = gene_universe[['chr', 'start', 'end', 'gene']]

# Intersect windows to find proximal gene-gene pairs
gene_pairs = fn.intersect_bed_dfs(gene_universe, gene_universe)
gene_pairs.columns = ['gene1', 'gene2']

# Save
if split_by_chromosome == True:
    gene_pairs = fn.get_gene_chrom(gene_pairs, gene_col = 'gene1')
    
    out_dir = '.'.join(outfile.split('.')[:-1])
    print(out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    chromosomes = gene_pairs['chr'].unique()
    for chromosome in chromosomes:
        print(chromosome)
        outfile_chrom = '%s/chrom%s.txt' % (out_dir, chromosome)
        gene_pairs_chrom = gene_pairs[gene_pairs['chr'] == chromosome]
        gene_pairs_chrom[['gene1', 'gene2']].to_csv(outfile_chrom, sep = '\t', index = False)
    
else:
    gene_pairs.to_csv(outfile, sep = '\t', index = False)
