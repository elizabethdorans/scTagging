import os
import pandas as pd
import argparse
from glob import glob
import pybedtools
import sys
sys.path.insert(1, "..")
import functions as fn

parser = argparse.ArgumentParser()

parser.add_argument('--peaks_bedfile')
parser.add_argument('--coactivity_scores_file')
parser.add_argument('--outfile')
parser.add_argument('--focal_feature', default = 'peak')

args = parser.parse_args()

peaks_bedfile = args.peaks_bedfile
coactivity_scores_file = args.coactivity_scores_file
outfile = args.outfile
focal_feature = args.focal_feature

# Create output directory
outdir = os.path.dirname(outfile)
if not os.path.exists(outdir):
    os.makedirs(outdir)
    
annot_name = 'number_nearby_genes'

# Load peaks
peaks_X_annot = pd.read_csv(peaks_bedfile, sep = '\t', header = None)
peaks_X_annot = peaks_X_annot.rename(columns = {3: 'peak'})[['peak']]
peaks_X_annot['peak'] = 'chr' + peaks_X_annot['peak']

# Load gene neighbors
neighbors = fn.read_scores(coactivity_scores_file)
neighbors = neighbors.rename(columns = {'cis_genes': annot_name})[[focal_feature, annot_name]]

# Merge with peaks, fillna(0) - peaks with no nearby genes
peaks_X_annot = peaks_X_annot.merge(neighbors, how = 'left')
peaks_X_annot[annot_name] = peaks_X_annot[annot_name].fillna(0).astype(int)

# Output annotation
peaks_X_annot.to_csv(outfile, sep = "\t", index = False)
print("Saved to %s!" % outfile)