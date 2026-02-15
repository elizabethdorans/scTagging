import os
import pandas as pd
import argparse
from glob import glob
import pybedtools
import sys
sys.path.insert(1, "../..")
import functions as fn

parser = argparse.ArgumentParser()

parser.add_argument('--peaks_bedfile')
parser.add_argument('--annot_bedfile')
parser.add_argument('--annot_name')
parser.add_argument('--outfile')
parser.add_argument('--proportion_bp_overlap',
                   action = 'store_true')

args = parser.parse_args()

peaks_bedfile = args.peaks_bedfile
annot_bedfile = args.annot_bedfile
annot_name = args.annot_name
outfile = args.outfile
proportion_bp_overlap = args.proportion_bp_overlap

# Create output directory
outdir = os.path.dirname(outfile)
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Load peaks
peaks_df = fn.bedfile_to_bedtool(peaks_bedfile)
peaks_bed = fn.df_to_bedtool(peaks_df)
print("Read bedfile for %s peaks!" % len(peaks_df))

# Load annot
annot_df = fn.bedfile_to_bedtool(annot_bedfile)
annot_bed = fn.df_to_bedtool(annot_df)
print("Read bedfile for annotation %s: %s!" % (annot_name, annot_bedfile))

# Intersect annot with peaks
if proportion_bp_overlap == True:
    print("Computing proportion of basepair overlap!")
    peaks_X_annot = peaks_bed.intersect(annot_bed, wao = True).to_dataframe().iloc[:, [0, 1, 2, 6]].rename(columns = {'thickStart': 'bp_overlap'})
    peaks_X_annot = peaks_X_annot.groupby(['chrom', 'start', 'end']).agg(sum).reset_index()
    peaks_X_annot['peak_length'] = peaks_X_annot['end'] - peaks_X_annot['start']
    peaks_X_annot[annot_name] = peaks_X_annot['bp_overlap'] / peaks_X_annot['peak_length']
    
else:
    ("Detecting binary overlap!")
    peaks_X_annot = peaks_bed.intersect(annot_bed, c = True).to_dataframe()
    # Define binary annot for any overlap
    peaks_X_annot['name'] = 1 * (peaks_X_annot['name'] >= 1)
    peaks_X_annot = peaks_X_annot.rename(columns = {'name': annot_name})

# Clean up columns
peaks_X_annot = fn.add_peak_column(peaks_X_annot, chrom_col = 'chrom')
peaks_X_annot = peaks_X_annot[['peak', annot_name]]

# Output annotation size
annot_size = peaks_X_annot[annot_name].mean()
print("Annotation size:", annot_size)
with open("%s.size" % outfile, 'w') as file:
    file.write('%f' % annot_size)

# Output annotation
peaks_X_annot.to_csv(outfile, sep = "\t", index = False)
print("Saved to %s!" % outfile)