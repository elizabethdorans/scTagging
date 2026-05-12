import pandas as pd
import numpy as np
import sys
import argparse
import os
from pathlib import Path

sys.path.insert(1, str(Path(__file__).resolve().parents[1]))
import functions as fn

parser = argparse.ArgumentParser()

parser.add_argument('--annot_file')
parser.add_argument('--annot_name')
parser.add_argument('--coaccessibility_file')
parser.add_argument('--stratified_coaccessibility_outfile')
parser.add_argument('--background_coaccess_file')
parser.add_argument('--feature', default = 'peak')

args = parser.parse_args()

annot_file = args.annot_file
annot_name = args.annot_name
coaccessibility_file = args.coaccessibility_file
stratified_coaccessibility_outfile = args.stratified_coaccessibility_outfile
background_coaccess_file = args.background_coaccess_file
feature = args.feature

if annot_name is None:
    raise ValueError("Please provide --annot_name matching the annotation column in --annot_file.")

# Read annotation
annot = pd.read_csv(annot_file, sep = '\t')
if annot_name not in annot.columns:
    raise ValueError("%s is not a column in %s." % (annot_name, annot_file))
all_peaks = annot[feature].tolist()

# Read coaccessibility
coaccess = pd.read_csv(coaccessibility_file, sep = '\t')
coaccess = coaccess[(coaccess['%s1' % feature].isin(all_peaks)) & (coaccess['%s2' % feature].isin(all_peaks))]

# Add peak-peak correlations for same-feature pairs
print(len(coaccess[coaccess['%s1' % feature] == coaccess['%s2' % feature]]))
if len(coaccess[coaccess['%s1' % feature] == coaccess['%s2' % feature]]) == 0:
    print("Adding self correlations!")
    peak_self_correlations = pd.DataFrame({'%s1' % feature: all_peaks,
                                                   '%s2' % feature: all_peaks,
                                                   'correlation': 1})
    
    coaccess = pd.concat([coaccess, peak_self_correlations]).reset_index(drop = True)
else:
    print("Self correlations already in data!")

# Annotate peak2 in peak-peak correlations
coaccess = coaccess.merge(annot, left_on = '%s2' % feature, right_on = feature)
coaccess = fn.square_column(coaccess, 'correlation')

# Subtract noise
coaccess_noise = pd.read_csv(background_coaccess_file,
                             sep = '\t')['noise'].iloc[0]
coaccess['correlation_squared'] = (coaccess['correlation_squared'] - coaccess_noise)

# Compute coaccessibility scores
print("Computing coaccessibility!")
coaccess['correlation_squared_annot'] = coaccess[annot_name] * coaccess['correlation_squared']
coaccess_scores = coaccess[['%s1' % feature, 'correlation_squared_annot']].groupby('%s1' % feature
    ).agg(sum).reset_index(
    ).rename(columns = {'%s1' % feature: feature, 'correlation_squared_annot': annot_name})

# Output scores
coaccess_scores.to_csv(stratified_coaccessibility_outfile, sep = "\t", index = False)
print("Saved scores to %s!" % stratified_coaccessibility_outfile)
