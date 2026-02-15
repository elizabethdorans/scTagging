import pandas as pd
import numpy as np
import sys
import argparse
import os
sys.path.insert(1, "../..")
import functions as fn

parser = argparse.ArgumentParser()

parser.add_argument('--dataset')
parser.add_argument('--cell_type')
parser.add_argument('--annot_name')
parser.add_argument('--annot_dir')
parser.add_argument('--coaccess_scores_dir')
parser.add_argument('--background_coaccess_file')
parser.add_argument('--feature', default = 'peak')
parser.add_argument('--method', default = 'archr')

args = parser.parse_args()

print(args)

dataset = args.dataset
cell_type = args.cell_type
annot_name = args.annot_name
annot_dir = args.annot_dir
coaccess_scores_dir = args.coaccess_scores_dir
background_coaccess_file = args.background_coaccess_file
feature = args.feature
method = args.method

# Create output directory
if not os.path.exists(coaccess_scores_dir):
    os.makedirs(coaccess_scores_dir)
    
outfile = "%s/%s_X_%ss_%s_%s.tsv" % (coaccess_scores_dir, annot_name, feature, dataset, cell_type)
print(outfile)

# Read annotation
annot_file = "%s/%s_X_%ss_%s_%s.annot" % (annot_dir, annot_name, feature, dataset, cell_type)
annot = pd.read_csv(annot_file, sep = '\t')

all_peaks = annot[feature].tolist()

# Read coaccessibility
if feature == 'peak':
    coaccess = fn.read_coaccessibility(dataset, cell_type, method = method, subfolder = "../../")
    coaccess = fn.process_coaccessibility(coaccess, dataset, cell_type, path_to_main_folder = "../../")
elif feature == 'gene':
    coaccess = fn.read_coexpression(dataset, cell_type, method = method, subfolder = "../../")
    coaccess = fn.process_coexpression(coaccess)
    
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
coaccess_scores.to_csv(outfile, sep = "\t", index = False)
print("Saved scores to %s!" % outfile)