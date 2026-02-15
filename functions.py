import pandas as pd
import numpy as np
import scipy
from matplotlib import pyplot as plt
import pybedtools
from glob import glob
import statsmodels.api as sm
from statsmodels.stats.meta_analysis import combine_effects
import os
import seaborn as sns
import pyreadr
from matplotlib.patches import Rectangle
from scipy import sparse

def sig_star(p):
    if p < 0.05:
        if p < 0.001:
            star = "***"
        elif p < 0.01:
            star = "**"
        else:
            star = "*"
    else:
        star = ""
    return star


def get_coaccess_file(dataset, cell_type, method = "archr", subfolder = "", folder_suffix = "", file_suffix = "_dist1e+06_k100_knniter500_overlapcutoff0.8", addl_file_suffix = ""):
    if method == "archr":
        return "%sarchr_%s_%s%s/coaccessibility%s%s.tsv" % (subfolder, dataset, cell_type, folder_suffix, file_suffix, addl_file_suffix)
    elif method == "signac":
        return "%ssignac_%s_%s/coaccessibilities_distance_1e+06.tsv" % (subfolder, dataset, cell_type)
    else:
        raise Exception("Method must be 'archr' or 'signac'")


def get_coaccess_scores_file(dataset, cell_type, subfolder = "", folder_suffix = "", file_suffix = "", addl_file_suffix = "_corrected_via_downsampling_metacells"):
    return "%sarchr_%s_%s%s/coaccessibility_scores%s%s.tsv" % (subfolder, dataset, cell_type, folder_suffix, file_suffix, addl_file_suffix)


def get_signac_coaccess_scores_file(dataset, cell_type, subfolder = "", folder_suffix = "", file_suffix = "_corrected_via_downsampling_metacells", addl_file_suffix = ""):
    return "%ssignac_%s_%s%s/coaccessibility_scores%s%s.tsv" % (subfolder, dataset, cell_type, folder_suffix, file_suffix, addl_file_suffix)

def get_coactivity_file(dataset, cell_type, method = "archr", subfolder = "", folder_suffix = "", file_suffix = "_dist_1e+06"):
    if method == "archr":
        return "%s%s_%s_%s%s/peak_gene_links_%s%s.tsv" % (subfolder,  method, dataset, cell_type, folder_suffix, cell_type, file_suffix)
    elif method == "signac":
        return "signac_%s_%s/pgl/signac_peak_gene_links.tsv" % (dataset, cell_type)


def get_coactivity_scores_file(dataset, cell_type, method = "archr", subfolder = "", folder_suffix = "", file_suffix = "", addl_file_suffix = "_corrected_via_downsampling_metacells"):
    return "%s%s_%s_%s%s/coactivity_scores%s%s.tsv" % (subfolder, method, dataset, cell_type, folder_suffix, file_suffix, addl_file_suffix)

def get_coactivity_scores_by_gene_file(dataset, cell_type, subfolder = "", folder_suffix = "", file_suffix = ""):
    return "%sarchr_%s_%s%s/coactivity_scores_by_gene%s.tsv" % (subfolder, dataset, cell_type, folder_suffix, file_suffix)

def get_coexpression_file(dataset, cell_type, subfolder = "", method = 'archr', file_suffix = "_dist_1e+06_k100_knniter500_overlapcutoff0.8"):
    return "%s%s_%s_%s/coexpression%s.tsv" % (subfolder, method, dataset, cell_type, file_suffix)

def get_coexpression_scores_file(dataset, cell_type, method = "archr", file_suffix = "_dist_1e+06_corrected_via_downsampling_metacells", addl_file_suffix = ""):
    return "%s_%s_%s/coexpression_scores%s%s.tsv" % (method, dataset, cell_type, file_suffix, addl_file_suffix)

def get_gene_coactivity_scores_file(dataset, cell_type, method = 'archr', file_suffix = "_dist_1e+06_corrected_via_downsampling_metacells", subfolder = "", addl_file_suffix = ""):
    return "%s%s_%s_%s/gene_coactivity_scores%s%s.tsv" % (subfolder, method, dataset, cell_type, file_suffix, addl_file_suffix)


def remove_blacklist_peaks(df, blacklist_peaks, col = "peak"):
    df = df[~df[col].isin(blacklist_peaks)]
    return df

def get_blacklist_peaks_file(dataset, cell_type, path_to_main_folder = ""):
    return "%sarchr_%s_%s/peaks_X_promoters.txt" % (path_to_main_folder, dataset, cell_type)

def get_blacklist_peaks(blacklist_peaks_file):
    with open(blacklist_peaks_file) as bpf:
        blacklist = [l.rstrip() for l in bpf.readlines()]
    blacklist = ["chr" + peak for peak in blacklist]
    return blacklist

def remove_sex_chromosome_peaks(df, col = "peak"):
    return df[~(df[col].str.startswith("chrX")) & ~(df[col].str.startswith("chrY"))]


def read_correlations(file):
    return pd.read_csv(file, sep = "\t")


def read_scores(file, rename_suffix = None, scored_unit = 'peak'):
    scores = pd.read_csv(file, sep = "\t")
    if rename_suffix != None:
        scores = scores.rename(columns = dict(zip(
            scores.columns, 
            [scored_unit] + ["%s_%s" % (c, rename_suffix) for c in scores.columns if c != scored_unit]))
        )
    return scores

def read_and_process_scores(dataset, cell_type, metric = "coaccess", rename_suffix = "coaccess"):
    if metric == "coaccess" or metric == "coaccessibility":
        scores_file = get_coaccess_scores_file(dataset, cell_type)
    elif metric == "coactivity":
        scores_file = get_coactivity_scores_file(dataset, cell_type)
    scores = read_scores(scores_file, rename_suffix = rename_suffix)
    blacklist_peaks_file = get_blacklist_peaks_file(dataset, cell_type)
    blacklist_peaks = get_blacklist_peaks(blacklist_peaks_file)
    scores = remove_blacklist_peaks(scores, blacklist_peaks)
    scores = remove_sex_chromosome_peaks(scores)
    return scores


def coaccess_coactivity_plot(merged, coacc_col, coact_col, num_bins_pre_dup = None, num_bins = None, num_peaks_per_bin = None, 
                             outfile = None, return_df = False, color_col = None, return_corr_only = False, unit = "peaks", report_corr = True, report_unbinned_corr = True, group_by_coacc_col = False,
                            best_fit_line = False, figsize = (6, 4), axis = None, legend_fontsize = 9, paper_figure = False, across_text = True, text_at_top = False, legend = True):
    merged = merged.sort_values(by = coacc_col)
    if group_by_coacc_col == True:
        merged_agg = merged[['peak', coacc_col, coact_col]].groupby(coacc_col).agg(np.mean).reset_index()
    else:
        if num_bins == None:
            num_bins = int(np.round(len(merged) / num_peaks_per_bin))
        if num_bins != len(merged):
            merged["bin"] = pd.qcut(merged[coacc_col], num_bins, duplicates = "drop")
            merged_agg = merged.groupby("bin").aggregate("mean", numeric_only = True).reset_index()
            print("Number of bins: %s" % len(merged_agg))
        else:
            merged_agg = merged
    x, y = merged_agg[coacc_col], merged_agg[coact_col]
    if color_col != None:
        color = merged_agg[color_col]
    else:
        color = "#7294d4"
    if return_corr_only == True:
        return merged_agg[[coacc_col, coact_col]].corr()[coacc_col].iloc[1]
    else:
        if axis == None:
            fig, axis = plt.subplots(figsize = figsize)
        axis.scatter(x, y, c = color, lw = 2)
        if report_corr == True:
            y_coord = (axis.get_ylim()[1]) if text_at_top == True else np.min(y)
            #y_coord = np.max(y) if text_at_top == True else np.min(y)
            
            if report_unbinned_corr == True:
                suffix = (" across individual %s" % unit) if across_text == True else ""
                axis.text(np.max(x), y_coord,
                     r"${\it{r} = %.2f}$%s" % (np.round(merged[[coacc_col, coact_col]].corr()[coacc_col].iloc[1], 2), suffix),
                    #transform=ax.transAxes,
                    fontsize = 11, color = "black",
                        ha = 'right')
            else:
                axis.text(np.max(x), y_coord,
                     r"${\it{r} = %.2f}$" % np.round(merged_agg[[coacc_col, coact_col]].corr()[coacc_col].iloc[1], 2),
                    #transform=ax.transAxes,
                    fontsize = 11, color = "black", fontweight = "bold",
                        ha = 'right')
        axis.set_ylabel("Co-activity score")
        axis.set_xlabel("Co-accessibility score")
        if legend == True:
            if group_by_coacc_col == True:
                axis.legend(["Bin (by %s)" % (coacc_col)])
            elif num_peaks_per_bin != None:
                axis.legend(["Bin (%s %s)" % (num_peaks_per_bin, unit)], fontsize = legend_fontsize)
            elif num_bins != None:
                number_units = int(np.floor(len(merged) / len(merged_agg)))
                if number_units > 10000: 
                    number_units = "%sk" % (round(number_units, -3) // 1000)
                axis.legend(["Bin (%s bins, %s %s per bin)" % (len(merged_agg), number_units, unit)], fontsize = legend_fontsize, handletextpad = 0, loc = 'upper left')
        current_ymax = merged_agg[coact_col].max()
        axis.set_ylim(axis.get_ylim()[0], current_ymax * 1.2)
        if best_fit_line == True:
            axis.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), color = "black")
        if color_col != None:
            colorbar = plt.colorbar()
            colorbar.set_label(color_col)
        if outfile != None:
            plt.savefig(outfile, fmt = "pdf")
        if return_df == True:
            return merged_agg


def get_mean(merged, metric, condition):
    return merged["%s_%s" % (metric, condition)].mean()


def get_loco_df(merged, chrom):
    return merged[~merged["peak"].str.startswith("chr%s" % chrom)]


def get_standard_errors(jk1, jk2):
    n = len(jk1) - 1
    jk_diff = jk1 - jk2
    se1 = jk1.std() * np.sqrt(n)
    se2 = jk2.std() * np.sqrt(n)
    se_diff = jk_diff.std() * np.sqrt(n)
    return se1, se2, se_diff


def get_p_from_est_and_se(est, se, null = 0, return_z = False):
    est = est - null
    z = est / se
    p = 2 * scipy.stats.norm.sf(abs(z))
    if return_z == True:
        return z
    else:
        return p


def compute_mean_jackknife_se_p(merged, metric,
                                h0_mean = 0,
                                print_summary = False):
    mean = merged[metric].mean()
    diff = mean - h0_mean
    
    # Jackknife
    jk = np.empty(22)
    for i in range(22):
        chrom = i + 1
        tmp = get_loco_df(merged, chrom)
        jk[i] = tmp[metric].mean()

    n = len(jk) - 1
    se = jk.std() * np.sqrt(n)
    
    jk_diff = jk - h0_mean
    
    p_diff = get_p_from_est_and_se(diff, se)
    
    if print_summary == True:
        print("mean = %s, SE = %s, p (diff. from %s) = %s" % (mean, se, h0_mean, p_diff))

    return mean, se, p_diff


def get_correlation(merged, quantity1, quantity2, condition):
    return merged[["%s_%s" % (quantity1, condition), "%s_%s" % (quantity2, condition)]].corr().iloc[0,1]


def compare_mean_across_conditions(merged, metric, condition1, condition2,
                                  print_summary = False):
    mean1 = get_mean(merged, metric, condition1)
    mean2 = get_mean(merged, metric, condition2)
    diff = mean1 - mean2
    
    # Jackknife
    jk1 = np.empty(22)
    jk2 = np.empty(22)
    for i in range(22):
        chrom = i + 1
        tmp = get_loco_df(merged, chrom)
        jk1[i] = get_mean(tmp, metric, condition1)
        jk2[i] = get_mean(tmp, metric, condition2)

    se1, se2, se_diff = get_standard_errors(jk1, jk2)
    p_diff = get_p_from_est_and_se(diff, se_diff)
    
    if print_summary == True:
        print("mean1 = %s, mean2 = %s, p_diff = %s" % (mean1, mean2, p_diff))

    return mean1, se1, mean2, se2, diff, p_diff


def get_correlation(merged, quantity1, quantity2, condition):
    return merged[["%s_%s" % (quantity1, condition), "%s_%s" % (quantity2, condition)]].corr().iloc[0,1]


def compare_correlation_across_conditions(merged, quantity1, quantity2, condition1, condition2,
                                          print_summary = False):
    corr1 = get_correlation(merged, quantity1, quantity2, condition1)
    corr2 = get_correlation(merged, quantity1, quantity2, condition2)
    diff = corr1 - corr2

    # Jackknife
    jk1 = np.empty(22)
    jk2 = np.empty(22)
    for i in range(22):
        chrom = i + 1
        tmp = get_loco_df(merged, chrom)
        jk1[i] = get_correlation(tmp, quantity1, quantity2, condition1)
        jk2[i] = get_correlation(tmp, quantity1, quantity2, condition2)

    se1, se2, se_diff = get_standard_errors(jk1, jk2)
    p_diff = get_p_from_est_and_se(diff, se_diff)
    
    if print_summary == True:
        print("corr1 = %s, corr2 = %s, p_diff = %s" % (corr1, corr2, p_diff))
    
    return corr1, se1, corr2, se2, diff, p_diff


def plot_across_conditions(res, ylabel, condition1_label, condition2_label, title,
                          print_summary = False):
    point_est1, se1, point_est2, se2, p_diff = res[0], res[1], res[2], res[3], res[5]
    fig, ax = plt.subplots(figsize = (3,4))
    plt.bar([condition1_label, condition2_label], [point_est1, point_est2], color = "#7294d4")
    plt.errorbar([condition1_label, condition2_label], [point_est1, point_est2], ls = "",
                  xerr = None, yerr = [se1, se2], 
                  color = "black", capsize = 4)
    plt.xlabel("Condition")
    plt.ylabel(ylabel)
    plt.title(title)
    if print_summary == True:
        print("Diff = %s, p-value for difference = %s" % (point_est1 - point_est2, p_diff))
        
        
def plot_across_conditions_multi(point_estimates, standard_errors, labels, title = "", xlabel = "Condition", ylabel = "Estimate", figsize = (3,4), xtick_rotation = 0, xtick_ha = "center", pvals = None, null = 0, v_space = 0, box_index = None, axis = None):
    if axis == None:
        fig, ax = plt.subplots(figsize = figsize)
    else:
        ax = axis
    ax.bar(labels, [e - null for e in point_estimates], bottom = null, color = "#7294d4")
    if np.max(standard_errors) > 0:
        ax.errorbar(labels, point_estimates, ls = "",
                      xerr = None, yerr = standard_errors,
                      color = "black", capsize = 4)
    ax.set_xticks(labels, labels, rotation = xtick_rotation, ha = xtick_ha)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if pvals != None:
        for i in range(len(labels)):
            if point_estimates[i] > null:
                ax.text(labels[i], point_estimates[i] + standard_errors[i], sig_star(pvals[i]), color = "green", ha = "center", va = "bottom")
            else:
                ax.text(labels[i], point_estimates[i] - standard_errors[i], sig_star(pvals[i]), color = "red", ha = "center", va = "top")
    ax.set_xlim(-0.75, len(labels) - 0.5)
    ax.set_ylim(min(ax.get_ylim()[0], min([e - s for e, s in zip(point_estimates, standard_errors)]) - v_space),
             max(ax.get_ylim()[1], max([e + s for e, s in zip(point_estimates, standard_errors)])) + v_space)
    if box_index is not None:
        ax.add_patch(Rectangle((box_index - 0.5, 0), 1.0, 
                                1.025 * (point_estimates[box_index] + standard_errors[box_index]),
                               fill=False, edgecolor='black', lw=3))
    if axis == None:
        return fig, ax


def plot_correlation_across_conditions(merged, quantity1, quantity2, condition1, condition2,
                                       condition1_label = None, condition2_label = None, title = None):
    res = compare_correlation_across_conditions(merged, quantity1, quantity2, condition1, condition2)
    ylabel = "Value"
    if condition1_label == None or condition2_label == None:
        condition1_label, condition2_label = condition1, condition2
    if title == None:
        title = "corr(%s, %s)" % (quantity1, quantity2)
    plot_across_conditions(res, ylabel, condition1_label, condition2_label, title)


def plot_mean_across_conditions(merged, metric, condition1, condition2,
                                condition1_label = None, condition2_label = None, title = None):
    res = compare_mean_across_conditions(merged, metric, condition1, condition2)
    ylabel = "Average"
    if condition1_label == None or condition2_label == None:
        condition1_label, condition2_label = condition1, condition2
    if title == None:
        title = "Mean %s" % metric
    plot_across_conditions(res, ylabel, condition1_label, condition2_label, title)


def add_bed_columns(data, element_column = "peak", element_gene_links = False, chr_prefix = True):
    df = data.copy()
    prev_cols = df.columns.to_list()
    df[["chr", "start", "end"]] = df[element_column].str.split("-", expand = True)
    if chr_prefix == True:
        if not df["chr"].iloc[0].startswith("chr"):
            df["chr"] = "chr" + df["chr"].astype(str)
    elif chr_prefix == False:
        if df["chr"].iloc[0].startswith("chr"):
            df["chr"] = df["chr"].str.removeprefix("chr")
    df["chr"] = df["chr"]
    df[["start", "end"]] = df[["start", "end"]].astype(int)
    if element_gene_links == True:
        df["chr_gene"] = df["chr"].astype(str) + "_" + df["gene"]
        df = df[["chr_gene", "start", "end"] + prev_cols]
    else:
        df = df[["chr", "start", "end"] + prev_cols]
    return df

def bedfile_to_bedtool(bedfile):
    bed_df = pd.read_csv(bedfile, sep = "\t", header = None, usecols = range(3))
    if not str(bed_df[0].iloc[0]).startswith('chr'):
        bed_df[0] = 'chr' + bed_df[0].astype(str)
    return bed_df

def df_to_bedtool(bed_df, merge = False):
    if merge == True:
        return pybedtools.BedTool.from_dataframe(bed_df).sort().merge()
    else:
        return pybedtools.BedTool.from_dataframe(bed_df).sort()


def plot_scores(labels, mean, se, p_values = None, hline_y = None, xlabel = "Dataset + cell type", ylabel = "Mean score"):
    fig, ax = plt.subplots(figsize = (5,4))
    plt.scatter(labels, mean, color = "#7294d4", edgecolor = "black", zorder = 3)
    plt.errorbar(labels, mean, yerr = se, ls = "", capsize = 4, color = "black", zorder = 2)
    if hline_y != None:
        plt.hlines(xmin = ax.get_xlim()[0], xmax = ax.get_xlim()[1], y = hline_y, ls = "--", color = "black")
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    ax.set_xticklabels(labels, rotation = 45, ha = "right")
    if p_values is not None:
        sig_stars = [sig_star(p) for p in p_values]
        for i in range(len(sig_stars)):
            plt.text(labels[i], mean[i] + se[i], sig_stars[i], color = "black", ha = "center")
    plt.tight_layout()
            
def p_from_est_and_se(est, se, h0_mean = 0):
    z = (est - h0_mean) / se
    p = 2 * scipy.stats.norm.sf(abs(z))
    return p


def intersect_bed_dfs(bed_df1, bed_df2, loj = False):
    columns_df1, columns_df2 = bed_df1.columns.tolist(), bed_df2.columns.tolist()
    ncol_df1, ncol_df2 = len(columns_df1), len(columns_df2)
    intersected_df = pybedtools.BedTool.from_dataframe(bed_df1).intersect(
        pybedtools.BedTool.from_dataframe(bed_df2), 
        wa = True, wb = True, loj = loj
    ).to_dataframe()
    if len(intersected_df) == 0:
        return pd.DataFrame(columns = [columns_df1 + columns_df2])
    else:
        intersected_df.columns = columns_df1 + columns_df2
        columns_idx = [*list(range(3, ncol_df1)), *list(range(ncol_df1 + 3, ncol_df1 + ncol_df2))]
        intersected_df = intersected_df.drop_duplicates()
        return intersected_df.iloc[:, columns_idx]
        return intersected_df

def element_distance(df, element1_col, element2_col):
    tmp = df.copy()
    tmp = add_bed_columns(tmp, element_column = element1_col)
    tmp["loc1"] = tmp[["start", "end"]].mean(axis = 1)
    tmp = tmp.drop(columns = ["chr", "start", "end"])
    tmp = add_bed_columns(tmp, element_column = element2_col)
    tmp["loc2"] = tmp[["start", "end"]].mean(axis = 1)
    tmp = tmp.drop(columns = ["chr", "start", "end"])
    tmp["distance"] = (tmp["loc2"] - tmp["loc1"]).abs()
    return tmp["distance"]

def add_peak_column(df, chrom_col = 'chr', start_col = 'start', end_col = 'end'):
    df['peak'] = df[chrom_col] + '-' + df[start_col].astype(str) + '-' + df[end_col].astype(str)
    return df

def element_gene_distance(df, element_col = "peak", gene_col = "gene", tss_file = "/n/groups/price/elizabeth/data/gene_coords/STRANDED_TSS/tss.txt"):
    tmp = df.copy()
    tss = pd.read_csv(tss_file, sep = "\t")
    tss.columns = ['gene', 'gene_chr', 'coord']
    tss['gene_chr'] = ('chr' + tss['gene_chr'].astype(str))
    tss['coord'] = tss['coord'].astype(int)
    tmp = tmp.merge(tss, left_on = gene_col, right_on = "gene", how = "left")
    tmp = add_bed_columns(tmp, element_column = element_col)
    tmp["loc"] = tmp[["start", "end"]].mean(axis = 1)
    tmp["distance"] = (tmp["loc"] - tmp["coord"]).abs()
    tmp.loc[tmp['chr'] != tmp['gene_chr'], 'distance'] = np.nan
    tmp['distance'] = tmp['distance'].astype(float)
    tmp = tmp.drop(columns = ["chr", "gene_chr", "start", "end"])
    return tmp["distance"]

def gene_gene_distance(df, gene1_col = "gene1", gene2_col = "gene2", tss_file = "/n/groups/price/elizabeth/data/gene_coords/STRANDED_TSS/tss.txt"):
    tmp = df.copy()
    tss = pd.read_csv(tss_file, sep = "\t")
    tss.columns = ['gene', 'chr', 'coord']
    tss['coord'] = tss['coord'].astype(int)
    tmp = tmp.merge(tss[['gene', 'coord']], left_on = gene1_col, right_on = "gene", how = "left"
                   ).rename(columns = {'coord': 'gene1_coord'}).drop(columns = 'gene')
    tmp = tmp.merge(tss[['gene', 'coord']], left_on = gene2_col, right_on = "gene", how = "left"
                   ).rename(columns = {'coord': 'gene2_coord'}).drop(columns = 'gene')
    tmp["distance"] = (tmp["gene1_coord"] - tmp["gene2_coord"]).abs()
    return tmp['distance']

def get_gene_universe_file():
    return "/n/groups/price/elizabeth/data/gene_coords/STRANDED_TSS/tss.txt"

def get_gene_universe():
    gene_universe = pd.read_csv(get_gene_universe_file(), sep = '\t')
    gene_universe['chr'] = gene_universe['chr'].astype(str)
    return gene_universe

def read_coaccessibility(dataset, cell_type, method = "archr", folder_suffix = "", subfolder = ""):
    coaccess_file = get_coaccess_file(dataset, cell_type, method = method, folder_suffix = folder_suffix, subfolder = subfolder)
    print(coaccess_file)
    coaccess = read_correlations(coaccess_file)
    return coaccess

def read_coexpression(dataset, cell_type, method = 'archr', file_suffix = "_dist_1e+06_k100_knniter500_overlapcutoff0.8", subfolder = ""):
    coexpression_file = get_coexpression_file(dataset, cell_type, subfolder = subfolder, method = method, file_suffix = file_suffix)
    print(coexpression_file)
    coexpression = read_correlations(coexpression_file)
    coexpression = coexpression.rename(columns = {'gene': 'gene1', 'peak': 'gene2'})
    return coexpression

def process_coexpression(coexpression):
    coexpression = restrict_to_gene_universe(coexpression, df_gene_col = 'gene1')
    coexpression = restrict_to_gene_universe(coexpression, df_gene_col = 'gene2')
    coexpression = coexpression.drop_duplicates()
    return coexpression

def get_gene_chrom(df, gene_col = 'gene'):
    gene_universe = get_gene_universe()
    gene_universe = gene_universe.rename(columns = {'gene': gene_col})
    return df.merge(gene_universe[[gene_col, 'chr']])

def process_coaccessibility(coaccess, dataset, cell_type, folder_suffix = "", path_to_main_folder = ""):
    blacklist_peaks_file = get_blacklist_peaks_file(dataset, cell_type, path_to_main_folder = path_to_main_folder)
    blacklist_peaks = get_blacklist_peaks(blacklist_peaks_file)
    coaccess = remove_blacklist_peaks(coaccess, blacklist_peaks, col = "peak1")
    coaccess = remove_blacklist_peaks(coaccess, blacklist_peaks, col = "peak2")
    coaccess = remove_sex_chromosome_peaks(coaccess, col = "peak1")
    coaccess = remove_sex_chromosome_peaks(coaccess, col = "peak2")
    return coaccess

def read_archr_peak_gene_links(dataset, cell_type, subfolder = "", folder_suffix = "", file_suffix = "_dist_1e+06"):
    peak_gene_links_file = get_coactivity_file(dataset, cell_type, subfolder = subfolder, folder_suffix = folder_suffix, file_suffix = file_suffix)
    print(peak_gene_links_file)
    peak_gene_links = read_correlations(peak_gene_links_file)
    return peak_gene_links

def read_signac_peak_gene_links(dataset, cell_type, subfolder = ""):
    peak_gene_links_file = "%ssignac_%s_%s/pgl/signac_peak_gene_links.tsv" % (subfolder, dataset, cell_type)
    print(peak_gene_links_file)
    peak_gene_links = read_correlations(peak_gene_links_file)
    peak_gene_links = peak_gene_links.rename(columns = {'Score': 'Correlation'})
    return peak_gene_links


def process_peak_gene_links(peak_gene_links, dataset, cell_type, path_to_main_folder = ""):
    blacklist_peaks_file = get_blacklist_peaks_file(dataset, cell_type, path_to_main_folder = path_to_main_folder)
    blacklist_peaks = get_blacklist_peaks(blacklist_peaks_file)
    peak_gene_links = remove_blacklist_peaks(peak_gene_links, blacklist_peaks)
    peak_gene_links = remove_sex_chromosome_peaks(peak_gene_links)
    peak_gene_links = restrict_to_gene_universe(peak_gene_links)
    return peak_gene_links


def crispr_with_link_bed_columns(bedfile = "EPCrisprBenchmark_ensemble_data_GRCh38_table.tsv"):
    crispr = pd.read_csv(bedfile, sep = "\t",
                         names = ["chr", "start", "end", "gene", "crispr"])
    crispr["chr_gene"] = crispr["chr"] + "_" + crispr["gene"]
    crispr = crispr[["chr_gene", "start", "end", "crispr"]]
    crispr["crispr"] = 1 * crispr["crispr"]
    return crispr

def gwas_with_link_bed_columns(bedfile = "gwas_evaluation.tsv"):
    gwas = pd.read_csv(bedfile, sep = "\t")
    gwas = add_bed_columns(gwas, 'snp', element_gene_links = True)
    return gwas.drop(columns = ['snp', 'gene'])

def eqtl_with_link_bed_columns(bedfile = "/n/groups/price/elizabeth/new_pgboost_revisions/pgboost_training_sets/max_pip_49_susie_eqtl_pip-gt0.5-lt0.01.tsv.gz"):
    eqtl = pd.read_csv(bedfile, sep = "\t")
    eqtl = add_bed_columns(eqtl, 'snp', element_gene_links = True)
    return eqtl.drop(columns = ['snp', 'gene', 'pip'])


def merge_peak_gene_links_with_gold(peak_gene_links, gold):
    # Read in ArchR results and add peak-gene link bed columns
    peak_gene_links = add_bed_columns(peak_gene_links, element_gene_links = True)
    # Merge with gold (take link as positive if any overlapping gold element is a positive)
    merged = intersect_bed_dfs(peak_gene_links, gold)
    merged = merged.groupby(["peak", "gene"]).agg(max).reset_index()
    return merged

def merge_archr_peak_gene_links_with_crispr(dataset, cell_type, crispr, method = "archr"):
    # Read in ArchR results and add peak-gene link bed columns
    if method == "signac":
        peak_gene_links = read_signac_peak_gene_links(dataset, cell_type)
    elif method == "archr":
        peak_gene_links = read_archr_peak_gene_links(dataset, cell_type)
    peak_gene_links = process_peak_gene_links(peak_gene_links, dataset, cell_type)
    peak_gene_links = add_bed_columns(peak_gene_links, element_gene_links = True)
    # Merge with CRISPR (take link as positive if any overlapping CRISPR element is a validated enhancer)
    merged = intersect_bed_dfs(peak_gene_links, crispr)
    merged = merged.groupby(["peak", "gene"]).agg(max).reset_index()
    # Clean up and return
    merged = merged.rename(columns = {"Correlation": "peak_gene_correlation", "Score": "peak_gene_correlation"})
    if method == "archr":
        merged = merged.drop(columns = ["FDR"])
    return merged

def restrict_to_gene_universe(df, df_gene_col = 'gene', gene_universe_file = "/n/groups/price/elizabeth/data/gene_coords/STRANDED_TSS/tss.txt"):
    gene_universe = pd.read_csv(gene_universe_file, sep = "\t")
    df = df[df[df_gene_col].isin(gene_universe["gene"])]
    return df.reset_index(drop = True)


def bootstrap_p_diff_mean_from_boot_estimates(est1, est2, boot1, boot2):
    diff_est = est1 - est2
    diff_boot = boot1 - boot2
    diff_bootstrap_se = np.std(diff_boot)
    z = diff_est / diff_bootstrap_se
    p_diff = 2 * scipy.stats.norm.sf(abs(z))
    return p_diff, diff_est, diff_bootstrap_se


def square_column(df, col):
    df["%s_squared" % col] = df[col] **2
    return df

def get_peak_length(peaks_series):
    coords = peaks_series.str.split("-", expand = True)
    length = coords[2].astype(int) - coords[1].astype(int)
    return length


def get_weighted_jackknife_estimate_and_se(theta, held_out_estimates, m, n):
    theta_jack = sum(theta - held_out_estimates) + sum((m * held_out_estimates) / n)
    h = n / m
    tau = (h * theta) - ((h - 1) * held_out_estimates)
    var_jack = np.mean((tau - theta_jack)**2 / (h - 1))
    se_jack = np.sqrt(var_jack)
    return theta_jack, se_jack

def block_jackknife(df, block_column, function, weighted = True, held_out_estimates_only = False, num_estimands = 1, *args, **kwargs):
    # Define total number of data points
    n = len(df)
    
    # Define blocks and g (number of blocks)
    blocks = np.sort(df[block_column].unique())
    g = len(blocks)
    
    # Get overall estimate
    thetas = function(df, **kwargs)
    
    # Get estimates and number of data points while holding out each block
    held_out_estimates = np.zeros((g, num_estimands))
    m = np.zeros(g)
    for i, block in enumerate(blocks):
        tmp = df[df[block_column] != block].reset_index(drop = True)
        if weighted == True:
            m[i] = n - len(tmp)
        else:
            m[i] = n / g
        held_out_estimates[i] = function(tmp, **kwargs)
    if held_out_estimates_only == True:
        return thetas, held_out_estimates, m
    
    ses = np.zeros(num_estimands)
    for i in range(num_estimands):
        if weighted == True:
            # Get jackknife estimate and SE
            theta_jack, se_jack = get_weighted_jackknife_estimate_and_se(thetas[i], held_out_estimates[:,i], m, n)
        else:
            # Get jackknife estimate
            se_jack = np.std(held_out_estimates[:,1]) * np.sqrt(g - 1)
        ses[i] = se_jack
    # Return jackknife estimate and SE
    if num_estimands == 1:
        return thetas[0], ses[0], held_out_estimates, m
    else:
        return thetas, ses, held_out_estimates, m


def read_scores_chromwise(dataset, cell_type, label, score_type = 'coaccessibility'):
    scores_files = glob("archr_%s_%s/%s_scores%s/*" % (dataset, cell_type, score_type, label))
    print(scores_files[0])
    scores = pd.read_csv(scores_files[0], sep = "\t")
    for file in scores_files[1:]:
        tmp = pd.read_csv(file, sep = "\t")
        scores = pd.concat([scores, tmp])
    scores = scores.rename(columns = {"uncorrected": "%s%s" % (score_type, label),
                                      "corrected": "corrected_%s%s" % (score_type, label)})
    scores = scores.drop_duplicates()
    return scores

def column_average(df, avg_col, weight_col = None):
    if weight_col == None:
        return [df[avg_col].mean()]
    else:
        weights = df[weight_col] / sum(df[weight_col])
        return [(df[avg_col] * weights).mean()]
    
def column_average_multi(df, avg_cols, weight_col = None):
    if weight_col == None:
        return df[avg_cols].mean().to_list()
    else:
        weights = df[weight_col] / sum(df[weight_col])
        return (df[avg_cols] * weights).mean().to_list()

def read_coexpression_scores(dataset, cell_type):
    coexpression_scores_file = "archr_%s_%s/coexpression_scores_all_genes.tsv" % (dataset, cell_type)
    print(coexpression_scores_file)
    coexpression_scores = read_scores(coexpression_scores_file)
    return coexpression_scores

def read_gene_coactivity_scores(dataset, cell_type):
    scores_file = "archr_%s_%s/coactivity_scores_by_gene_all_peaks.tsv" % (dataset, cell_type)
    print(scores_file)
    scores = pd.read_csv(scores_file, sep = '\t')
    scores = scores.rename(columns = dict(zip(
        scores.columns, 
        ["gene"] + ["%s_gene_coactivity" % c for c in scores.columns if c != "gene"])))
    return scores

def drop_dups_two_columns(in_df, col1, col2):
    df = in_df.copy()
    df['col_pair_sorted'] = df.apply(lambda row: tuple(sorted([row[col1], row[col2]])), axis=1)
    df = df.drop_duplicates('col_pair_sorted')
    df = df.drop(columns = ['col_pair_sorted'])
    return df

# Related to S-LDSC analyses
def get_annotation_sd(df, annot_names = None):
    return df[annot_names].std()

def get_annotation_size(df, annot_names = None, proportion = False):
    if proportion == True:
        return (df[annot_names].mean())
    else:
        return (df[annot_names].sum())

def get_coefficients(df, predictors, target_predictor = None, covariates = [], outcome_variable = "coactivity_score", weights_col = None):
    Y = df[outcome_variable].to_numpy()
    X = df[predictors + covariates]
    X = sm.add_constant(X).to_numpy()
    if weights_col != None:
        W = sparse.diags(df[weights_col], format='dia')
        betas = np.linalg.inv(X.T @ W @ X) @ X.T @ W @ Y
    else:
        betas = np.linalg.inv(X.T @ X) @ X.T @ Y
    if target_predictor != None:
        return betas[1:][predictors.index(target_predictor)]
    elif len(covariates) > 0:
        return betas[1:(-1 * len(covariates))]
    else:
        return betas[1:]
    
def get_intercept(df, predictors, outcome_variable = "coactivity_score"):
    
    formula = "%s ~ %s" % (outcome_variable, " + ".join(predictors))
    model = sm.OLS.from_formula(formula, data = df).fit()
    return np.array([model.params[0]])

def get_annotation_h2(df, annot_name = None):
    return (df[annot_name] * df['per_peak_h2']).sum()


def get_annotation_per_peak_h2(df, annot_names = None, taus = None):
    for annot_name, tau in zip(annot_names, taus):
        df['%s_per_peak_h2' % annot_name] = df[annot_name] * tau
    h2_cols = ['%s_per_peak_h2' % a for a in annot_names]
    df['per_peak_h2'] = df[h2_cols].sum(axis = 1)
    return df

def get_annotation_h2_multiple(in_df, annot_names = None, taus = None):
    df = in_df.copy()
    df = get_annotation_per_peak_h2(df, annot_names = annot_names, taus = taus)
    h2s = np.zeros(len(annot_names))
    for i in range(len(annot_names)):
        h2s[i] = get_annotation_h2(df, annot_names[i])
    return h2s

def r2_predictor_configurations(data, outcome_variable, predictor_sets):
    r2_list = []
    for predictor_set in predictor_sets:
        formula = '%s ~ %s' % (outcome_variable, ' + '.join(predictor_set))
        model = sm.OLS.from_formula(formula, data = data).fit()
        r2_list.append(model.rsquared)
    return r2_list

def meta_analysis(effect_sizes, standard_errors):
    result = combine_effects(effect_sizes, standard_errors ** 2).summary_frame()
    meta_est = result.loc["random effect", "eff"]
    meta_se = result.loc["random effect", "sd_eff"]
    return meta_est, meta_se
          
def meta_analysis_multi(effect_sizes, standard_errors):
    # Index 0 should index meta-analysis units (e.g. datasets)
    # Index 1 should index quantities to meta-analyze
    num_quantities = effect_sizes.shape[1]
    meta_est = np.zeros(num_quantities)
    meta_se = np.zeros(num_quantities)
    for i in range(num_quantities):
        meta_est[i], meta_se[i] = meta_analysis(effect_sizes[:,i], standard_errors[:,i])
    return meta_est, meta_se

def make_dir(path_to_dir):
    if not os.path.exists(path_to_dir):
        os.makedirs(path_to_dir, exist_ok = True)
        
def get_background_corr_sq(dataset, cell_type, metric, path_to_main_folder = ".", file_suffix = ""):
    file = '%s/archr_%s_%s/background_%s_via_downsampling_metacells/k100_knniter500_overlapcutoff0.8%s.txt' % (path_to_main_folder, dataset, cell_type, metric, file_suffix)
    background = pd.read_csv(file, sep = '\t')['noise'].iloc[0]
    return background

def get_rna_matrix(dataset, cell_type, suffix = "", path_to_main_folder = "", remove_nonexpressed_genes = True):
    rna_matrix_file = "%sarchr_%s_%s/RNA_metacell_matrices/RNA_metacell_matrix_k100_knniter500_overlapcutoff0.8%s.rds" % (path_to_main_folder, dataset, cell_type, suffix)
    rna_matrix = pyreadr.read_r(rna_matrix_file)[None].T
    gene_coords = get_gene_universe()
    gene_coords['chr'] = gene_coords['chr'].astype(str)
    gene_univ = set(gene_coords['gene'])
    rna_matrix = rna_matrix[[gene for gene in rna_matrix.columns if gene in gene_univ]]
    if remove_nonexpressed_genes == True:
        rna_matrix = rna_matrix[rna_matrix.columns[rna_matrix.sum() > 0]]
    return rna_matrix

def get_atac_matrix(dataset, cell_type, suffix = "", path_to_main_folder = "", remove_nonactive_peaks = True):
    atac_matrix_file = "%sarchr_%s_%s/atac_metacell_matrices/atac_metacell_matrix_k100_knniter500_overlapcutoff0.8%s.rds" % (path_to_main_folder, dataset, cell_type, suffix)
    atac_matrix = pyreadr.read_r(atac_matrix_file)[None].T
    blacklist_peaks = get_blacklist_peaks(get_blacklist_peaks_file(dataset, cell_type))
    atac_matrix = atac_matrix[[p for p in atac_matrix.columns if p not in blacklist_peaks]]
    atac_matrix = atac_matrix[[p for p in atac_matrix.columns if not p.startswith('chrX') and not p.startswith('chrY') and not p.startswith('chrM')]]
    if remove_nonactive_peaks == True:
        atac_matrix = atac_matrix[atac_matrix.columns[atac_matrix.sum() > 0]]
    return atac_matrix



def standardize_df_columns_to_matrix(input_df):
    colnames = input_df.columns.tolist()
    rownames = input_df.index.tolist()
    matrix_values = input_df.values
    matrix_values_std = (matrix_values - matrix_values.mean(axis = 0)) / (matrix_values.std(axis = 0))
    return matrix_values_std, rownames, colnames

def array_sumstats(arr):
    print("Mean: %s" % np.mean(arr))
    print("Median: %s" % np.median(arr))
    print("Minimum: %s" % np.min(arr))
    print("Maximum: %s" % np.max(arr))
    print("1st, 2nd, 3rd quartiles: %s, %s, %s" %
          (np.round(np.quantile(arr, 0.25), 5),
           np.round(np.quantile(arr, 0.5), 5),
           np.round(np.quantile(arr, 0.75), 5))
         )

def rsq_linalg(X, y):
    beta = np.linalg.inv(X.T @ X) @ X.T @ y
    preds = X @ beta
    rss = np.sum((y - preds)**2)
    y_mean = np.mean(y)
    tss = np.sum((y - y_mean)**2)
    r_squared = 1 - (rss / tss)
    return r_squared


def plot_corr(mtx):
    plt.figure(figsize = (6,6))
    heatmap = sns.heatmap(mtx,
                          square = True,
                          linewidths = .5,
                          cmap = 'Purples',
                          cbar_kws = {'shrink': .4,
                                    'ticks' : [-1, -.5, 0, 0.5, 1]},
                                     #'label': r"$r$"},
                          vmin = 0,
                          vmax = 1,
                          annot = True,
                          annot_kws = {"size": 10},
                         fmt='.2f')
    plt.xticks(rotation = 0, fontsize = 10)
    plt.yticks(rotation = 90, ha = "right", fontsize = 10)
    #plt.text(x = 4.25, y = 0.8, s = r"$r$", fontsize = 11)
    
def merge_gold_and_candidate_links(gold_standard_file, candidate_links, gold_merge_how = "left", replace_nonscored = "0"):
    # Prepare gold standard in bed format
    gold_standard = pd.read_csv(gold_standard_file, sep = '\t')
    gold_standard['chr_gene'] = gold_standard['chr'] + '_' + gold_standard['gene']
    gold_standard = gold_standard[['chr_gene', 'start', 'end', 'gold']] 
    
    # Prepare candidate links in bed format
    candidate_links = add_bed_columns(candidate_links, element_gene_links = True)
    candidate_links = candidate_links[['chr_gene', 'start', 'end', 'peak', 'gene']]
    
    # Merge (and take maximum gold value per element)
    if gold_merge_how == "inner":
        gold_X_candidate = intersect_bed_dfs(candidate_links, gold_standard)
    elif gold_merge_how == "left":
        gold_X_candidate = intersect_bed_dfs(candidate_links, gold_standard, loj = True)
        gold_X_candidate['gold'] = gold_X_candidate['gold'].str.replace('.', replace_nonscored, regex = True)
    
    gold_X_candidate['gold'] = gold_X_candidate['gold'].astype(int)
    
    # Take maximum value of 'gold' per element-gene link
    gold_X_candidate = gold_X_candidate.groupby(['peak', 'gene']).agg(max).reset_index()
    
    return gold_X_candidate


def restrict_analysis_to_positive_negative_genes(df, restrict_col = 'gene'):
    print("%s links pre-subset" % len(df))
    positive_genes = df[df["gold"] == 1][restrict_col].unique()
    negative_genes = df[df["gold"] == 0][restrict_col].unique()
    positive_negative_genes = set.intersection(set(positive_genes), set(negative_genes))
    print("Subsetting to %s positive-negative %ss from %s positive %ss and %s negative %ss" % (len(positive_negative_genes), restrict_col, len(positive_genes), restrict_col, len(negative_genes), restrict_col))
    df = df[df[restrict_col].isin(positive_negative_genes)]
    print("%s links post-subset" % len(df))
    return df


def read_all_peaks(dataset, cell_type, path_to_main_folder = "./"):
    print("Reading!")
    peaks_bedfile = "%sarchr_%s_%s/all_called_peaks_minus_promoters_and_sex_chromosomes.bed" % (path_to_main_folder, dataset, cell_type)
    print(peaks_bedfile)
    all_peaks = pd.read_csv(peaks_bedfile, sep = '\t', header = None)
    all_peaks = all_peaks.rename(columns = {3: 'peak'})[['peak']]
    all_peaks['peak'] = 'chr' + all_peaks['peak']
    return all_peaks

def adjust_start_coords(df, element_col = 'peak', adjustment = 1):
    df = add_bed_columns(df, element_column = element_col)
    df[element_col] = df['chr'].astype(str) + '-' + (df['start'].astype(int) + adjustment).astype(str) + '-' + df['end'].astype(str)
    df = df.drop(columns = ['chr', 'start', 'end'])
    return df

def load_table_dataset_labels():
    datasets = ['Xu K562', 'Satpathy K562', 'SHARE-seq LCL', 'Luecken BMMC', 'Luecken BMMC', 'Luecken BMMC', 'Luecken BMMC']
    cell_types = ['K562', 'lymphoblastoid cell line (LCL)', 'T', 'B', 'myeloid', 'erythroid', 'K562']
    table = pd.DataFrame({'Dataset': datasets, 'Cell type': cell_types})
    return table

def save_supp_table(df, label, path_to_main_folder = "./"):
    df.to_csv('%stables/%s.tsv' % (path_to_main_folder, label), sep = '\t', index = False)
    
def load_dataset_label_dict():
    label_dict = {
        'Xu K562': 'Xu K562',
        'McGinnisDataset10 K562': 'Satpathy K562',
        'shareseq lcl': 'SHARE-seq LCL',
        'neurips T': 'Luecken T',
        'neurips B': 'Luecken B',
        'neurips mono': 'Luecken myeloid',
        'neurips eryth': 'Luecken erythroid'
    }
    return label_dict

def make_chr_col(df, peak_column = "peak", delimiter = "-"):
    df['chr'] = df[peak_column].str.split(delimiter).str[0]
    return df