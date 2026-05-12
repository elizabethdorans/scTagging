from pathlib import Path

import numpy as np
import pandas as pd
import pybedtools
import scipy
import statsmodels.api as sm
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from scipy import sparse


REPO_DIR = Path(__file__).resolve().parent
DEFAULT_GENE_TSS_FILE = str(REPO_DIR / "gene_TSS.txt")


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


def read_scores(file, rename_suffix=None, scored_unit="peak"):
    scores = pd.read_csv(file, sep="\t")
    if rename_suffix is not None:
        scores = scores.rename(columns=dict(zip(
            scores.columns,
            [scored_unit] + [
                "%s_%s" % (c, rename_suffix)
                for c in scores.columns
                if c != scored_unit
            ],
        )))
    return scores


def get_p_from_est_and_se(est, se, null=0, return_z=False):
    est = est - null
    z = est / se
    p = 2 * scipy.stats.norm.sf(abs(z))
    if return_z:
        return z
    return p


def plot_across_conditions_multi(
    point_estimates,
    standard_errors,
    labels,
    title="",
    xlabel="Condition",
    ylabel="Estimate",
    figsize=(3, 4),
    xtick_rotation=0,
    xtick_ha="center",
    pvals=None,
    null=0,
    v_space=0,
    box_index=None,
    axis=None,
):
    if axis is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        ax = axis
    ax.bar(labels, [e - null for e in point_estimates], bottom=null, color="#7294d4")
    if np.max(standard_errors) > 0:
        ax.errorbar(
            labels,
            point_estimates,
            ls="",
            xerr=None,
            yerr=standard_errors,
            color="black",
            capsize=4,
        )
    ax.set_xticks(labels, labels, rotation=xtick_rotation, ha=xtick_ha)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if pvals is not None:
        for i in range(len(labels)):
            if point_estimates[i] > null:
                ax.text(
                    labels[i],
                    point_estimates[i] + standard_errors[i],
                    sig_star(pvals[i]),
                    color="green",
                    ha="center",
                    va="bottom",
                )
            else:
                ax.text(
                    labels[i],
                    point_estimates[i] - standard_errors[i],
                    sig_star(pvals[i]),
                    color="red",
                    ha="center",
                    va="top",
                )
    ax.set_xlim(-0.75, len(labels) - 0.5)
    ax.set_ylim(
        min(
            ax.get_ylim()[0],
            min([e - s for e, s in zip(point_estimates, standard_errors)]) - v_space,
        ),
        max(
            ax.get_ylim()[1],
            max([e + s for e, s in zip(point_estimates, standard_errors)]),
        ) + v_space,
    )
    if box_index is not None:
        ax.add_patch(Rectangle(
            (box_index - 0.5, 0),
            1.0,
            1.025 * (point_estimates[box_index] + standard_errors[box_index]),
            fill=False,
            edgecolor="black",
            lw=3,
        ))
    if axis is None:
        return fig, ax


def bedfile_to_bedtool(bedfile):
    bed_df = pd.read_csv(bedfile, sep="\t", header=None, usecols=range(3))
    if not str(bed_df[0].iloc[0]).startswith("chr"):
        bed_df[0] = "chr" + bed_df[0].astype(str)
    return bed_df


def df_to_bedtool(bed_df, merge=False):
    bed = pybedtools.BedTool.from_dataframe(bed_df).sort()
    if merge:
        return bed.merge()
    return bed


def intersect_bed_dfs(bed_df1, bed_df2, loj=False):
    columns_df1, columns_df2 = bed_df1.columns.tolist(), bed_df2.columns.tolist()
    ncol_df1, ncol_df2 = len(columns_df1), len(columns_df2)
    intersected_df = pybedtools.BedTool.from_dataframe(bed_df1).intersect(
        pybedtools.BedTool.from_dataframe(bed_df2),
        wa=True,
        wb=True,
        loj=loj,
    ).to_dataframe()
    if len(intersected_df) == 0:
        return pd.DataFrame(columns=columns_df1 + columns_df2)
    intersected_df.columns = columns_df1 + columns_df2
    columns_idx = [
        *list(range(3, ncol_df1)),
        *list(range(ncol_df1 + 3, ncol_df1 + ncol_df2)),
    ]
    intersected_df = intersected_df.drop_duplicates()
    return intersected_df.iloc[:, columns_idx]


def add_peak_column(df, chrom_col="chr", start_col="start", end_col="end"):
    df["peak"] = (
        df[chrom_col] + "-" + df[start_col].astype(str) + "-" + df[end_col].astype(str)
    )
    return df


def get_gene_universe_file():
    return DEFAULT_GENE_TSS_FILE


def get_gene_universe():
    gene_universe = pd.read_csv(get_gene_universe_file(), sep="\t")
    gene_universe["chr"] = gene_universe["chr"].astype(str)
    return gene_universe


def get_gene_chrom(df, gene_col="gene"):
    gene_universe = get_gene_universe()
    gene_universe = gene_universe.rename(columns={"gene": gene_col})
    return df.merge(gene_universe[[gene_col, "chr"]])


def square_column(df, col):
    df["%s_squared" % col] = df[col] ** 2
    return df


def get_weighted_jackknife_estimate_and_se(theta, held_out_estimates, m, n):
    theta_jack = sum(theta - held_out_estimates) + sum((m * held_out_estimates) / n)
    h = n / m
    tau = (h * theta) - ((h - 1) * held_out_estimates)
    var_jack = np.mean((tau - theta_jack) ** 2 / (h - 1))
    se_jack = np.sqrt(var_jack)
    return theta_jack, se_jack


def block_jackknife(
    df,
    block_column,
    function,
    weighted=True,
    held_out_estimates_only=False,
    num_estimands=1,
    *args,
    **kwargs,
):
    n = len(df)
    blocks = np.sort(df[block_column].unique())
    g = len(blocks)
    thetas = function(df, **kwargs)

    held_out_estimates = np.zeros((g, num_estimands))
    m = np.zeros(g)
    for i, block in enumerate(blocks):
        tmp = df[df[block_column] != block].reset_index(drop=True)
        if weighted:
            m[i] = n - len(tmp)
        else:
            m[i] = n / g
        held_out_estimates[i] = function(tmp, **kwargs)
    if held_out_estimates_only:
        return thetas, held_out_estimates, m

    ses = np.zeros(num_estimands)
    for i in range(num_estimands):
        if weighted:
            _, se_jack = get_weighted_jackknife_estimate_and_se(
                thetas[i],
                held_out_estimates[:, i],
                m,
                n,
            )
        else:
            se_jack = np.std(held_out_estimates[:, i]) * np.sqrt(g - 1)
        ses[i] = se_jack
    if num_estimands == 1:
        return thetas[0], ses[0], held_out_estimates, m
    return thetas, ses, held_out_estimates, m


def get_annotation_sd(df, annot_names=None):
    return df[annot_names].std()


def get_annotation_size(df, annot_names=None, proportion=False):
    if proportion:
        return df[annot_names].mean()
    return df[annot_names].sum()


def get_coefficients(
    df,
    predictors,
    target_predictor=None,
    covariates=[],
    outcome_variable="coactivity_score",
    weights_col=None,
):
    Y = df[outcome_variable].to_numpy()
    X = df[predictors + covariates]
    X = sm.add_constant(X).to_numpy()
    if weights_col is not None:
        W = sparse.diags(df[weights_col], format="dia")
        betas = np.linalg.inv(X.T @ W @ X) @ X.T @ W @ Y
    else:
        betas = np.linalg.inv(X.T @ X) @ X.T @ Y
    if target_predictor is not None:
        return betas[1:][predictors.index(target_predictor)]
    if len(covariates) > 0:
        return betas[1:(-1 * len(covariates))]
    return betas[1:]


def get_annotation_h2(df, annot_name=None):
    return (df[annot_name] * df["per_peak_h2"]).sum()


def get_annotation_per_peak_h2(df, annot_names=None, taus=None):
    for annot_name, tau in zip(annot_names, taus):
        df["%s_per_peak_h2" % annot_name] = df[annot_name] * tau
    h2_cols = ["%s_per_peak_h2" % a for a in annot_names]
    df["per_peak_h2"] = df[h2_cols].sum(axis=1)
    return df


def get_annotation_h2_multiple(in_df, annot_names=None, taus=None):
    df = in_df.copy()
    df = get_annotation_per_peak_h2(df, annot_names=annot_names, taus=taus)
    h2s = np.zeros(len(annot_names))
    for i in range(len(annot_names)):
        h2s[i] = get_annotation_h2(df, annot_names[i])
    return h2s
