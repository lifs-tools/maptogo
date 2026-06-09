# Systematic regulation score

# Dear Bo, dear SLING group,
#
# Thank you very much for the mail. I am very delighted that you consider me for hosting the data analysis workshop. With please, I accept your generous invitation :-) For the moment, my calendar looks very empty for April next year. Therefore, I am fine with both weeks. Once decided, please let us know. Looking forward to see you guys next year in person, again.
#
# Cheers,
# Dominik


from dataframe import *
import numpy as np
from time import time
import random
import numpy as np
from scipy import stats
from typing import Dict, List, Tuple
from term_sets_bp import term_sets
from scipy.stats import norm

print(f"loaded {len(term_sets)} terms")






from scipy import stats
import numpy as np
from statsmodels.stats.multitest import multipletests

def preprocess(df, min_cnt_per_grp=2):
    df_work = df.copy().astype(float)
    df_grp = df_work.groupby(level=(0,))
    no_var_features = (d := df_grp.var().min())[(d <= 0) | (np.isnan(d)) | (np.isinf(d))].index
    too_few_features = (d := df_grp.count().min())[d < min_cnt_per_grp].index
    all_removed = set(no_var_features) | set(too_few_features)
    df_clean = df_work.drop(columns=list(all_removed))
    df_clean[df_clean <= 0] = np.nan
    return df_clean, all_removed



def log2_dataset(df):
    df_work = df.copy()
    df_work.iloc[:, :] = np.where(~np.isnan(df_work.values), np.log2(df_work.values), np.nan)
    return df_work

import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

def pathway_network_enrichment(df, term_sets, conditions, min_cnt_per_grp=2):
    df_clean, _ = preprocess(df, min_cnt_per_grp)
    X_df = log2_dataset(df_clean).replace([np.inf, -np.inf], np.nan)
    X = X_df.values.astype(np.float64)
    feature_idx = {f: i for i, f in enumerate(X_df.columns)}

    conditions = np.asarray(conditions)
    groups = np.unique(conditions)

    # ------------------------------------------------------------
    # differential expression DE (feature-wise)
    # ------------------------------------------------------------
    group_means = np.array([np.nanmean(X[conditions == g], axis=0) for g in groups])
    overall_mean = np.nanmean(X, axis=0)
    ss_between = np.sum((group_means - overall_mean) ** 2, axis=0)
    ss_within = np.array([np.nanvar(X[conditions == g], axis=0) for g in groups]).sum(axis=0)
    de_score = ss_between / (ss_between + ss_within + 1e-12)

    # DE direction: linear regression slope of group means vs condition label
    x = groups.astype(float) - groups.astype(float).mean()
    feature_direction = (x @ group_means) / (x @ x)

    # ------------------------------------------------------------
    # correlation matrix
    # ------------------------------------------------------------
    X_corr = X_df.fillna(X_df.mean()).values.astype(np.float64)
    X_std = (X_corr - X_corr.mean(axis=0, keepdims=True)) / X_corr.std(axis=0, keepdims=True)
    corr_abs = np.abs((X_std.T @ X_std) / (X_std.shape[0] - 1))
    np.fill_diagonal(corr_abs, np.nan)
    node_conn_abs = np.nanmean(corr_abs, axis=1)

    # reward pairs of biomolecules that have similar shifts over conditions (node_conn_abs)
    # and that actually change over conditions (de_score)
    combined_score_abs = node_conn_abs * de_score

    # ------------------------------------------------------------
    # pathway mapping
    # ------------------------------------------------------------
    term_sets_idx = {
        pw: np.array([feature_idx[b] for b in biomols if b in feature_idx], dtype=np.int64)
        for pw, biomols in term_sets.items()
    }
    term_sets_idx = {k: v for k, v in term_sets_idx.items() if len(v) >= 5}

    # global node statistics (for z-score baseline)
    mu_abs = np.mean(combined_score_abs)
    sigma_abs = np.std(combined_score_abs)

    # ------------------------------------------------------------
    # state assignment
    # ------------------------------------------------------------
    def assign_state(pval, direction, threshold = 0.05, dir_threshold = 0.1):
        if pval >= threshold:
            return "unchanged"
        if np.abs(direction) < dir_threshold:
            return "modulated"
        return "enriched" if direction > 0 else "suppressed"

    # ------------------------------------------------------------
    # scoring
    # ------------------------------------------------------------
    results = []
    for pw, idx in term_sets_idx.items():
        k = len(idx)
        S_abs = corr_abs[np.ix_(idx, idx)]
        node_scores_abs = np.nanmean(S_abs, axis=1) * de_score[idx]
        if len(node_scores_abs) == 0:
            continue

        mean_conn_abs = np.mean(node_scores_abs)
        pathway_signal_abs = np.mean(combined_score_abs[idx])
        pathway_direction = np.nanmedian(feature_direction[idx])

        se = sigma_abs / np.sqrt(k)
        z = (mean_conn_abs - mu_abs) / se
        pval = float(np.exp(norm.logsf(z)))

        if np.isnan(pval):
            print(np.isnan(ss_within).sum())
            exit()

        results.append([
            pw,
            float(mean_conn_abs),
            float(pathway_signal_abs),
            float(z),
            float(pval),
            int(k),
            float(pathway_direction),
        ])

    # ------------------------------------------------------------
    # FDR correction + state assignment
    # ------------------------------------------------------------
    if len(results) == 0:
        return []

    pvals = np.array([r[4] for r in results])
    _, qvals, _, _ = multipletests(pvals, method="fdr_bh")

    for r, q in zip(results, qvals):
        state = assign_state(q, r[6])
        r.append(float(q))   # index 7: qval
        r.append(state)      # index 8: state

    results.sort(key=lambda x: x[7], reverse=True)
    return results















results = pathway_network_enrichment(data_frame, term_sets, [0, 0, 0, 1, 1, 1])
for r in results:
    print(r[0], r[7], r[6], r[8])



























exit()


import numpy as np

def pathway_network_enrichment(df, term_sets, conditions, min_cnt_per_grp=2):

    df_clean, _ = preprocess(df, min_cnt_per_grp)
    X_df = log2_dataset(df_clean).replace([np.inf, -np.inf], np.nan).fillna(0)
    X = X_df.values.astype(np.float64)

    feature_idx = {f: i for i, f in enumerate(X_df.columns)}

    term_sets_idx = {
        pw: sorted(feature_idx[b] for b in biomols if b in feature_idx)
        for pw, biomols in term_sets.items()
    }
    term_sets_idx = {k: np.array(v) for k, v in term_sets_idx.items() if len(v) >= 5}

    conditions = np.asarray(conditions)
    groups = np.unique(conditions)

    # ------------------------------------------------------------
    # differential expression (feature-wise)
    # ------------------------------------------------------------

    group_means = np.array([
        X[conditions == g].mean(axis=0)
        for g in groups
    ])

    overall_mean = X.mean(axis=0)

    ss_between = np.sum((group_means - overall_mean) ** 2, axis=0)

    ss_within = np.array([
        np.var(X[conditions == g], axis=0)
        for g in groups
    ]).sum(axis=0) + 1e-12

    de_score = ss_between / (ss_between + ss_within)  # per-feature signal

    # ------------------------------------------------------------
    # correlation matrix
    # ------------------------------------------------------------


    X = X_df.values.astype(np.float64)
    cv = X.std(axis = 0) / X.mean(axis = 0)
    n = len(cv)

    # standardize
    X = (X - X.mean(axis=0, keepdims=True)) / X.std(axis=0, keepdims=True)

    # fast correlation
    corr_abs = np.abs((X.T @ X) / (X.shape[0] - 1))
    np.fill_diagonal(corr_abs, np.nan)

    # node connectivity
    node_conn_abs = np.nanmean(corr_abs, axis=1)


    # combined gene score: expression × connectivity
    cv_de_scores = de_score #* cv / (cv.mean() + cv)
    combined_score_abs = node_conn_abs * cv_de_scores

    # ------------------------------------------------------------
    # pathway mapping
    # ------------------------------------------------------------

    term_sets_idx = {
        pw: np.array(
            [feature_idx[b] for b in biomols if b in feature_idx],
            dtype=np.int64
        )
        for pw, biomols in term_sets.items()
    }

    term_sets_idx = {
        k: v for k, v in term_sets_idx.items()
        if len(v) >= 5
    }

    # global node statistics (for z-score baseline)
    mu_abs = np.mean(combined_score_abs)
    sigma_abs = np.std(combined_score_abs)

    # ------------------------------------------------------------
    # scoring
    # ------------------------------------------------------------

    results = []

    for pw, idx in term_sets_idx.items():

        k = len(idx)
        S_abs = corr_abs[np.ix_(idx, idx)]

        node_scores_abs = np.nanmean(S_abs, axis=1)
        node_scores_abs = node_scores_abs * cv_de_scores[idx]

        if len(node_scores_abs) == 0:
            continue

        mean_conn_abs = np.mean(node_scores_abs)

        # pathway functional signal (DE × connectivity)
        pathway_signal_abs = np.mean(combined_score_abs[idx])

        # node-based z-score
        se = sigma_abs / np.sqrt(k)
        z = (mean_conn_abs - mu_abs) / se
        pval = np.exp(norm.logsf(z))

        # if pw.startswith("Acyl chain remodelling"):
        #     print(pw, mean_conn_abs, z, pval)

        results.append([
            pw,
            float(mean_conn_abs),
            float(pathway_signal_abs),
            float(z),
            float(pval),
            int(k),
        ])

    # ------------------------------------------------------------
    # FDR correction
    # ------------------------------------------------------------

    if len(results) == 0:
        return []

    pvals = np.array([r[4] for r in results])
    _, qvals, _, _ = multipletests(pvals, method="fdr_bh")

    for r, q in zip(results, qvals):
        r.append(float(q))

    results.sort(key=lambda x: x[6], reverse = True)

    return results
