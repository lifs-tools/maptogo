# Systematic regulation score

from dataframe import *
import numpy as np
from time import time
from term_sets_pw import term_sets
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from collections import Counter

print(f"loaded {len(term_sets)} terms")



def preprocess(df, min_cnt_per_grp = 2):
    df_work = df.copy().astype(float)
    df_grp = df_work.groupby(level = (0, ))
    no_var_features = (d := df_grp.var().min())[(d <= 0) | (np.isnan(d)) | (np.isinf(d))].index
    too_few_features = (d := df_grp.count().min())[d < min_cnt_per_grp].index
    all_removed = set(no_var_features) | set(too_few_features)
    df_clean = df_work.drop(columns = list(all_removed))
    df_clean[df_clean <= 0] = np.nan
    df_work = df_clean.replace([np.inf, -np.inf], np.nan)
    df_work.iloc[:, :] = np.where(~np.isnan(df_work.values), np.log2(df_work.values), np.nan)

    return df_work, all_removed




def pathway_network_enrichment(X_df, term_sets, conditions, min_cnt_per_grp = 2, min_cnt_per_term = 5):

    for i, c in Counter(conditions).items():
        if c < min_cnt_per_grp:
            raise Exception(f"Group {i} has less than {min_cnt_per_grp} samples / observations.")

    X = X_df.values.astype(np.float64)
    feature_idx = {f: i for i, f in enumerate(X_df.columns)}

    conditions = np.asarray(conditions)
    groups = np.unique(conditions)

    # ------------------------------------------------------------
    # differential expression DE (feature-wise)
    # ------------------------------------------------------------
    group_means = np.array([np.nanmean(X[conditions == g], axis = 0) for g in groups])
    overall_mean = np.nanmean(X, axis = 0)
    ss_between = np.sum((group_means - overall_mean) ** 2, axis = 0)
    ss_within = np.array([np.nanvar(X[conditions == g], axis = 0) for g in groups]).sum(axis = 0)
    de_score = ss_between / (ss_between + ss_within)


    # DE direction: linear regression slope of group means vs condition label
    x = groups.astype(float) - groups.astype(float).mean()
    feature_direction = (x @ group_means) / (x @ x)

    # ------------------------------------------------------------
    # correlation matrix
    # ------------------------------------------------------------
    corr_abs = X_df.corr().abs().values
    np.fill_diagonal(corr_abs, np.nan)              # remove self correlations

    # reward pairs of biomolecules that have similar shifts over conditions (node_conn_abs)
    # and that actually change over conditions (de_score)
    node_conn_abs = np.nanmean(corr_abs, axis = 1)
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
        if k < min_cnt_per_term: continue

        node_scores_abs = np.nanmean(corr_abs[np.ix_(idx, idx)], axis = 1) * de_score[idx]
        mean_conn_abs = np.mean(node_scores_abs)
        pathway_direction = np.nanmedian(feature_direction[idx])

        se = sigma_abs / np.sqrt(k)
        z = (mean_conn_abs - mu_abs) / se
        pval = float(np.exp(norm.logsf(z)))

        results.append([
            pw,
            float(mean_conn_abs),
            float(pval),
            int(k),
            float(pathway_direction),
        ])

    # ------------------------------------------------------------
    # FDR correction + state assignment
    # ------------------------------------------------------------
    if len(results) == 0:
        return []

    pvals = np.array([r[2] for r in results])
    _, qvals, _, _ = multipletests(pvals, method = "fdr_bh")

    for r, q in zip(results, qvals):
        state = assign_state(q, r[4])
        r.append(float(q))   # index 7: qval
        r.append(state)      # index 8: state

    results.sort(key = lambda x: x[5], reverse = True)
    return results


data_frame_processed, _ = preprocess(data_frame)
results = pathway_network_enrichment(data_frame_processed, term_sets, [0, 0, 0, 1, 1, 1])
for r in results:
    print(r[0], r[5], r[4], r[6])

