# Systematic regulation score

#from dataframe import *
import numpy as np
import pandas as pd
from time import time
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
from collections import Counter
import matplotlib.pyplot as plt
from EnrichmentDataStructure import *

#print(f"loaded {len(term_sets)} terms")



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



def pathway_network_enrichment(X_df, term_sets_idx, conditions, min_cnt_per_grp = 2, min_cnt_per_term = 20):
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
    overall_mean = np.nanmean(X, axis=0)

    ss_between = 0
    ss_within = 0

    for g in groups:
        Xg = X[conditions == g]
        num_g = Xg.shape[0]
        n_valid = np.sum(~np.isnan(Xg), axis = 0)  # per-feature valid count
        mean_g = np.nanmean(Xg, axis = 0)
        var_g = np.where(n_valid > 2, np.nanvar(Xg, axis = 0, ddof = 1), 0.0)
        ss_between += num_g * (mean_g - overall_mean) ** 2
        ss_within += (num_g - 1) * var_g

    de_score = ss_between / (ss_between + ss_within)

    # DE direction: linear regression slope of group means vs condition label
    x = groups.astype(float) - groups.astype(float).mean()
    feature_direction = (x @ group_means) / (x @ x)

    # ------------------------------------------------------------
    # correlation matrix
    # ------------------------------------------------------------
    corr_abs = X_df.corr().abs().values             # abs: reward anti-correlation, too
    np.fill_diagonal(corr_abs, np.nan)              # remove self correlations

    # reward pairs of biomolecules that have similar shifts over conditions (biomol_conn_abs)
    # and that actually change over conditions (de_score)
    biomol_conn_abs = np.nanmean(corr_abs, axis = 1)
    combined_score_abs = biomol_conn_abs * de_score


    # ------------------------------------------------------------
    # pathway mapping
    # ------------------------------------------------------------
    # term_sets_idx = {
    #     pw: np.array([feature_idx[b] for b in biomols if b in feature_idx], dtype = np.int64)
    #     for pw, biomols in term_sets.items()
    # }
    term_sets_idx = {k: v for k, v in term_sets_idx.items() if len(v) >= min_cnt_per_term}

    # global statistics (for z-score baseline)
    mu_abs = np.mean(combined_score_abs)
    sigma_abs = np.std(combined_score_abs, ddof = 1)

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

        biomol_scores_abs = np.nanmean(corr_abs[np.ix_(idx, idx)], axis = 1) * de_score[idx]
        pathway_direction = np.nanmedian(feature_direction[idx])

        # this value needs to be normally distributed over all terms according
        # to the Central Limit Therorem (CLT) to apply the p-value computation
        mean_conn_abs = np.mean(biomol_scores_abs)

        # z-score
        z = (mean_conn_abs - mu_abs) * np.sqrt(k) / sigma_abs
        # one sided test (greater)
        pvalue = float(np.exp(norm.logsf(z)))

        results.append([
            pw,
            float(mean_conn_abs),
            float(pvalue),
            int(k),
            float(pathway_direction),
        ])


    # ------------------------------------------------------------
    # FDR correction + state assignment
    # ------------------------------------------------------------
    if len(results) == 0:
        return []

    pvalues = np.array([r[2] for r in results])
    _, qvals, _, _ = multipletests(pvalues, method = "fdr_bh")

    for r, q in zip(results, qvals):
        state = assign_state(q, r[4])
        r.append(float(q))   # index 5: qval
        r.append(state)      # index 6: state

    results.sort(key = lambda x: x[5], reverse = True)
    return results




# data_frame = pd.read_excel("AOValves-comlete.xlsx")
# data_frame = data_frame.set_index("Condition").iloc[:, 1:]
# data_frame = data_frame.iloc[list(range(10)) + list(range(20, 30))]
# #conditions = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
# conditions = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

from itertools import groupby

data_frame = pd.read_excel("MK.xlsx", index_col = "Day").iloc[:, 2:]
conditions = [i for i, (_, g) in enumerate(groupby(data_frame.index)) for _ in g]



data_frame_processed, _ = preprocess(data_frame)
feature_idx = {f: i for i, f in enumerate(data_frame_processed.columns)}



# possible domains: "Biological process", "Cellular component", "Disease", "Metabolic and signalling pathway", "Molecular function", "Phenotype", "Physical or chemical properties", "Enzymatic activity (Swiss-Prot)", "Enzymatic activity (Swiss-Prot + TrEMBL)", or any combination
domains = {"Metabolic and signalling pathway"}

omics_lists = [
    (bm := list(data_frame_processed.columns)),
    bm,
    [], # up-regulated_lipids_list
    [], # down-regulated_lipids_list
    bm,
    bm,
    [], # up-regulated_proteins_list
    [], # down-regulated_proteins_list
    [],
    [],
    [], # up-regulated_metabolites_list
    [], # down-regulated_metabolites_list
    [], # background_transcripts_list
    [], # regulated_transcripts_list
    [], # up-regulated_transcripts_list
    [], # down-regulated_transcripts_list
]

print("Loading ontology")
ontology = EnrichmentOntology(f"Data/ontology_10090.zst")

## parse all biomolocule entries
print("Parse all biomolecule entries")
try:
    (
        target_set,
        lipidome,
        regulated_lipids,
        upregulated_lipids,
        downregulated_lipids,
        proteome,
        regulated_proteins,
        upregulated_proteins,
        downregulated_proteins,
        metabolome,
        regulated_metabolites,
        upregulated_metabolites,
        downregulated_metabolites,
        transcriptome,
        regulated_transcripts,
        upregulated_transcripts,
        downregulated_transcripts,
        background_list,
    ) = check_user_input(
        False,
        [True, True, False, False],
        omics_lists,
        ontology,
        MOLECULE_HANDLING_REMOVE,
        MOLECULE_HANDLING_REMOVE,
    )
except Exception as error_message:
    print(error_message)
    exit(-1)



## set enrichment background
print("Get term associated biomolecules")
(
    search_terms,
    all_parent_nodes,
) = ontology.set_background(
    lipid_dict = lipidome,
    protein_set = proteome,
    metabolite_set = metabolome,
    transcript_set = transcriptome,
    use_bounded_fatty_acyls = False,
)



term_sets_idx = {
    pw.name: np.array([feature_idx[b] for b in biomols if b in feature_idx], dtype = np.int64)
    for pw, biomols in search_terms.items() if len(set(ontology.get_domains(pw.domains)) & domains) > 0
}






results = pathway_network_enrichment(
    data_frame_processed,
    term_sets_idx,
    conditions,
)
for r in results:
    print(r[0], r[5], r[4], r[6])

