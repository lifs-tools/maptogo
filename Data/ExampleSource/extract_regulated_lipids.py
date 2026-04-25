import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests


def extract_significant(df, start_data, variable, condition_1, condition_2, pval_ths = 0.05, log_fc_ths = 1, filters = None, up_down = "both"):

    if filters != None:
        for filter_column, filter_value in filters:
            df = df[df[filter_column] == filter_value]

    data = df.iloc[:, start_data:]

    print(len(data.columns))
    print(list(data.columns))
    print("-" * 20)

    keep_lipids = (data[df[variable] == condition_1].mean().notna()) & (data[df[variable] == condition_2].mean().notna())
    data = data.loc[:, keep_lipids]
    data[data == 0] = np.nan

    keep_lipids = (data[df[variable] == condition_1].notna().sum(axis = 0) > 1) & (data[df[variable] == condition_2].notna().sum(axis = 0) > 1)
    data = data.loc[:, keep_lipids]

    c1_data = data[df[variable] == condition_1]
    c2_data = data[df[variable] == condition_2]


    t, p_corr = ttest_ind(np.log2(c1_data), np.log2(c2_data), nan_policy = "omit", equal_var = "False")
    p_corr = multipletests(p_corr, method = "fdr_bh")[1]
    p_corr = pd.Series(p_corr, index = data.columns)

    log_fc = np.log2(c2_data.mean() / c1_data.mean())
    thresholds = (p_corr <= pval_ths)
    if up_down == "both": thresholds &= (np.abs(log_fc) >= log_fc_ths)
    elif up_down == "down": thresholds &= (log_fc <= -log_fc_ths)
    else: thresholds &= (log_fc >= log_fc_ths)

    #print(list(data.loc[: , thresholds].columns))

    print(log_fc[log_fc.index.isin(["Cer 34:2", "Cer 34:1", "Cer 36:1", "Cer 40:2", "Cer 40:1", "Cer 42:3", "Cer 42:2", "Cer 42:1", "DG 30:0", "DG 32:0", "DG 32:1", "DG 34:0", "DG 34:1", "DG 34:2", "DG 36:1", "DG 36:2", "DG 36:4", "DG 38:5", "DG 43:6", "DG 46:1", "DG 48:1", "DG 50:1", "LPA 16:0", "LPA 18:3", "LPA 18:0", "PA 32:0", "PA 32:1", "PA 34:0", "PA 34:1", "PA 34:2", "PA 35:2", "PA 36:2", "PA 36:1", "PA 38:4"])])
    print(p_corr[p_corr.index.isin(["Cer 34:2", "Cer 34:1", "Cer 36:1", "Cer 40:2", "Cer 40:1", "Cer 42:3", "Cer 42:2", "Cer 42:1", "DG 30:0", "DG 32:0", "DG 32:1", "DG 34:0", "DG 34:1", "DG 34:2", "DG 36:1", "DG 36:2", "DG 36:4", "DG 38:5", "DG 43:6", "DG 46:1", "DG 48:1", "DG 50:1", "LPA 16:0", "LPA 18:3", "LPA 18:0", "PA 32:0", "PA 32:1", "PA 34:0", "PA 34:1", "PA 34:2", "PA 35:2", "PA 36:2", "PA 36:1", "PA 38:4"])])
    #print(log_fc[log_fc.index.str.startswith("PA")])


df = pd.read_excel("Simplex2016.xlsx", "Lipids")
extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "both")
# print("=" * 20)
# extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "up")
# print("=" * 20)
# extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "down")

#df = pd.read_excel("Simplex2016.xlsx", "Proteins")
# extract_significant(df, 5, "Condition", "Dmso", "Rosi", log_fc_ths = 1, up_down = "down")
# extract_significant(df, 5, "Condition", "Dmso", "Rosi", log_fc_ths = 1, up_down = "both")


#df = pd.read_excel("Simplex2016.xlsx", "Metabolites")
#extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "both")
# print("=" * 20)
# extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "up")
# print("=" * 20)
# extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "down")




# df = pd.read_excel("AOValves-All.xlsx")
# extract_significant(df, 2, "Condition", "Healthy", "Calcified", log_fc_ths = np.log2(1.5), up_down = "both")

# df = pd.read_excel("AOValves-Proteomics.xlsx")
# extract_significant(df, 2, "Condition", "Healthy", "Calcified", up_down = "down")


# df = pd.read_excel("AOValves-Metabolites.xlsx")
# extract_significant(df, 2, "Condition", "Healthy", "Calcified", log_fc_ths = np.log2(1.5), up_down = "down")


# df = pd.read_excel("Mouse.xlsx")
# extract_significant(df, 4, "Treatment", "Unst", "5CRP")




# df = pd.read_excel("Contraceptives.xlsx")
# extract_significant(df, 3, "Condition", "Female no CC", "Female CC", log_fc_ths = 0.5)




# df = pd.read_excel("MK_Proteomics.xlsx")
# df.columns = [col.strip() for col in df.columns]
# extract_significant(df, 3, "Day", 0, 7, log_fc_ths = 2, up_down = "down")

# df = pd.read_excel("MK_Lipidomics.xlsx")
# df.columns = [col.strip() for col in df.columns]
# extract_significant(df, 3, "Day", 1, 7, log_fc_ths = np.log2(2), up_down = "down")




# df = pd.read_excel("Heart-reperfusion-Metabolomics-data.xlsx")
# extract_significant(df, 5, "Time", "0h", "I2h", log_fc_ths = 0.5,  filters = [("State", "MI"), ("Tissue", "Heart"), ("Group", "MI progression")])

