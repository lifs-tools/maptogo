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

    keep_lipids = (data[df[variable] == condition_1].notna().sum(axis = 0) > 1) & (data[df[variable] == condition_2].notna().sum(axis = 0) > 1)
    data = data.loc[:, keep_lipids]

    c1_data = data[df[variable] == condition_1]
    c2_data = data[df[variable] == condition_2]

    t, p_corr = ttest_ind(c1_data, c2_data, nan_policy = "omit")
    p_corr = pd.Series(multipletests(p_corr, method = "fdr_bh")[1], index = data.columns)

    log_fc = np.log2(c2_data.mean() / c1_data.mean())
    thresholds = (p_corr <= pval_ths)
    if up_down == "both": thresholds &= (np.abs(log_fc) >= log_fc_ths)
    elif up_down == "down": thresholds &= (log_fc <= -log_fc_ths)
    else: thresholds &= (log_fc >= log_fc_ths)

    print(list(data.loc[: , thresholds].columns))

# df = pd.read_excel("Simplex2016.xlsx", "Lipids")
# extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "both")
# print("=" * 20)
# extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "up")
# print("=" * 20)
# extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "down")

# df = pd.read_excel("Simplex2016.xlsx", "Proteins")
# extract_significant(df, 5, "Condition", "Dmso", "Rosi", log_fc_ths = 1, up_down = "down")


df = pd.read_excel("Simplex2016.xlsx", "Metabolites")
extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "both")
print("=" * 20)
extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "up")
print("=" * 20)
extract_significant(df, 3, "Condition", "Dmso", "Rosi", up_down = "down")




# df = pd.read_excel("AOValves-All.xlsx")
# extract_significant(df, 2, "Condition", "Healthy", "Calcified", log_fc_ths = np.log2(1.5), up_down = "up")

#df = pd.read_excel("AOValves-Proteomics.xlsx")
#extract_significant(df, 2, "Condition", "Healthy", "Calcified")

# df = pd.read_excel("Mouse.xlsx")
# extract_significant(df, 4, "Treatment", "Unst", "5CRP")

# df = pd.read_excel("Contraceptives.xlsx")
# extract_significant(df, 3, "Condition", "Female no CC", "Female CC", log_fc_ths = 0.5)

# df = pd.read_excel("MK_Proteomics.xlsx")
# df.columns = [col.strip() for col in df.columns]
# extract_significant(df, 3, "Day", 0, 7, log_fc_ths = 2) #, up_down = "down")

# df = pd.read_excel("Heart-reperfusion-Metabolomics-data.xlsx")
# extract_significant(df, 5, "Time", "0h", "I2h", log_fc_ths = 0.5,  filters = [("State", "MI"), ("Tissue", "Heart"), ("Group", "MI progression")])

