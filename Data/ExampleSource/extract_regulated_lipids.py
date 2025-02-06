import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def extract_significant(df, start_data, variable, condition_1, condition_2, pval_ths = 0.05, fc_ths = 1, filters = None):

    if filters != None:
        for filter_column, filter_value in filters:
            df = df[df[filter_column] == filter_value]

    data = df.iloc[:, start_data:]

    print(list(data.columns))
    print("=" * 20)
    keep_lipids = (data[df[variable] == condition_1].mean().notna()) & (data[df[variable] == condition_2].mean().notna())
    data = data.loc[:, keep_lipids]

    keep_lipids = (data[df[variable] == condition_1].notna().sum(axis = 0) > 1) & (data[df[variable] == condition_2].notna().sum(axis = 0) > 1)
    data = data.loc[:, keep_lipids]

    c1_data = data[df[variable] == condition_1]
    c2_data = data[df[variable] == condition_2]

    t, p = ttest_ind(c1_data, c2_data, nan_policy = "omit")
    p_corr = pd.Series(multipletests(p, method = "fdr_bh")[1], index = data.columns)

    fc = np.abs(np.log2(c2_data.mean() / c1_data.mean()))
    print(list(data.loc[: , (p_corr < pval_ths) & (fc > fc_ths)].columns))

# df = pd.read_excel("AOValves.xlsx")
# extract_significant(df, 2, "Condition", "Healthy", "Calcified")

df = pd.read_excel("Mouse.xlsx")
extract_significant(df, 4, "Treatment", "Unst", "5CRP")

#df = pd.read_excel("Contraceptives.xlsx")
#extract_significant(df, 3, "Condition", "Female no CC", "Female CC")

# df = pd.read_excel("MK.xlsx")
# df.columns = [col.strip() for col in df.columns]
# extract_significant(df, 6, "Day", 1, 3)
