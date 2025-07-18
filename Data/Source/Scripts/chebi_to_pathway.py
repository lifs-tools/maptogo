import pandas as pd
import zipfile


def a(x):
    return '\t'.join(x)


species = {
    'Saccharomyces cerevisiae': 4932,
    'Escherichia coli': 562,
    'Drosophila melanogaster': 7227,
    'Rattus norvegicus': 10116,
    'Bos taurus': 9913,
    'Homo sapiens': 9606,
    'Caenorhabditis elegans': 6239,
    'Pseudomonas aeruginosa': 287,
    'Arabidopsis thaliana': 3702,
    'Mus musculus': 10090,
}


# 'De Novo Triacylglycerol Biosynthesis TG'
# 'Cardiolipin Biosynthesis CL'
# 'Triacylglycerol Metabolism TG'
# 'Phosphatidylethanolamine Biosynthesis PE'
# 'Phosphatidylcholine Biosynthesis PC'

with zipfile.ZipFile("Data/pathbank_all_metabolites.csv.zip", "r") as zip_ref:
    with zip_ref.open("pathbank_all_metabolites.csv") as f:
        df = pd.read_csv(f, dtype = {"ChEBI ID": str})
df = df[df["ChEBI ID"].notna() & (df["ChEBI ID"] != "")]

df = df[~df["Pathway Name"].str.startswith('De Novo Triacylglycerol Biosynthesis TG')]
df = df[~df["Pathway Name"].str.startswith('Cardiolipin Biosynthesis CL')]
df = df[~df["Pathway Name"].str.startswith('Triacylglycerol Metabolism TG')]
df = df[~df["Pathway Name"].str.startswith('Phosphatidylethanolamine Biosynthesis PE')]
df = df[~df["Pathway Name"].str.startswith('Phosphatidylcholine Biosynthesis PC')]


s = df[["PathBank ID", "Species", "ChEBI ID"]].groupby(["ChEBI ID", "Species"]).agg(a).reset_index()
s["Species"] = s.apply(lambda row: species[row["Species"]], axis = 1)

with open("Data/chebi_to_pathway.csv", "wt") as fout:
    for i, row in s.iterrows():
        fout.write(f"{row["ChEBI ID"]}\t{row["Species"]}\t{row["PathBank ID"]}\n")
