import pandas as pd
import zipfile

with zipfile.ZipFile("Data/pathbank_all_proteins.csv.zip", "r") as zip_ref:
    with zip_ref.open("pathbank_all_proteins.csv") as f:
        df = pd.read_csv(f)

df = df[~df["Uniprot ID"].isna()]
df = df[~df["Pathway Name"].isna()]
df = df[~df["PathBank ID"].isna()]

df = df[~df["Pathway Name"].str.startswith('De Novo Triacylglycerol Biosynthesis TG')]
df = df[~df["Pathway Name"].str.startswith('Cardiolipin Biosynthesis CL')]
df = df[~df["Pathway Name"].str.startswith('Triacylglycerol Metabolism TG')]
df = df[~df["Pathway Name"].str.startswith('Phosphatidylethanolamine Biosynthesis PE')]
df = df[~df["Pathway Name"].str.startswith('Phosphatidylcholine Biosynthesis PC')]

pathbank_ids = list(df["PathBank ID"])
pathway_names = list(df["Pathway Name"])
uniprot_ids = list(df["Uniprot ID"])

pathbank = {}
for pathbank_id, pathway_name, uniprot_id in zip(pathbank_ids, pathway_names, uniprot_ids):
    if pathbank_id not in pathbank: pathbank[pathbank_id] = [pathway_name, set()]
    pathbank[pathbank_id][1].add(uniprot_id)

with open("Data/pathbank_to_uniprot.csv", "wt") as out:
    for r_id, r_data in pathbank.items():
        out.write(f"{r_id}\t{r_data[0]}\t{'\t'.join(r_data[1])}\n")
