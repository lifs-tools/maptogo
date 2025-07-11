import gzip

hgnc_to_mondo = {}
with gzip.open("Data/causal_gene_to_disease_association.all.tsv.gz") as infile:
    for i, line in enumerate(infile):
        if i == 0: continue
        tokens = line.decode("utf8").split("\t")
        if len(tokens) < 8 or len(tokens[7]) == 0: continue
        if tokens[0] not in hgnc_to_mondo: hgnc_to_mondo[tokens[0]] = set()
        hgnc_to_mondo[tokens[0]].add(tokens[7])

with gzip.open("Data/correlated_gene_to_disease_association.all.tsv.gz") as infile:
    for i, line in enumerate(infile):
        if i == 0: continue
        tokens = line.decode("utf8").split("\t")
        if len(tokens) < 8 or len(tokens[7]) == 0: continue
        if tokens[0] not in hgnc_to_mondo: hgnc_to_mondo[tokens[0]] = set()
        hgnc_to_mondo[tokens[0]].add(tokens[7])

with open("Data/hgnc_to_mondo.csv", "wt") as out:
    for hgnc, mondo in hgnc_to_mondo.items():
        out.write(f"{hgnc}\t{'\t'.join(mondo)}\n")
