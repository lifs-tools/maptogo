import gzip

with open("Data/uniprot_to_go.csv", "wt") as out:
    with gzip.open('Data/uniprot.csv.gz','r') as fin:
        for i, line in enumerate(fin):
            if i == 0: continue
            if i % 10000 == 0: print(f"uniprot_to_go: {i}")
            tokens = line.decode().strip().split("\t")
            if len(tokens) < 4: continue

            if len(tokens[3]) > 0 and tokens[3].find("[") > -1:
                go_terms = []
                for entry in tokens[3].split("["):
                    go_term = entry.split("]")[0]
                    if len(go_term) != 10 or go_term[:3] != "GO:": continue
                    go_terms.append(go_term)
                tokens[3] = ",".join(go_terms)

            out.write(f"{'\t'.join(tokens[:4])}\n")







