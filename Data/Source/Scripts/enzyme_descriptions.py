import gzip

enzymes = {}

with gzip.open("Data/uniprot.csv.gz", "rt") as input_stream:
    line = input_stream.readline()
    for line in input_stream:
        tokens = line.strip().split("\t")
        if len(tokens) < 2: continue
        gene_name = tokens[2]
        description = tokens[12].replace("FUNCTION: ", "") if len(tokens) > 12 else ""
        if gene_name not in enzymes or len(description) > len(enzymes[gene_name]):
            enzymes[gene_name] = description


with open("Data/enzyme_descriptions.csv", "wt") as output_stream:
    for gene_name, description in enzymes.items():
        output_stream.write(f"{gene_name}\t{description}\n")
