import gzip

i = 0
uniprot_to_go = {}
with gzip.open("Data/goa_uniprot_gcrp.gpa.gz", "rt") as input_stream:
    for line in input_stream:
        if (i := i + 1) % 1e6 == 0: print(i)
        if line[0] == "!": continue
        tokens = line.strip().split("\t")
        if len(tokens) < 5 or tokens[1] not in uniport_accs: continue

        if tokens[1] not in uniprot_to_go: uniprot_to_go[tokens[1]] = {tokens[3]}
        else: uniprot_to_go[tokens[1]].add(tokens[3])

print("len", len(uniprot_to_go))
with gzip.open("Data/goa_uniprot.csv.gz", "wt") as output_stream:
    for acc_id, uniprot_ids in uniprot_to_go.items():
        output_stream.write(f"{acc_id}\t{'\t'.join(uniprot_ids)}\n")
