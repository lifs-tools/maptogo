import gzip

uniprot_to_term = {}

with gzip.open("Data/uniprot.csv.gz") as infile:
    for i, line in enumerate(infile):
        if i == 0: continue
        tokens = line.decode("utf8").split("\t")
        if len(tokens) < 6 or len(tokens[5]) == 0: continue

        for t in tokens[5].split('";'):
            t = t.strip('"').split("; ")
            if len(t) < 2 or not t[1].startswith("phenotype"): continue
            if tokens[0] not in uniprot_to_term: uniprot_to_term[tokens[0]] = []
            uniprot_to_term[tokens[0]].append(f"OMIM:{t[0]}")


class Term:
    def __init__(self, _id, _name, _relations = None, _synonyms = None, _namespace = "", _xref = None):

        if _relations == None: _relations = set()
        if _synonyms == None: _synonyms = []
        self.id = _id
        self.name = _name
        self.relations = _relations
        self.synonyms = _synonyms
        self.namespace = _namespace
        self.xref = _xref

ontology_terms = []
with open("Data/doid-base.obo", "rt") as obo:
    term_id, name, relations, synonyms, namespace, xref = "", "", set(), [], "", []
    for i, line in enumerate(obo):
        line = line.strip()
        if len(line) == 0: continue
        if i > 0 and i % 10000 == 0: print("DOID:", i)

        if line == "[Term]":
            if len(term_id) > 0 and len(name) > 0 and len(xref) > 0:
                term = Term(term_id, name, relations, synonyms, namespace, xref)
                ontology_terms.append(term)
            term_id, name, relations, synonyms, namespace, xref = "", "", set(), [], "", []

        elif line[:3] == "id:":
            term_id = line[3:].strip(" ")

        elif line[:5] == "name:":
            name = line[5:].strip(" ")

        elif line[:5] == "xref:":
            xref.append(line[5:].strip(" "))


        elif line[:5] == "is_a:":
            relation = line[5:].split(" ! ")[0].strip(" ")
            relations.add(relation)

        elif line[:8] == "synonym:":
            synonyms.append(line.split("\"")[1])

        elif line[:11] == "namespace: ":
            namespace = line[11:]


    if len(term_id) > 0 and len(name) > 0 and len(xref) > 0:
        term = Term(term_id, name, relations, synonyms, namespace, xref)
        ontology_terms.append(term)


omim_to_doid = {}
for term in ontology_terms:
    omim = ""
    for xref in term.xref:
        if xref[:4] == "MIM:": omim = "O" + xref

    if omim != "":
        omim_to_doid[omim] = term.id


for uniprot, terms in uniprot_to_term.items():
    uniprot_to_term[uniprot] = [omim_to_doid[term] for term in terms if term in omim_to_doid]
uniprot_to_term = {k: v for k, v in uniprot_to_term.items() if len(v) > 0}


with open("Data/uniprot_to_doid.csv", "wt") as out:
    for k, v in uniprot_to_term.items():
        out.write(f"{k}\t{'\t'.join(v)}\n")

