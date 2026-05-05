import gzip

class Term:
    def __init__(self, _id, _name, _relations_conjugate, _inchikey = None):
        self.id = _id
        self.name = _name
        self.relations_conjugate = _relations_conjugate
        self.inchikey = _inchikey


term_ids, name, relations_conjugate, inchikey = [], "", set(), None
ontology_terms = {}
with gzip.open("Data/chebi.obo.gz", "rt") as obo:
    for i, line in enumerate(obo.read().split("\n")):
        if len(line) == 0: continue

        if line == "[Term]":
            if len(term_ids) > 0 and len(name) > 0:
                term = Term(term_ids, name, relations_conjugate, inchikey)
                for term_id in term_ids: ontology_terms[term_id] = term

            term_ids, name, relations_conjugate, inchikey = [], "", set(), None

        elif line[:3] == "id:":
            term_ids.append(line[3:].strip(" ").replace("CHEBI:", ""))

        elif line[:5] == "name:":
            name = line[5:].strip(" ").replace("CHEBI:", "")

        elif line[:35] == "relationship: is_conjugate_base_of ":
            relation = line[35:].strip(" ")
            relations_conjugate.add(relation.replace("CHEBI:", ""))

        elif line.startswith("alt_id:"):
            term_ids.append(line.replace("alt_id: CHEBI:", "").strip())

        elif line.startswith("property_value: http://purl.obolibrary.org/obo/chebi/inchikey"):
            inchikey = line[61:].split("\"")[1].split("-")[0]


if len(term_ids) > 0 and len(name) > 0:
    for term_id in term_ids: ontology_terms[term_id] = Term(term_id, name, relations_conjugate, inchikey)

for term_id, term in ontology_terms.items():
    for relation in term.relations_conjugate:
        ontology_terms[relation].relations_conjugate.add(term_id)

added = set()
with open("Data/chebi.csv", "wt") as out:
    for term_id, term in ontology_terms.items():
        if term.id[0] in added: continue
        added.add(term.id[0])
        out.write(f"{'\t'.join(term.id)}\t{term.name}\t{'\t'.join(term.relations_conjugate)}\n")
