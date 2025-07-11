class Term:
    def __init__(self, _id, _name, _relations):
        self.id = _id
        self.name = _name
        self.relations = _relations


term_id, name, relations = "", "", set()
ontology_terms = {}
with open("Data/chebi_lite.obo", "rt") as obo:
    for i, line in enumerate(obo):
        line = line.strip()
        if len(line) == 0: continue

        if line == "[Term]":
            if len(term_id) > 0 and len(name) > 0:
                term = Term(term_id, name, relations)
                ontology_terms[term_id] = term

            term_id, name, relations = "", "", set()

        elif line[:3] == "id:":
            term_id = line[3:].strip(" ").replace("CHEBI:", "")

        elif line[:5] == "name:":
            name = line[5:].strip(" ").replace("CHEBI:", "")

        elif line[:35] == "relationship: is_conjugate_base_of ":
            relation = line[35:].strip(" ")
            relations.add(relation.replace("CHEBI:", ""))


if len(term_id) > 0 and len(name) > 0:
    term = Term(term_id, name, relations)
    ontology_terms[term_id] = term

for term_id, term in ontology_terms.items():
    for relation in term.relations:
        ontology_terms[relation].relations.add(term_id)

with open("Data/chebi.csv", "wt") as out:
    for term_id, term in ontology_terms.items():
        out.write(f"{term_id}\t{term.name}\t{'\t'.join(term.relations)}\n")
