import gzip

class Term:
    def __init__(self, _id, _name, _relations_conjugate, _relations, _inchikey = None):
        self.id = _id
        self.name = _name
        self.relations_conjugate = _relations_conjugate
        self.relations = _relations
        self.inchikey = _inchikey


term_id, name, relations_conjugate, relations, inchikey = "", "", set(), set(), None
ontology_terms = {}
with gzip.open("Data/chebi.obo.gz", "rt") as obo:
    for i, line in enumerate(obo.read().split("\n")):
        if len(line) == 0: continue

        if line == "[Term]":
            if len(term_id) > 0 and len(name) > 0:
                term = Term(term_id, name, relations_conjugate, relations, inchikey)
                ontology_terms[term_id] = term

            term_id, name, relations_conjugate, relations, inchikey = "", "", set(), set(), None

        elif line[:3] == "id:":
            term_id = line[3:].strip(" ").replace("CHEBI:", "")

        elif line[:5] == "name:":
            name = line[5:].strip(" ").replace("CHEBI:", "")

        elif line[:35] == "relationship: is_conjugate_base_of ":
            relation = line[35:].strip(" ")
            relations_conjugate.add(relation.replace("CHEBI:", ""))

        # elif line.startswith("is_a:"):
        #     relations.add(line[6:].replace("CHEBI:", ""))

        elif line.startswith("property_value: http://purl.obolibrary.org/obo/chebi/inchikey"):
            inchikey = line[61:].split("\"")[1].split("-")[0]


if len(term_id) > 0 and len(name) > 0:
    ontology_terms[term_id] = Term(term_id, name, relations_conjugate, relations, inchikey)

for term_id, term in ontology_terms.items():
    for relation in term.relations_conjugate:
        ontology_terms[relation].relations_conjugate.add(term_id)

    for relation in term.relations:
        if relation in ontology_terms:
            if ontology_terms[relation].inchikey != None and ontology_terms[relation].inchikey == term.inchikey:
                ontology_terms[relation].relations.add(term_id)


with open("Data/chebi.csv", "wt") as out:
    for term_id, term in ontology_terms.items():
        term.relations |= term.relations_conjugate
        out.write(f"{term_id}\t{term.name}\t{'\t'.join(term.relations)}\n")
