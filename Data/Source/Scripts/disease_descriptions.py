class Term:
    def __init__(self, _id, _name, _relations = None, _synonyms = None, _namespace = None, _xref = None, upper = False, _categories = None):
        if _relations == None: _relations = set()
        if _synonyms == None: _synonyms = set()
        if _namespace == None: _namespace = set()
        if _xref == None: _xref = set()
        if _categories == None: _categories = set()

        self.id = set(_id) if type(_id) in {list, set} else {_id}
        self.name = _name[0].upper() + _name[1:] if upper else _name
        self.relations = set(_relations) if type(_relations) in {list, set} else {_relations}
        self.synonyms = set(_synonyms) if type(_synonyms) in {list, set} else {_synonyms}
        self.namespace = set(_namespace) if type(_namespace) in {list, set} else {_namespace}
        self.xref = set(_xref) if type(_xref) in {list, set} else {_xref}
        self.categories = set(_categories) if type(_categories) in {list, set} else {_categories}
        self.description = ""



def get_terms(term_file, prefix, upper = False, with_relationship = False):
    ontology_terms = {}

    opener = open if term_file.lower()[-3:] != ".gz" else gzip.open

    with opener(term_file, "rt") as obo:
        term_id, name, relations, synonyms, namespace, is_obsolete, xref, description = "", "", set(), [], "", False, set(), ""
        for i, line in enumerate(obo):
            line = line.strip()
            if len(line) == 0: continue
            if i > 0 and i % 10000 == 0: print(term_file, i)

            if line == "[Term]":
                if len(term_id) > 0 and len(name) > 0 and not is_obsolete and term_id.startswith(prefix):
                    ontology_terms[term_id] = Term(term_id, name, relations, synonyms, namespace, xref, upper)
                    ontology_terms[term_id].description = description
                term_id, name, relations, synonyms, is_obsolete, namespace, xref, description = "", "", set(), [], False, "", set(), ""

            elif line.startswith("id:"):
                term_id = line[3:].strip(" ")

            elif line.startswith("name:"):
                name = line[5:].strip(" ")

            elif line.startswith("is_a:"):
                relation = line[5:].split(" ! ")[0].strip(" ").split(" ")[0]
                relations.add(relation)

            elif with_relationship and line.startswith("relationship: "):
                relation = prefix + line[14:].split(prefix)[1].split(" ! ")[0]
                relations.add(relation)

            elif line.startswith( "synonym:"):
                synonyms.append(line.split("\"")[1])

            elif line.startswith("namespace: "):
                namespace = line[11:]

            elif line.startswith("def: "):
                description = t[1] if len((t := line.split("\""))) > 2 else ""

            elif line.startswith("xref: "):
                xref.add(line[6:].split(" ")[0])

            elif line.startswith("is_obsolete:"):
                is_obsolete = True

        if len(term_id) > 0 and len(name) > 0 and not is_obsolete and term_id.startswith(prefix):
            ontology_terms[term_id] = Term(term_id, name, relations, synonyms, namespace, xref, upper)
            ontology_terms[term_id].description = description

    return ontology_terms


doid_terms = get_terms("Data/doid-base.obo", "DOID:", upper = True)
mondo_terms = get_terms("Data/mondo.obo", "MONDO:", upper = True)

diseases = {}
for _, term in doid_terms.items():
    name = term.name.lower()
    if name not in diseases: diseases[name] = (term.name, term.description)

for _, term in mondo_terms.items():
    name = term.name.lower()
    if name not in diseases: diseases[name] = (term.name, term.description)

with open("Data/disease_descriptions.csv", "wt") as output_stream:
    for _, (name, description) in diseases.items():
        output_stream.write(f"{name}\t{description}\n")

