import gzip

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


    def get_first_id(self):
        return list(self.id)[0]

    def joining(self, a_set):
        return "|".join(sorted(list(a_set)))

    def to_string(self):
        return f"{self.joining(self.id)}\t{self.name}\t{self.joining(self.relations)}\t{self.joining(self.synonyms)}\t{self.joining(self.namespace)}\t{self.joining(self.categories)}\n"


    def copy(self):
        return Term(set(self.id), self.name, set(self.relations), set(self.synonyms), set(self.namespace), set(self.xref), _categories = set(self.categories))


    def merge(self, copy):
        self.id |= copy.id
        self.relations |= copy.relations
        self.synonyms |= copy.synonyms
        self.namespace |= copy.namespace
        self.xref |= copy.xref
        self.categories |= copy.categories


def get_terms(term_file, prefix, upper = False):
    ontology_terms = {}

    opener = open if term_file.lower()[-3:] != ".gz" else gzip.open

    with opener(term_file, "rt") as obo:
        term_id, name, relations, synonyms, namespace, is_obsolete, xref = "", "", set(), [], "", False, set()
        for i, line in enumerate(obo):
            line = line.strip()
            if len(line) == 0: continue
            if i > 0 and i % 100000 == 0: print(term_file, i)

            if line == "[Term]":
                if len(term_id) > 0 and len(name) > 0 and not is_obsolete and term_id.startswith(prefix):
                    ontology_terms[term_id] = Term(term_id, name, relations, synonyms, namespace, xref, upper)
                term_id, name, relations, synonyms, is_obsolete, namespace, xref = "", "", set(), [], False, "", set()

            elif line.startswith("id:"):
                term_id = line[3:].strip(" ")

            elif line.startswith("name:"):
                name = line[5:].strip(" ")

            elif line.startswith("is_a:"):
                relation = line[5:].split(" ! ")[0].strip(" ").split(" ")[0]
                relations.add(relation)

            elif line.startswith( "synonym:"):
                synonyms.append(line.split("\"")[1])

            elif line.startswith("namespace: "):
                namespace = line[11:]

            elif line.startswith("xref: "):
                xref.add(line[6:].split(" ")[0])

            elif line.startswith("is_obsolete:"):
                is_obsolete = True

        if len(term_id) > 0 and len(name) > 0 and not is_obsolete and term_id.startswith(prefix):
            ontology_terms[term_id] = Term(term_id, name, relations, synonyms, namespace, xref, upper)

    return ontology_terms

chebi_obo = get_terms("Data/chebi.obo.gz", "CHEBI:")
chebi_synonyms = {}
for term_id, term in chebi_obo.items():
    for synonym in term.synonyms:
        chebi_synonyms[synonym] = term_id

with gzip.open("../CHEBI_synonyms.csv.gz", "wt") as out:
    for synonym, term_id in chebi_synonyms.items():
        out.write(f"{synonym}\t{term_id}\n")
