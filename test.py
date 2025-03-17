from EnrichmentDataStructure import EnrichmentOntology, current_path, SessionEntry
#from pygoslin.domain.LipidFaBondType import LipidFaBondType
#from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import LipidParser

lipid_parser = LipidParser()
session = SessionEntry()
ontology = EnrichmentOntology(f"{current_path}/Data/ontology_10090.gz", lipid_parser = lipid_parser)

test_data = {}
lipid_dict = {}
protein_set = set()
metabolite_set = set()
with open(f"{current_path}/Test/test_data.csv", "rt") as infile:
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 3 or (len(tokens[0]) and tokens[0][0] == '#') or tokens[0] not in {"M", "L", "P"}: continue
        if tokens[1] not in test_data:
            test_data[tokens[1]] = set()
            if tokens[0] == "L": lipid_dict[tokens[1]] = lipid_parser.parse(tokens[1])
            if tokens[0] == "M": metabolite_set.add(tokens[1])
            if tokens[0] == "P": protein_set.add(tokens[1])
        test_data[tokens[1]].add(tokens[2])


ontology.set_background(session, lipid_dict = lipid_dict, protein_set = protein_set, metabolite_set = metabolite_set)
unit_tests = 0
for molecule_name, term_ids in test_data.items():
    for term_id in term_ids:
        if term_id not in session.search_terms:
            print("ERROR:", molecule_name, term_id)
            exit(-1)
        unit_tests += 1

print(f"Test passed ({unit_tests} unit tests) without errors")
