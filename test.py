from EnrichmentDataStructure import EnrichmentOntology, current_path, SessionEntry
#from pygoslin.domain.LipidFaBondType import LipidFaBondType
#from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import LipidParser

lipid_parser = LipidParser()
session = SessionEntry()
ontology = EnrichmentOntology(f"{current_path}/Data/ontology_10090.gz", lipid_parser = lipid_parser)

test_data = {}
lipid_dict = {}
with open(f"{current_path}/Test/test_data.csv", "rt") as infile:
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 2: continue

        if tokens[0] not in test_data:
            test_data[tokens[0]] = set()
            lipid_dict[tokens[0]] = lipid_parser.parse(tokens[0])
        test_data[tokens[0]].add(tokens[1])



ontology.set_background(session, lipid_dict = lipid_dict)
unit_tests = 0
for lipid_name, term_ids in test_data.items():
    for term_id in term_ids:
        if term_id not in session.search_terms:
            print("ERROR:", lipid_name, term_id)
            exit(-1)
        unit_tests += 1

print(f"Test passed ({unit_tests} unit tests) without errors")
