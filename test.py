from EnrichmentDataStructure import EnrichmentOntology, current_path, SessionEntry
from pygoslin.parser.Parser import LipidParser
import time

lipid_parser = LipidParser()
session = SessionEntry()

test_data = {}
lipid_dict = {}
protein_set = set()
metabolite_set = set()
with open(f"{current_path}/Test/test_data.csv", "rt") as infile:
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 4 or \
            (len(tokens[0]) and tokens[0][0] == '#') or \
            tokens[0] not in {"9606", "10090"} or \
            tokens[1] not in {"M", "L", "P"}: continue
        if tokens[0] not in test_data: test_data[tokens[0]] = {}
        if tokens[2] not in test_data[tokens[0]]:
            test_data[tokens[0]][tokens[2]] = set()
            if tokens[1] == "L": lipid_dict[tokens[2]] = lipid_parser.parse(tokens[2])
            if tokens[1] == "M": metabolite_set.add(tokens[2])
            if tokens[1] == "P": protein_set.add(tokens[2])
        test_data[tokens[0]][tokens[2]].add(tokens[3])

unit_tests = 0
for organism, organism_test_data in test_data.items():
    print(f"testing organism {organism}")
    start_time = time.time()
    ontology = EnrichmentOntology(f"{current_path}/sandbox/Data/ontology_{organism}.gz", lipid_parser = lipid_parser)
    print(f"time: {time.time() - start_time}\n")
    ontology.set_background(session, lipid_dict = lipid_dict, protein_set = protein_set, metabolite_set = metabolite_set)
    for molecule_name, term_ids in organism_test_data.items():
        for term_id in term_ids:
            if term_id not in session.search_terms:
                print("ERROR:", molecule_name, term_id)
                exit(-1)
            unit_tests += 1

print(f"Test passed ({unit_tests} unit tests) without errors")
