from EnrichmentDataStructure import EnrichmentOntology, current_path, SessionEntry
#from pygoslin.domain.LipidFaBondType import LipidFaBondType
#from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import LipidParser

lipid_parser = LipidParser()
session = SessionEntry()
ontology = EnrichmentOntology(f"{current_path}/Data/ontology_10090.gz", lipid_parser = lipid_parser)
ontology.set_background(session, lipid_dict = {"12-HETE": lipid_parser.parse("12-HETE")})

term_ids = ["GO:0006979", "GO:0001659", "GO:0004602", "GO:0005829", "GO:0002862", "SMP0063595", "SMP0120510"]

for term_id in term_ids:
    if term_id not in session.search_terms:
        print("ERROR:", term_id)
        exit(-1)
