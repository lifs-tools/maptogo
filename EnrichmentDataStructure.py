import gzip
from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidFaBondType import LipidFaBondType
import scipy.stats as stats
import logging
import json
import ctypes
import pathlib
import time

current_path = pathlib.Path(__file__).parent.resolve()

logger = logging.getLogger(__name__)
try:
    fisher_exact = ctypes.CDLL(f"{current_path}/assets/fisher_exact.so")
    fisher_exact.exact_fisher.restype = ctypes.c_double
except Exception as e:
    logger.error("Unable to load c++ library file(s)")



class SessionEntry:
    def __init__(self):
        self.time = time.time()
        self.data = None
        self.data_loaded = False
        self.search_terms = {}
        self.num_background = 0



class OntologyTerm:
    def __init__(self, _term_id, _name, _relations, _domain = ""):
        self.term_id = _term_id
        self.name = _name
        self.relations = list(_relations)
        self.domain = _domain



class OntologyResult:
    def __init__(
        self, _term, _number_background, _number_background_events, _targets, _pvalue, _source_terms, fisher = [0, 0, 0, 0]
    ):
        self.term = _term
        self.number_background = _number_background
        self.number_background_events = _number_background_events
        self.targets = _targets
        self.pvalue = _pvalue
        self.source_terms = _source_terms
        self.fisher_data = fisher

        if min(fisher) < 0:
            print("ERROR: ", self.term.name, self.term.term_id, fisher)



class EnrichmentOntology:
    def __init__(self, file_name, lipid_parser):
        self.ontology_terms = {}
        self.lipids = {}
        self.lipid_classes = {}
        self.carbon_chains = {}
        self.proteins = {}
        self.clean_protein_ids = set()
        self.metabolites = {}
        self.clean_metabolite_ids = set()
        self.domains = set()
        self.lipid_parser = lipid_parser

        try:
            with gzip.open(file_name) as input_stream:
                term_id = ""
                name = ""
                is_lipid_species = False
                is_lipid_class = False
                is_carbon_chain = False
                is_protein = False
                is_metabolite = False
                is_any = False
                domain, relations, synonyms = "", set(), []

                for line in input_stream:
                    line = line.decode("utf8").strip().strip(" ")

                    if line == "[Term]":
                        if term_id != "" and name != "": # and name.find("inding") == -1:
                            if term_id in self.ontology_terms:
                                raise Exception(
                                    "Term id '%s' in ontology already defined."
                                    % term_id
                                )

                            term = OntologyTerm(term_id, name, relations, domain)
                            if is_any:
                                if is_lipid_species:
                                    if name not in self.lipids: self.lipids[name] = term

                                elif is_lipid_class:
                                    if name not in self.lipid_classes: self.lipid_classes[name] = []
                                    self.lipid_classes[name].append(term)
                                    for synonym in synonyms:
                                        if synonym not in self.lipid_classes: self.lipid_classes[synonym] = []
                                        self.lipid_classes[synonym].append(term)

                                elif is_carbon_chain:
                                    if name not in self.carbon_chains: self.carbon_chains[name] = term
                                    for synonym in synonyms:
                                        self.carbon_chains[synonym] = term

                                elif is_protein:
                                    if name not in self.proteins: self.proteins[term_id] = term

                                elif is_metabolite:
                                    if name not in self.metabolites: self.metabolites[term_id] = term

                            if len(domain) > 0:
                                self.domains.add(domain)

                            self.ontology_terms[term_id] = term

                        term_id = ""
                        name = ""
                        is_lipid_species = False
                        is_lipid_class = False
                        is_carbon_chain = False
                        is_protein = False
                        is_metabolite = False
                        is_any = False
                        domain, relations, synonyms = "", set(), []

                    elif line[:4] == "id: ":
                        term_id = line[4:]

                    elif line[:6] == "name: ":
                        name = line[6:]

                    elif line[:6] == "is_a: ":
                        relation = line[6:].split(" ! ")[0]
                        relations.add(relation)
                        if relation[:3] == "LS:":
                            if relation == "LS:0000001": is_lipid_class = True
                            elif relation == "LS:0000002": is_lipid_species = True
                            elif relation == "LS:0000003": is_carbon_chain = True
                            elif relation == "LS:0000004": is_protein = True
                            elif relation == "LS:0000005": is_metabolite = True
                            is_any = True

                    elif line[:14] == "relationship: ":
                        line = line[14:]
                        pos = line.find(" ")
                        if pos > -1:
                            relations.add(line[pos + 1 :].split(" ! ")[0])

                    elif line[:11] == "namespace: ":
                        domain = line[11:]

                    elif line[:8] == "synonym:":
                        tokens = line.split('"')
                        if len(tokens) >= 2:
                            synonyms.append(tokens[1])

                    elif line == "is_obsolete: true":
                        term_id, name = "", ""

                if term_id != "" and name != "" and name.find("inding") == -1:
                    if term_id in self.ontology_terms:
                        raise Exception(
                            "Term id '%s' in ontology already defined." % term_id
                        )

                    term = OntologyTerm(term_id, name, relations, domain)

                    if is_lipid_species:
                        if name not in self.lipids: self.lipids[name] = term

                    elif is_lipid_class:
                        if name not in self.lipid_classes: self.lipid_classes[name] = term
                        for synonym in synonyms:
                            if synonym not in self.lipid_classes: self.lipid_classes[synonym] = []
                            self.lipid_classes[synonym].append(term)

                    elif is_carbon_chain:
                        if name not in self.carbon_chains: self.carbon_chains[name] = term
                        for synonym in synonyms:
                            if synonym not in self.carbon_chains: self.carbon_chains[synonym] = []
                            self.carbon_chains[synonym] = term

                    elif is_protein:
                        if name not in self.proteins: self.proteins[name] = term

                    elif is_metabolite:
                        if name not in self.metabolites: self.metabolites[name] = term

                    if len(domain) > 0:
                        self.domains.add(domain)

                    self.ontology_terms[term_id] = term

                if "external" in self.domains:
                    self.domains.remove("external")

        except Exception as e:
            logger.error(e)

        self.clean_protein_ids = set([key.replace("UNIPROT:", "") for key in self.proteins.keys()])
        self.clean_metabolite_ids = set([key.replace("CHEBI:", "") for key in self.metabolites.keys()])



    def recursive_event_adding(self, session, molecule_input_name, term_id, visited_terms, source_term = None, path = None):
        if term_id not in self.ontology_terms: return

        if source_term == None:
            source_term = term_id
            path = [term_id]
        term = self.ontology_terms[term_id]

        if term.domain not in self.domains:
            for parent_term_id in term.relations:
                if parent_term_id not in visited_terms:
                    visited_terms.add(parent_term_id)
                    path.append(parent_term_id)
                    self.recursive_event_adding(
                        session,
                        molecule_input_name,
                        parent_term_id,
                        visited_terms,
                        source_term,
                        path,
                    )
                    path.pop()

        if source_term != term_id:
            if term_id not in session.search_terms: session.search_terms[term_id] = {}
            session.search_terms[term_id][molecule_input_name] = list(path)



    def set_background(self, session, lipid_dict = {}, protein_set = set(), metabolite_set = set()):
        session.search_terms = {}
        session.num_background = len(lipid_dict) + len(protein_set) + len(metabolite_set)

        for lipid_input_name, lipid in lipid_dict.items():
            lipid.lipid.sort_fatty_acyl_chains() # only effects lipids on molecular species level or lower
            lipid_name = lipid.get_lipid_string()
            lipid_name_class = lipid.get_extended_class()

            visited_terms = set()
            if lipid_name in self.lipids:
                term = self.lipids[lipid_name]
                visited_terms.add(term.term_id)
                self.recursive_event_adding(
                    session, lipid_input_name, term.term_id, visited_terms
                )

            if lipid_name_class in self.lipid_classes:
                for term in self.lipid_classes[lipid_name_class]:
                    if term.term_id not in visited_terms:
                        visited_terms.add(term.term_id)
                        self.recursive_event_adding(
                            session, lipid_input_name, term.term_id, visited_terms
                        )

            if lipid != None and lipid.lipid != None:
                for fa in lipid.lipid.fa_list:
                    if fa.num_carbon == 0: continue

                    fa_string = (
                        "L"
                        if fa.lipid_FA_bond_type
                        in {LipidFaBondType.LCB_REGULAR, LipidFaBondType.LCB_EXCEPTION}
                        else "C"
                    ) + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                    if fa_string in self.carbon_chains:
                        visited_terms.add(self.carbon_chains[fa_string].term_id)
                        self.recursive_event_adding(
                            session,
                            lipid_input_name,
                            self.carbon_chains[fa_string].term_id,
                            visited_terms,
                        )

        for protein_input_name in protein_set:
            visited_terms = set()
            if protein_input_name not in self.proteins: continue
            term = self.proteins[protein_input_name]
            visited_terms.add(term.term_id)
            self.recursive_event_adding(
                session, protein_input_name, term.term_id, visited_terms
            )

        for metabolite_input_name in metabolite_set:
            visited_terms = set()
            if metabolite_input_name not in self.metabolites: continue
            term = self.metabolites[metabolite_input_name]
            visited_terms.add(term.term_id)
            self.recursive_event_adding(
                session, metabolite_input_name, term.term_id, visited_terms
            )



    def enrichment_analysis(self, session, target_set, enrichment_domains, term_regulation = "two-sided"):
        if session.num_background == 0 or len(enrichment_domains) == 0: return []

        try: # C++ implementation, just way faster
            result_list = [None] * len(session.search_terms)
            side = 0 if term_regulation == "two-sided" else (1 if term_regulation == "less" else 2)
            for i, (term_id, term_metabolites) in enumerate(session.search_terms.items()):
                if (
                    len(term_metabolites) == 0
                    or term_id not in self.ontology_terms
                    or self.ontology_terms[term_id].domain not in enrichment_domains
                ): continue

                target_number = len(term_metabolites.keys() & target_set)
                p_hyp = fisher_exact.exact_fisher(target_number, len(term_metabolites), len(target_set), session.num_background, side)
                result_list[i] = OntologyResult(
                    self.ontology_terms[term_id],
                    session.num_background,
                    len(term_metabolites),
                    len(target_set),
                    p_hyp,
                    term_metabolites,
                    [
                        target_number,
                        len(term_metabolites) - target_number,
                        len(target_set) - target_number,
                        session.num_background - len(term_metabolites) - len(target_set) + target_number,
                    ],
                )

        except Exception as e:
            logger.error("C++ implementation of fisher exact test failed.")
            result_list = [None] * len(session.search_terms)
            for i, (term_id, term_metabolites) in enumerate(session.search_terms.items()):
                if (
                    len(term_metabolites) == 0
                    or term_id not in self.ontology_terms
                    or self.ontology_terms[term_id].domain not in enrichment_domains
                ): continue

                target_number = len(term_metabolites.keys() & target_set)

                a, b, c, d = (
                    target_number,
                    len(term_metabolites) - target_number,
                    len(target_set) - target_number,
                    session.num_background - len(term_metabolites) - len(target_set) + target_number,
                )
                p_hyp = stats.fisher_exact([[a, b], [c, d]], alternative = term_regulation)[1]
                result_list[i] = OntologyResult(
                    self.ontology_terms[term_id],
                    session.num_background,
                    len(term_metabolites),
                    len(target_set),
                    p_hyp,
                    term_metabolites,
                    [a, b, c, d]
                )

        return [result for result in result_list if result != None]
