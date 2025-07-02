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
from collections import defaultdict

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
        self.use_bounded_fatty_acyls = True
        self.ontology = None



class OntologyTerm:
    def __init__(self, _term_id, _name, _relations, _domain = "", _xrefs = None):
        self.term_id = _term_id
        self.name = _name
        self.relations = list(_relations)
        self.domain = _domain
        self.xrefs = list(_xrefs) if _xrefs != None else []

    def to_string(self):
        return f"{self.term_id}\t{self.name}\t{'|'.join(self.relations)}\t{self.domain}\t{'|'.join(self.xrefs)}\n"


class OntologyResult:
    def __init__(
        self, _term, _number_background, _number_background_events, _targets, _pvalue, _source_terms, fisher = [0, 0, 0, 0]
    ):
        self.term = _term
        self.number_background = _number_background
        self.number_background_events = _number_background_events
        self.targets = _targets
        self.pvalue = _pvalue
        self.pvalue_corrected = _pvalue
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
        self.metabolite_names = {}

        try:
            with gzip.open(file_name) as input_stream:
                input_content = [l for l in input_stream.read().decode("utf8").split("\n") if len(l) > 0]
                term_id = ""
                name = ""

                for line in input_content:
                    tokens = line.strip(" ").split("\t")
                    if len(tokens) < 6: continue
                    term_id, name, relations, synonyms, domain, xrefs = tokens
                    xrefs = [xref for xref in xrefs.split("|") if xref in self.ontology_terms and self.ontology_terms[xref].domain == domain]

                    is_lipid_species = False
                    is_lipid_class = False
                    is_carbon_chain = False
                    is_protein = False
                    is_metabolite = False
                    is_any = False
                    relations = relations.split("|")
                    synonyms = synonyms.split("|")
                    for relation in relations:
                        if relation[:9] == "LS:000000":
                            match relation[9]:
                                case "1": is_lipid_class = True
                                case "2": is_lipid_species = True
                                case "3": is_carbon_chain = True
                                case "4": is_protein = True
                                case "5": is_metabolite = True
                            is_any = is_lipid_class | is_lipid_species | is_carbon_chain | is_protein | is_metabolite

                    term = OntologyTerm(term_id, name, relations, domain, xrefs)
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
                            if term_id not in self.proteins: self.proteins[term_id] = term

                        elif is_metabolite:
                            if term_id not in self.metabolites: self.metabolites[term_id] = term

                    if len(domain) > 0 and domain != "external":
                        self.domains.add(domain)

                    self.ontology_terms[term_id] = term

        except Exception as e:
            logger.error(e)

        self.clean_protein_ids = set([key.replace("UNIPROT:", "") for key in self.proteins.keys()])
        self.clean_metabolite_ids = set([key.replace("CHEBI:", "") for key in self.metabolites.keys()])
        self.metabolite_names = {term.name: term for _, term in self.metabolites.items()}

        try:
            for line in open("Data/additional_links.csv").read().split("\n"):
                if len(line) < 2: continue
                tokens = line.split("\t")
                if len(tokens) < 2 or tokens[0] not in self.ontology_terms or tokens[1] not in self.ontology_terms: continue
                self.ontology_terms[tokens[0]].relations.append(tokens[1])

        except Exception as e:
            logger.error(e)
            pass


    def recursive_event_adding(self, session_and_molecule_input_name, visited_terms, path):
        term_id = path[-1]
        term = self.ontology_terms[term_id]
        accepting_state = term.domain != ""

        if term.domain not in self.domains or len(term.xrefs) > 0:
            for is_xref, terms in zip([False, True], [term.relations, term.xrefs]):
                for parent_term_id in terms:
                    if is_xref and term.domain != "": accepting_state = False
                    if parent_term_id not in visited_terms and parent_term_id in self.ontology_terms:
                        visited_terms.add(parent_term_id)
                        path.append(parent_term_id)
                        self.recursive_event_adding(
                            session_and_molecule_input_name, visited_terms, path,
                        )
                        path.pop()

        if accepting_state and path[0] != term_id:
            session, molecule_input_name = session_and_molecule_input_name
            if term_id not in session.search_terms: session.search_terms[term_id] = {molecule_input_name: list(path)}
            else: session.search_terms[term_id][molecule_input_name] = list(path)


    def set_background(self, session, lipid_dict = {}, protein_set = set(), metabolite_set = set()):
        session.search_terms = {}
        session.num_background = len(lipid_dict) + len(protein_set) + len(metabolite_set)

        def fill_path(lipid_term, lipid_name, lipid_input_name):
            path = []
            if lipid_term == None:
                if lipid_input_name != lipid_name: path.append(lipid_input_name)
                path.append(lipid_name)
            else:
                if lipid_input_name != lipid_term.name: path.append(lipid_input_name)
                path.append(lipid_term.term_id)
            return path

        for lipid_input_name, lipid in lipid_dict.items():
            lipid.lipid.sort_fatty_acyl_chains() # only effects lipids on molecular species level or lower
            lipid_name = lipid.get_lipid_string()
            lipid_name_class = lipid.get_extended_class()

            visited_terms = set()
            lipid_term = self.lipids[lipid_name] if lipid_name in self.lipids else None
            if lipid_term != None:
                path = fill_path(lipid_term, lipid_name, lipid_input_name)
                visited_terms.add(lipid_term.term_id)
                self.recursive_event_adding(
                    (session, lipid_input_name), visited_terms, path,
                )

            if lipid_name_class in self.lipid_classes:
                for term in self.lipid_classes[lipid_name_class]:
                    if term.term_id not in visited_terms:
                        visited_terms.add(term.term_id)
                        path = fill_path(lipid_term, lipid_name, lipid_input_name)
                        path.append(term.term_id)
                        self.recursive_event_adding(
                            (session, lipid_input_name), visited_terms, path,
                        )

            if lipid != None and lipid.lipid != None and session.use_bounded_fatty_acyls:
                for fa in lipid.lipid.fa_list:
                    if fa.num_carbon == 0: continue

                    if fa.lipid_FA_bond_type in {LipidFaBondType.LCB_REGULAR, LipidFaBondType.LCB_EXCEPTION}:
                        lcb_string = "L" + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                        if lcb_string in self.carbon_chains:
                            lcb_term_id = self.carbon_chains[lcb_string].term_id
                            visited_terms.add(lcb_term_id)
                            path = fill_path(lipid_term, lipid_name, lipid_input_name)
                            path.append(lcb_term_id)
                            self.recursive_event_adding(
                                (session, lipid_input_name), visited_terms, path,
                            )

                    else:
                        fa_string = "FA " + fa.to_string(lipid.lipid.info.level)
                        if fa_string in self.lipids:
                            fa_term = self.lipids[fa_string]
                            fa_term_id = fa_term.term_id
                            visited_terms.add(fa_term_id)
                            path = fill_path(lipid_term, lipid_name, lipid_input_name)
                            path.append(fa_term_id)
                            self.recursive_event_adding(
                                (session, lipid_input_name), visited_terms, path,
                            )

                        fa_string = "C" + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                        if fa_string in self.carbon_chains:
                            fa_term_id = self.carbon_chains[fa_string].term_id
                            visited_terms.add(fa_term_id)
                            path = fill_path(lipid_term, lipid_name, lipid_input_name)
                            path.append(fa_term_id)
                            self.recursive_event_adding(
                                (session, lipid_input_name), visited_terms, path,
                            )

        for protein_input_name in protein_set:
            visited_terms = set()
            if protein_input_name not in self.proteins: continue
            term = self.proteins[protein_input_name]
            visited_terms.add(term.term_id)
            self.recursive_event_adding(
                (session, protein_input_name), visited_terms, [term.term_id]
            )

        for metabolite_input_name in metabolite_set:
            visited_terms = set()
            if metabolite_input_name in self.metabolites:
                term = self.metabolites[metabolite_input_name]

            elif metabolite_input_name.lower() in self.metabolite_names:
                term = self.metabolite_names[metabolite_input_name.lower()]

            else: continue

            visited_terms.add(term.term_id)
            self.recursive_event_adding(
                (session, metabolite_input_name), visited_terms, [term.term_id]
            )


    def enrichment_analysis(self, session, target_set, enrichment_domains, term_regulation = "two-sided"):
        session.result = []
        if session.num_background == 0 or len(enrichment_domains) == 0: return

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

        session.result = [result for result in result_list if result != None]
