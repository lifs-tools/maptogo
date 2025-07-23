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
import traceback
import pickle
import os
from collections import deque, defaultdict
SKIP_LOADING = False

current_path = pathlib.Path(__file__).parent.resolve()
logger = logging.getLogger(__name__)
try:
    fisher_exact = ctypes.CDLL(f"{current_path}/assets/fisher_exact.so")
    fisher_exact.exact_fisher.restype = ctypes.c_double
except Exception as e:
    logger.error("Unable to load c++ library file(s)")



def time_elapsed(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        logger.info(f"Time elapsed for function '{func.__name__}': {end_time - start_time}s")
        return result
    return wrapper



class SessionEntry:
    def __init__(self):
        self.time = time.time()
        self.data = None
        self.results = []
        self.data_loaded = False
        self.search_terms = {}
        self.all_parent_nodes = {}
        self.num_background = 0
        self.use_bounded_fatty_acyls = True
        self.ontology = None
        self.domains = None
        self.background_lipids = None
        self.regulated_lipids = None
        self.background_proteins = None
        self.regulated_proteins = None
        self.background_metabolites = None
        self.regulated_metabolites = None
        self.background_transcripts = None
        self.regulated_transcripts = None
        self.background_list = []



class OntologyTerm:
    def __init__(self, _term_id, _name, _relations, _domain = None, _categories = None):
        if _domain == None: _domain = set()
        if _categories == None: _categories = set()

        self.term_id = sorted(list(_term_id) if type(_term_id) in {list, set} else list(_term_id.split("|")))
        self.name = _name
        self.relations = sorted([r for r in _relations if r != ""])
        self.domain = (set(_domain) if type(_domain) in {list, set} else set(_domain.split("|"))) - {"", "external"}
        self.categories = (set(_categories) if type(_categories) in {list, set} else set(_categories.split("|"))) - {""}

    def get_term_id(self, space = False):
        if space: " | ".join(sorted(list(self.term_id)))
        return "|".join(sorted(list(self.term_id)))



class OntologyResult:
    def __init__(
        self,
        _term,
        _pvalue,
        _source_terms,
        _fisher = None,
    ):
        if _fisher == None: _fisher = [0, 0, 0, 0]
        self.term = _term
        self.pvalue = _pvalue
        self.pvalue_corrected = _pvalue
        self.source_terms = _source_terms
        self.fisher_data = _fisher

        if min(_fisher) < 0:
            logger.error(f"ERROR: {self.term.name} / {self.term.get_term_id()} / {_fisher}")



class EnrichmentOntology:
    def __init__(self, file_name, ontology_name, lipid_parser):
        self.lipid_parser = lipid_parser
        self.ontology_terms = {}
        self.lipids = {}
        self.lipid_classes = {}
        self.carbon_chains = {}
        self.transcripts = {}
        self.proteins = {}
        self.clean_protein_ids = set()
        self.metabolites = {}
        self.clean_metabolite_ids = set()
        self.domains = set()
        self.metabolite_names = {}
        self.ontology_name = ontology_name
        self.molecule_lookup = {}
        self.reviewed_proteins = set()

        if SKIP_LOADING: return

        try:
            with gzip.open(file_name) as input_stream:
                input_content = [l for l in input_stream.read().decode("utf8").split("\n") if len(l) > 0]
                term_id = ""
                name = ""

                for line in input_content:
                    tokens = line.strip(" ").split("\t")
                    if len(tokens) < 5: continue
                    term_id = tokens[0]
                    name = tokens[1]
                    relations = tokens[2]
                    synonyms = tokens[3]
                    domain = tokens[4]
                    categories = set() if len(tokens) <= 5 else tokens[5]

                    is_lipid_species = False
                    is_lipid_class = False
                    is_carbon_chain = False
                    is_reviewed_protein = False
                    is_metabolite = False
                    is_stable_transcript = False
                    is_unreviewed_protein = False
                    is_any = False
                    relations = relations.split("|")
                    synonyms = synonyms.split("|")
                    for relation in relations:
                        if relation.startswith("LS:000000"):
                            match relation[9]:
                                case "1": is_lipid_class = True
                                case "2": is_lipid_species = True
                                case "3": is_carbon_chain = True
                                case "4": is_reviewed_protein = True
                                case "5": is_metabolite = True
                                case "6": is_stable_transcript = True
                                case "7": is_unreviewed_protein = True
                            is_any = is_lipid_class | is_lipid_species | is_carbon_chain | is_reviewed_protein | is_metabolite | is_stable_transcript | is_unreviewed_protein

                    term = OntologyTerm(term_id, name, relations, domain, categories)
                    str_term_id = term.get_term_id()
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

                        elif is_reviewed_protein:
                            if str_term_id not in self.proteins:
                                self.proteins[str_term_id] = term
                                self.reviewed_proteins.add(str_term_id)

                        elif is_metabolite:
                            if str_term_id not in self.metabolites: self.metabolites[str_term_id] = term

                        elif is_stable_transcript:
                            if str_term_id not in self.transcripts: self.transcripts[str_term_id] = term

                        elif is_unreviewed_protein:
                            if str_term_id not in self.proteins: self.proteins[str_term_id] = term

                    for d in domain.split("|"):
                        if len(d) > 0 and d != "external":
                            self.domains.add(d)

                    for t_id in term_id.split("|"):
                        self.ontology_terms[t_id] = term

        except Exception as e:
            logger.error("".join(traceback.format_tb(e.__traceback__)))
            logger.error(e)

        self.clean_protein_ids = set([key.replace("UNIPROT:", "") for key in self.proteins.keys()])
        self.clean_metabolite_ids = set([key.replace("CHEBI:", "") for key in self.metabolites.keys()])
        self.metabolite_names = {term.name: term for _, term in self.metabolites.items()}

        try:
            for line in open(f"{current_path}/Data/additional_links.csv").read().split("\n"):
                if len(line) < 2: continue
                tokens = line.split("\t")
                if len(tokens) < 2 or tokens[0] not in self.ontology_terms or tokens[1] not in self.ontology_terms: continue
                self.ontology_terms[tokens[0]].relations.append(tokens[1])

        except Exception as e:
            logger.error("".join(traceback.format_tb(e.__traceback__)))
            logger.error(e)

        # clean up ontology
        for term_id, term in self.ontology_terms.items():
            if not term.relations or type(term.relations[0]) == OntologyTerm: continue
            term.relations = sorted(
                list(
                    {self.ontology_terms[t] for t in term.relations if t in self.ontology_terms}
                ),
                key = lambda term: term.get_term_id()
            )

        logger.info(f"Loaded '{self.ontology_name}' with {len(self.ontology_terms)} terms.")


    @time_elapsed
    def set_background(self, session, lipid_dict = {}, protein_set = set(), metabolite_set = set(), transcript_set = {}):
        session.search_terms = defaultdict(list)
        session.all_parent_nodes = {}
        all_parent_nodes = session.all_parent_nodes
        session.num_background = len(lipid_dict) + len(protein_set) + len(metabolite_set) + len(transcript_set)
        all_paths = []

        def fill_path(lipid_term, lipid_name, lipid_input_name):
            path = []
            if lipid_term == None:
                if lipid_input_name != lipid_name: path.append(lipid_input_name)
                path.append(lipid_name)
            else:
                if lipid_input_name != lipid_term.name: path.append(lipid_input_name)
                path.append(lipid_term)
            return path

        for lipid_input_name, lipid in lipid_dict.items():
            if lipid is None: continue
            lipid.lipid.sort_fatty_acyl_chains() # only effects lipids on molecular species level or lower
            lipid_name = lipid.get_lipid_string()
            lipid_name_class = lipid.get_extended_class()
            parent_nodes = {}
            all_parent_nodes[lipid_input_name] = parent_nodes

            lipid_term = self.lipids[lipid_name] if lipid_name in self.lipids else None
            default_path = fill_path(lipid_term, lipid_name, lipid_input_name)
            if lipid_term != None:
                path = list(default_path)
                all_paths.append([lipid_input_name, tuple(path), parent_nodes])

            if lipid_name_class in self.lipid_classes:
                for term in self.lipid_classes[lipid_name_class]:
                    path = list(default_path)
                    path.append(term)
                    all_paths.append([lipid_input_name, tuple(path), parent_nodes])

            if session.use_bounded_fatty_acyls:
                for fa in lipid.lipid.fa_list:
                    if fa.num_carbon == 0: continue

                    if fa.lipid_FA_bond_type in {LipidFaBondType.LCB_REGULAR, LipidFaBondType.LCB_EXCEPTION}:
                        lcb_string = "L" + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                        if lcb_string in self.carbon_chains:
                            lcb_term = self.carbon_chains[lcb_string]
                            if lipid_term != lcb_term:
                                path = list(default_path)
                                path.append(lcb_term)
                                all_paths.append([lipid_input_name, tuple(path), parent_nodes])

                    else:
                        fa_string = "FA " + fa.to_string(lipid.lipid.info.level)
                        if fa_string in self.lipids:
                            fa_term = self.lipids[fa_string]
                            if lipid_term != fa_term:
                                path = list(default_path)
                                path.append(fa_term)
                                all_paths.append([lipid_input_name, tuple(path), parent_nodes])

                        fa_string = "C" + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                        if fa_string in self.carbon_chains:
                            fa_term = self.carbon_chains[fa_string]
                            if lipid_term != fa_term:
                                path = list(default_path)
                                path.append(fa_term)
                                all_paths.append([lipid_input_name, tuple(path), parent_nodes])

        for protein_input_name in protein_set:
            if protein_input_name not in self.proteins: continue
            parent_nodes = {}
            all_parent_nodes[protein_input_name] = parent_nodes
            all_paths.append([protein_input_name, (self.proteins[protein_input_name], ), parent_nodes])

        for metabolite_input_name in metabolite_set:
            if metabolite_input_name in self.metabolites:
                parent_nodes = {}
                all_parent_nodes[metabolite_input_name] = parent_nodes
                all_paths.append([metabolite_input_name, (self.metabolites[metabolite_input_name], ), parent_nodes])
                continue
            elif metabolite_input_name.lower() in self.metabolite_names:
                parent_nodes = {}
                all_parent_nodes[metabolite_input_name] = parent_nodes
                all_paths.append([metabolite_input_name, (self.metabolite_names[metabolite_input_name.lower()], ), parent_nodes])

        for transcript_input_name in transcript_set:
            transcript_name = transcript_input_name.split(".")[0]
            if transcript_name not in self.transcripts: continue
            parent_nodes = {}
            all_parent_nodes[transcript_input_name] = parent_nodes
            all_paths.append([transcript_input_name, (self.transcripts[transcript_name], ), parent_nodes])

        # run all registered molecules
        search_terms = session.search_terms
        for molecule_input_name, path, parent_nodes in all_paths:
            queue = deque([path[-1]])
            for i, p in enumerate(path): parent_nodes[p] = path[i - 1] if i > 0 else None
            while queue:
                term = queue.pop()
                if term.domain: search_terms[term].append(molecule_input_name)
                for relation_term in term.relations:
                    if relation_term not in parent_nodes:
                        parent_nodes[relation_term] = term
                        queue.append(relation_term)

        for term, term_molecules in search_terms.items():
            search_terms[term] = set(term_molecules)


    @time_elapsed
    def enrichment_analysis(self, session, target_set, enrichment_domains, term_regulation = "greater"):
        if len(target_set) == 0 or session.num_background < 2 or len(enrichment_domains) == 0: return []

        search_terms = session.search_terms
        enrichment_domains = set(enrichment_domains)
        result_list = [None] * len(search_terms)

        try: # C++ implementation, just way faster
            side = 0 if term_regulation == "two-sided" else (1 if term_regulation == "less" else 2)
            for i, (term, term_molecules) in enumerate(search_terms.items()):
                if not (term.domain & enrichment_domains): continue
                target_number = len(term_molecules & target_set)
                p_hyp = fisher_exact.exact_fisher(target_number, len(term_molecules), len(target_set), session.num_background, side)
                if p_hyp == 0: continue
                result_list[i] = OntologyResult(
                    term,
                    p_hyp,
                    term_molecules,
                    [
                        target_number,
                        len(term_molecules) - target_number,
                        len(target_set) - target_number,
                        session.num_background - len(term_molecules) - len(target_set) + target_number,
                    ],
                )

        except Exception as e:
            logger.error("".join(traceback.format_tb(e.__traceback__)))
            logger.error("C++ implementation of fisher exact test failed.")
            for i, (term, term_molecules) in enumerate(search_terms.items()):
                if not (term.domain & enrichment_domains): continue
                target_number = len(term_molecules & target_set)
                a, b, c, d = (
                    target_number,
                    len(term_molecules) - target_number,
                    len(target_set) - target_number,
                    session.num_background - len(term_molecules) - len(target_set) + target_number,
                )
                p_hyp = stats.fisher_exact([[a, b], [c, d]], alternative = term_regulation)[1]
                if p_hyp == 0: continue
                result_list[i] = OntologyResult(
                    term,
                    p_hyp,
                    term_molecules,
                    [a, b, c, d]
                )

        return [result for result in result_list if result != None]
