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

WITH_LOOKUP = True
WITH_LOOKUP_STORAGE = True
SKIP_LOADING = False

def time_elapsed(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"Time elapsed for function '{func.__name__}': {end_time - start_time}s")
        return result
    return wrapper

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
        self.results = []
        self.data_loaded = False
        self.search_terms = {}
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


class TracebackGraph:
    def __init__(self, path):
        self.root = path[0]
        self.nodes = {}
        self.current_node = None
        for p in path:
            if p in self.nodes:
                logger.warning("Traceback: got path {path} with '{p}' at least twice inside")
                continue
            self.nodes[p] = self.current_node
            self.current_node = p
        self.current_depth = len(path)

    def add_node(self, node):
        if node in self.nodes: return
        self.nodes[node] = self.current_node
        self.current_node = node
        self.current_depth += 1

    def step_back(self):
        if self.current_node == None: return
        self.current_node = self.nodes[self.current_node]
        self.current_depth -= 1

    def get_path(self, node):
        if node not in self.nodes: return []
        path = []
        current_node = node
        while current_node != None:
            path.append(current_node)
            c = current_node
            current_node = self.nodes[current_node]
        return path[::-1]

    def go_back(self, path_len):
        if self.current_node == None or path_len <= 0: return
        while self.current_depth > path_len:
            self.current_node = self.nodes[self.current_node]
            self.current_depth -= 1


class OntologyTerm:
    def __init__(self, _term_id, _name, _relations, _domain = None, _categories = None):
        if _domain == None: _domain = set()
        if _categories == None: _categories = set()

        self.term_ontology_id = -1
        self.term_id = sorted(list(_term_id) if type(_term_id) in {list, set} else list(_term_id.split("|")))
        self.name = _name
        self.relations = sorted([r for r in _relations if r != ""])
        self.domain = (set(_domain) if type(_domain) in {list, set} else set(_domain.split("|"))) - {"", "external"}
        self.categories = (set(_categories) if type(_categories) in {list, set} else set(_categories.split("|"))) - {""}

    def to_string(self):
        return f"{'|'.join(self.term_id)}\t{self.name}\t{'|'.join(self.relations)}\t{'|'.join(self.domain)}\n"

    def get_term_id(self, space = False):
        if space: " | ".join(sorted(list(self.term_id)))
        return "|".join(sorted(list(self.term_id)))



class OntologyResult:
    def __init__(
        self,
        _term_id,
        _term,
        _number_background,
        _number_background_events,
        _targets,
        _pvalue,
        _source_terms,
        _fisher = None,
        _leaf = True,
    ):
        if _fisher == None: _fisher = [0, 0, 0, 0]
        self.term_id = _term_id
        self.term = _term
        self.number_background = _number_background
        self.number_background_events = _number_background_events
        self.targets = _targets
        self.pvalue = _pvalue
        self.pvalue_corrected = _pvalue
        self.source_terms = _source_terms
        self.fisher_data = _fisher
        self.leaf = _leaf

        if min(_fisher) < 0:
            print("ERROR: ", self.term.name, self.term.get_term_id(), _fisher)



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

        if SKIP_LOADING: return

        pickle_file_name = file_name.replace(".gz", ".pkl.gz")
        pickle_file_exists = os.path.exists(pickle_file_name)
        if WITH_LOOKUP_STORAGE and pickle_file_exists:
            try:
                with gzip.open(pickle_file_name) as input_stream:
                    load_obj = pickle.load(input_stream)
                    self.ontology_terms = load_obj["ontology_terms"]
                    self.lipids = load_obj["lipids"]
                    self.lipid_classes = load_obj["lipid_classes"]
                    self.carbon_chains = load_obj["carbon_chains"]
                    self.transcripts = load_obj["transcripts"]
                    self.proteins = load_obj["proteins"]
                    self.clean_protein_ids = load_obj["clean_protein_ids"]
                    self.metabolites = load_obj["metabolites"]
                    self.clean_metabolite_ids = load_obj["clean_metabolite_ids"]
                    self.domains = load_obj["domains"]
                    self.metabolite_names = load_obj["metabolite_names"]
                    self.ontology_name = load_obj["ontology_name"]
                    self.molecule_lookup = load_obj["molecule_lookup"]

                    logger.info(f"Dump loaded '{self.ontology_name}' with {len(self.ontology_terms)} terms and {len(self.molecule_lookup)} lookups.")
                    return

            except Exception as e:
                logger.error("".join(traceback.format_tb(e.__traceback__)))
                logger.error(e)

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
                    is_protein = False
                    is_metabolite = False
                    is_stable_transcript = False
                    is_stable_gene = False
                    is_stable_protein = False
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
                                case "6": is_stable_transcript = True
                                case "7": is_stable_gene = True
                                case "8": is_stable_protein = True
                            is_any = is_lipid_class | is_lipid_species | is_carbon_chain | is_protein | is_metabolite | is_stable_transcript | is_stable_gene | is_stable_protein

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

                        elif is_protein:
                            if str_term_id not in self.proteins: self.proteins[str_term_id] = term

                        elif is_metabolite:
                            if str_term_id not in self.metabolites: self.metabolites[str_term_id] = term

                        elif is_stable_transcript or is_stable_gene or is_stable_protein:
                            if str_term_id not in self.transcripts: self.transcripts[str_term_id] = term

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
        for term in self.ontology_terms.values():
            term.relations = [t for t in term.relations if t in self.ontology_terms]

        logger.info(f"Loaded '{self.ontology_name}' with {len(self.ontology_terms)} terms.")

        if WITH_LOOKUP:
            try:
                # creating lookup table for molecules
                for molecule_dict in [self.proteins, self.metabolites, self.transcripts]:
                    for molecule_input_name, input_term in molecule_dict.items():
                        path = (input_term, )
                        graph = TracebackGraph(path)
                        self.molecule_lookup[path] = [set(), graph]
                        m_lookup = self.molecule_lookup[path]

                        queue = [(input_term, len(path))]
                        visited_terms = {input_term}
                        while queue:
                            term, path_len = queue.pop()
                            graph.go_back(path_len)
                            graph.add_node(term)
                            contains_no_domain = len(term.domain) == 0

                            if not contains_no_domain: m_lookup[0].add(term)

                            for parent_term_id in term.relations:
                                parent_term = self.ontology_terms[parent_term_id]
                                if parent_term not in visited_terms:
                                    visited_terms.add(parent_term)
                                    queue.append((parent_term, graph.current_depth))
            except Exception as e:
                logger.error("".join(traceback.format_tb(e.__traceback__)))
                logger.error(e)


        if WITH_LOOKUP and WITH_LOOKUP_STORAGE and not pickle_file_exists:
            try:
                with gzip.open(pickle_file_name, 'wb') as output_stream:
                    dump_obj = {
                        "ontology_terms": self.ontology_terms,
                        "lipids": self.lipids,
                        "lipid_classes": self.lipid_classes,
                        "carbon_chains": self.carbon_chains,
                        "transcripts": self.transcripts,
                        "proteins": self.proteins,
                        "clean_protein_ids": self.clean_protein_ids,
                        "metabolites": self.metabolites,
                        "clean_metabolite_ids": self.clean_metabolite_ids,
                        "domains": self.domains,
                        "metabolite_names": self.metabolite_names,
                        "ontology_name": self.ontology_name,
                        "molecule_lookup": self.molecule_lookup,
                    }
                    pickle.dump(dump_obj, output_stream)
                    logger.info(f"Stored {file_name.replace(".gz", ".pkl.gz")} dump")

            except Exception as e:
                logger.error("".join(traceback.format_tb(e.__traceback__)))
                logger.error(e)


    @time_elapsed
    def set_background(self, session, lipid_dict = {}, protein_set = set(), metabolite_set = set(), transcript_set = {}):
        session.search_terms = {}
        session.num_background = len(lipid_dict) + len(protein_set) + len(metabolite_set) + len(transcript_set)
        all_paths = set()

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

            lipid_term = self.lipids[lipid_name] if lipid_name in self.lipids else None
            if lipid_term != None:
                path = fill_path(lipid_term, lipid_name, lipid_input_name)
                all_paths.add((lipid_input_name, tuple(path)))

            if lipid_name_class in self.lipid_classes:
                for term in self.lipid_classes[lipid_name_class]:
                    path = fill_path(lipid_term, lipid_name, lipid_input_name)
                    path.append(term)
                    all_paths.add((lipid_input_name, tuple(path)))

            if lipid != None and lipid.lipid != None and session.use_bounded_fatty_acyls:
                for fa in lipid.lipid.fa_list:
                    if fa.num_carbon == 0: continue

                    if fa.lipid_FA_bond_type in {LipidFaBondType.LCB_REGULAR, LipidFaBondType.LCB_EXCEPTION}:
                        lcb_string = "L" + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                        if lcb_string in self.carbon_chains:
                            lcb_term = self.carbon_chains[lcb_string]
                            if lipid_term != lcb_term:
                                path = fill_path(lipid_term, lipid_name, lipid_input_name)
                                path.append(lcb_term)
                                all_paths.add((lipid_input_name, tuple(path)))

                    else:
                        fa_string = "FA " + fa.to_string(lipid.lipid.info.level)
                        if fa_string in self.lipids:
                            fa_term = self.lipids[fa_string]
                            if lipid_term != fa_term:
                                path = fill_path(lipid_term, lipid_name, lipid_input_name)
                                path.append(fa_term)
                                all_paths.add((lipid_input_name, tuple(path)))

                        fa_string = "C" + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                        if fa_string in self.carbon_chains:
                            fa_term = self.carbon_chains[fa_string]
                            if lipid_term != fa_term:
                                path = fill_path(lipid_term, lipid_name, lipid_input_name)
                                path.append(fa_term)
                                all_paths.add((lipid_input_name, tuple(path)))

        for protein_input_name in protein_set:
            if protein_input_name not in self.proteins: continue
            term = self.proteins[protein_input_name]
            all_paths.add((protein_input_name, (term, )))

        for metabolite_input_name in metabolite_set:
            if metabolite_input_name in self.metabolites:
                term = self.metabolites[metabolite_input_name]

            elif metabolite_input_name.lower() in self.metabolite_names:
                term = self.metabolite_names[metabolite_input_name.lower()]

            else: continue

            all_paths.add((metabolite_input_name, (term, )))

        for transcript_input_name in transcript_set:
            transcript_name = transcript_input_name.split(".")[0]
            if transcript_name not in self.transcripts: continue
            term = self.transcripts[transcript_name]
            all_paths.add((transcript_input_name, (term, )))

        # run all registered molecules
        search_terms = session.search_terms
        for molecule_input_name, path in all_paths:
            if path not in self.molecule_lookup:
                graph = TracebackGraph(path)
            else:
                g = self.molecule_lookup[path][1]
                for term in self.molecule_lookup[path][0]:
                    if term not in search_terms: search_terms[term] = {molecule_input_name: g}
                    else: search_terms[term][molecule_input_name] = g
                continue

            queue = [(path[-1], len(path))]
            visited_terms = {path[-1]}
            while queue:
                term, path_len = queue.pop()
                graph.go_back(path_len)
                graph.add_node(term)
                contains_no_domain = len(term.domain) == 0

                if not contains_no_domain:
                    if term not in search_terms: search_terms[term] = {molecule_input_name: graph}
                    else: search_terms[term][molecule_input_name] = graph

                for parent_term_id in term.relations:
                    parent_term = self.ontology_terms[parent_term_id]
                    if parent_term not in visited_terms:
                        visited_terms.add(parent_term)
                        queue.append((parent_term, graph.current_depth))



    @time_elapsed
    def enrichment_analysis(self, session, target_set, enrichment_domains, term_regulation = "greater"):
        if len(target_set) == 0 or session.num_background == 0 or len(enrichment_domains) == 0: return []

        enrichment_domains = set(enrichment_domains)
        try: # C++ implementation, just way faster
            result_list = [None] * len(session.search_terms)
            side = 0 if term_regulation == "two-sided" else (1 if term_regulation == "less" else 2)
            visited_terms = set()
            for i, (term, term_metabolites) in enumerate(session.search_terms.items()):
                if (
                    len(term_metabolites) == 0
                    or len(term.domain & enrichment_domains) == 0
                    or term in visited_terms
                ): continue
                visited_terms.add(term)
                target_number = len(term_metabolites.keys() & target_set)
                p_hyp = fisher_exact.exact_fisher(target_number, len(term_metabolites), len(target_set), session.num_background, side)
                if p_hyp == 0: continue
                result_list[i] = OntologyResult(
                    term.get_term_id(),
                    term,
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
            logger.error("".join(traceback.format_tb(e.__traceback__)))
            logger.error("C++ implementation of fisher exact test failed.")
            result_list = [None] * len(session.search_terms)
            visited_terms = set()
            for i, (term, term_metabolites) in enumerate(session.search_terms.items()):
                if (
                    len(term_metabolites) == 0
                    or len(term.domain & enrichment_domains) == 0
                    or term in visited_terms
                ): continue
                visited_terms.add(term)
                target_number = len(term_metabolites.keys() & target_set)

                a, b, c, d = (
                    target_number,
                    len(term_metabolites) - target_number,
                    len(target_set) - target_number,
                    session.num_background - len(term_metabolites) - len(target_set) + target_number,
                )
                p_hyp = stats.fisher_exact([[a, b], [c, d]], alternative = term_regulation)[1]
                if p_hyp == 0: continue
                result_list[i] = OntologyResult(
                    term.get_term_id(),
                    term,
                    session.num_background,
                    len(term_metabolites),
                    len(target_set),
                    p_hyp,
                    term_metabolites,
                    [a, b, c, d]
                )

        return [result for result in result_list if result != None]
