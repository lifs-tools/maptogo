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
import csv
from enum import Enum

try:
    @profile
    def test_dummy():
        pass

except Exception as e:
    def profile(func):
        def wrapper(*args, **kwargs):
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            logger.info(f"Time elapsed for function '{func.__name__}': {end_time - start_time}s")
            return result
        return wrapper


class TermType(Enum):
    # Undefined / non-inferable lipid level
    LIPID_CLASS = 1                 # special handling
    LIPID_SPECIES = 2               # special handling
    CARBON_CHAIN = 3                # special handling
    UNSPECIFIC_LIPID = 4            # special handling
    REVIEWED_PROTEIN = 5            # special handling
    UNREVIEWED_PROTEIN = 6          # special handling
    ENSEMBLE_PROTEIN = 7            # special handling
    METABOLITE = 8                  # special handling
    GENERIC_REACTION = 9
    ENSEMBLE_TRANSCRIPT = 10        # special handling
    ENSEMBLE_GENE = 11              # special handling
    GENE = 12
    INPUT_TERM = 90
    UNCLASSIFIED_TERM = 99


SKIP_LOADING = False

current_path = pathlib.Path(__file__).parent.resolve()
logger = logging.getLogger(__name__)
try:
    fisher_exact = ctypes.CDLL(f"{current_path}/assets/fisher_exact.so")
    fisher_exact.exact_fisher.restype = ctypes.c_double
    exact_fisher = fisher_exact.exact_fisher
except Exception as e:
    logger.error("Unable to load c++ library file(s)")


CHEBI_synonym_table_filename = f"{current_path}/Data/CHEBI_synonyms.csv.gz"
CHEBI_synonym_table = {}


if os.path.isfile(CHEBI_synonym_table_filename):
    logger.info("Reading ChEBI synonyms table")
    with gzip.open(CHEBI_synonym_table_filename, "rt") as input_stream:
        CHEBI_synonym_table = {tokens[0].lower(): tokens[1] for line in input_stream.read().split("\n") if len(line) > 0 and (tokens := line.split("\t")) and len(tokens) >= 2}
else:
    logger.warning("No ChEBI synonyms table found")


class SessionEntry:
    def __init__(self, cip = "0.0.0.0"):
        self.time = time.time()
        self.data = None
        self.results = []
        self.data_loaded = False
        self.search_terms = {}
        self.all_parent_nodes = {}
        self.num_background = 0
        self.use_bounded_fatty_acyls = False
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
        self.min_pvalue = "0.0001"
        self.max_pvalue = "0.05"
        self.ui = {}
        self.cip = cip



class OntologyTerm:
    def __init__(self, _term_id, _name, _relations, _domain = None, _categories = None):
        if _domain == None: _domain = set()
        if _categories == None: _categories = set()


        self.term_id = sorted(list(_term_id) if type(_term_id) in {list, set} else list(_term_id.split("|")))
        self.term_id_str = "|".join(sorted(list(self.term_id)))
        self.term_type = TermType.UNCLASSIFIED_TERM
        self.name = _name
        self.relations = sorted([r for r in _relations if r != ""])
        self.domain = (set(_domain) if type(_domain) in {list, set} else set(_domain.split("|"))) - {"", "external"}
        self.categories = (set(_categories) if type(_categories) in {list, set} else set(_categories.split("|"))) - {""}

    def get_term_id(self, space = False):
        return " | ".join(sorted(list(self.term_id)))






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
        self.lor = _fisher[0] * _fisher[3] / (_fisher[1] * _fisher[2]) if _fisher[1] > 0 and _fisher[2] > 0 else 0
        self.source_terms = _source_terms
        self.fisher_data = _fisher

        if min(_fisher) < 0:
            logger.error(f"ERROR: {self.term.name} / {self.term.term_id_str} / {_fisher}")




class EnrichmentOntology:
    @profile
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
        self.reviewed_proteins = set()
        self.unspecific_lipids = {}

        if SKIP_LOADING: return

        try:
            with gzip.open(file_name, "rt") as file_stream:
                lines = file_stream.read().split("\n")
                tokens_list = [l.split("\t") for l in lines if len(l) > 0]
                for tokens in tokens_list:
                    term_id, name, relations, synonyms, domain, categories = tokens
                    relations_split = relations.split("|")
                    synonyms = synonyms.split("|")
                    moea_pos = relations.find("MOEA:00000")
                    term_type = TermType(int(relations[moea_pos + 5 : moea_pos + 12])) if moea_pos > -1 else TermType.UNCLASSIFIED_TERM

                    term = OntologyTerm(term_id, name, relations_split, domain, categories)
                    term.term_type = term_type
                    str_term_id = term.term_id_str
                    match term_type:
                        case TermType.LIPID_CLASS: # lipid class
                            if name not in self.lipid_classes: self.lipid_classes[name] = []
                            self.lipid_classes[name].append(term)
                            for synonym in synonyms:
                                if synonym not in self.lipid_classes: self.lipid_classes[synonym] = []
                                self.lipid_classes[synonym].append(term)

                        case TermType.LIPID_SPECIES: # lipid species
                            if name not in self.lipids: self.lipids[name] = term

                        case TermType.CARBON_CHAIN: # carbon chain
                            if name not in self.carbon_chains: self.carbon_chains[name] = term
                            for synonym in synonyms:
                                self.carbon_chains[synonym] = term

                        case TermType.UNSPECIFIC_LIPID:
                            if name.lower() not in self.unspecific_lipids: self.unspecific_lipids[name.lower()] = term
                            for synonym in synonyms:
                                if synonym.lower() not in self.unspecific_lipids:
                                    self.unspecific_lipids[synonym.lower()] = term

                        case TermType.REVIEWED_PROTEIN: # reviewed protein
                            if str_term_id not in self.proteins:
                                self.proteins[str_term_id] = term
                                self.reviewed_proteins.add(str_term_id)

                        case TermType.UNREVIEWED_PROTEIN: # unreviewed protein
                            if str_term_id not in self.proteins: self.proteins[str_term_id] = term

                        case TermType.ENSEMBLE_PROTEIN | TermType.ENSEMBLE_TRANSCRIPT | TermType.ENSEMBLE_GENE: # ensemble
                            if str_term_id not in self.transcripts: self.transcripts[str_term_id] = term

                        case TermType.METABOLITE: # metabolite
                            if str_term_id not in self.metabolites: self.metabolites[str_term_id] = term


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
        self.metabolite_names = {term.name.lower(): term for _, term in self.metabolites.items()}
        for synonym, term_id in CHEBI_synonym_table.items():
            if term_id in self.metabolites:
                self.metabolite_names[synonym] = self.metabolites[term_id]

        # add additional knowlegde graph edges
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
                {self.ontology_terms[t] for t in term.relations if t in self.ontology_terms},
                key = lambda term: term.term_id_str
            )
        logger.info(f"Loaded '{self.ontology_name}' with {len(self.ontology_terms)} terms.")


    @profile
    def set_background(self, session, lipid_dict = {}, protein_set = set(), metabolite_set = set(), transcript_set = {}):
        session.search_terms = defaultdict(list)
        session.all_parent_nodes = {}
        all_parent_nodes = session.all_parent_nodes
        session.num_background = len(lipid_dict) + len(protein_set) + len(metabolite_set) + len(transcript_set)
        all_paths = []

        for lipid_input_name, lipid in lipid_dict.items():
            if lipid is None: continue
            elif type(lipid) == OntologyTerm:
                parent_nodes = {lipid: lipid_input_name, lipid_input_name: None}
                all_paths.append([lipid_input_name, lipid, parent_nodes])
                all_parent_nodes[lipid_input_name] = parent_nodes
                continue

            lipid.lipid.sort_fatty_acyl_chains() # only effects lipids on molecular species level or lower
            lipid_name = lipid.get_lipid_string()
            lipid_name_class = lipid.get_lipid_string(LipidLevel.CLASS)
            #lipid_name_class = lipid.get_extended_class()
            parent_nodes = {}
            all_parent_nodes[lipid_input_name] = parent_nodes
            lipid_term = self.lipids[lipid_name] if lipid_name in self.lipids else None

            start_term = None
            if not lipid_term:
                parent_nodes[lipid_input_name] = None
                if lipid_input_name != lipid_name:
                    parent_nodes[lipid_name] = lipid_input_name
                    start_term = lipid_name
                else:
                    start_term = lipid_input_name
            else:
                if lipid_input_name != lipid_term.name:
                    parent_nodes[lipid_input_name] = None
                    parent_nodes[lipid_term] = lipid_input_name
                else:
                    parent_nodes[lipid_term] = None
                start_term = lipid_term

            if lipid_term != None:
                all_paths.append([lipid_input_name, start_term, parent_nodes])

            if lipid_name_class in self.lipid_classes:
                for class_term in self.lipid_classes[lipid_name_class]:
                    parent_nodes[class_term] = start_term
                    all_paths.append([lipid_input_name, class_term, parent_nodes])

            if session.use_bounded_fatty_acyls:
                for fa in lipid.lipid.fa_list:
                    if fa.num_carbon == 0: continue

                    if fa.lipid_FA_bond_type in {LipidFaBondType.LCB_REGULAR, LipidFaBondType.LCB_EXCEPTION}:
                        lcb_string = "L" + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                        if lcb_string in self.carbon_chains:
                            lcb_term = self.carbon_chains[lcb_string]
                            if lipid_term != lcb_term:
                                parent_nodes[lipid_term] = start_term
                                all_paths.append([lipid_input_name, lcb_term, parent_nodes])
                        continue

                    fa_string = "FA " + fa.to_string(lipid.lipid.info.level)
                    if fa_string in self.lipids:
                        fa_term = self.lipids[fa_string]
                        if lipid_term != fa_term:
                            parent_nodes[fa_term] = start_term
                            all_paths.append([lipid_input_name, fa_term, parent_nodes])

                    fa_string = "C" + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                    if fa_string in self.carbon_chains:
                        c_term = self.carbon_chains[fa_string]
                        if lipid_term != c_term:
                            parent_nodes[c_term] = start_term
                            all_paths.append([lipid_input_name, c_term, parent_nodes])

        for protein_input_name in protein_set:
            if protein_input_name not in self.proteins: continue
            protein_term = self.proteins[protein_input_name]
            parent_nodes = {protein_term: None}
            all_parent_nodes[protein_input_name] = parent_nodes
            all_paths.append([protein_input_name, protein_term, parent_nodes])

        for metabolite_input_name in metabolite_set:
            if metabolite_input_name in self.metabolites:
                metabolite_term = self.metabolites[metabolite_input_name]
                parent_nodes = {metabolite_term: None}
                all_parent_nodes[metabolite_input_name] = parent_nodes
                all_paths.append([metabolite_input_name, metabolite_term, parent_nodes])
                continue

            elif metabolite_input_name.lower() in self.metabolite_names:
                metabolite_term = self.metabolite_names[metabolite_input_name.lower()]
                parent_nodes = {metabolite_term: None}
                all_parent_nodes[metabolite_input_name] = parent_nodes
                all_paths.append([metabolite_input_name, metabolite_term, parent_nodes])

        for transcript_input_name in transcript_set:
            transcript_name = transcript_input_name.split(".")[0]
            if transcript_name not in self.transcripts: continue
            transcript_term = self.transcripts[transcript_name]
            parent_nodes = {transcript_term: None}
            all_parent_nodes[transcript_input_name] = parent_nodes
            all_paths.append([transcript_input_name, transcript_term, parent_nodes])

        # run all registered molecules
        search_terms = session.search_terms
        for molecule_input_name, start_term, parent_nodes in all_paths:
            queue = [start_term]
            while queue:
                term = queue.pop()
                if term.domain: search_terms[term].append(molecule_input_name)
                for relation_term in term.relations:
                    if relation_term not in parent_nodes:
                        queue.append(relation_term)
                        parent_nodes[relation_term] = term

        for term, term_molecules in search_terms.items():
            search_terms[term] = set(term_molecules)


    @profile
    def enrichment_analysis(self, session, target_set, enrichment_domains, term_regulation = "greater"):
        if len(target_set) == 0 or session.num_background < 2 or len(enrichment_domains) == 0: return []

        search_terms, num_background = session.search_terms, session.num_background
        enrichment_domains = set(enrichment_domains)
        result_list = [None] * len(search_terms)
        len_target_set = len(target_set)

        try: # C++ implementation, just way faster
            side = 2 if term_regulation == "greater" else (1 if term_regulation == "less" else 0)
            for i, (term, term_molecules) in enumerate(search_terms.items()):
                if term.domain.isdisjoint(enrichment_domains): continue
                target_number = len(term_molecules & target_set)

                p_hyp = exact_fisher(target_number, len(term_molecules), len_target_set, num_background, side)
                if p_hyp == 0: continue
                result_list[i] = OntologyResult(
                    term,
                    p_hyp,
                    term_molecules,
                    [
                        target_number,
                        len(term_molecules) - target_number,
                        len_target_set - target_number,
                        num_background - len(term_molecules) - len_target_set + target_number,
                    ],
                )

        except Exception as e:
            logger.error("".join(traceback.format_tb(e.__traceback__)))
            logger.error("C++ implementation of fisher exact test failed.")
            for i, (term, term_molecules) in enumerate(search_terms.items()):
                if term.domain.isdisjoint(enrichment_domains): continue
                target_number = len(term_molecules & target_set)
                a, b, c, d = (
                    target_number,
                    len(term_molecules) - target_number,
                    len_target_set - target_number,
                    num_background - len(term_molecules) - len_target_set + target_number,
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
