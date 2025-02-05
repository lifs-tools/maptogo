import gzip
from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidFaBondType import LipidFaBondType
import scipy.stats as stats
import logging

logger = logging.getLogger(__name__)


class OntologyTerm:
    def __init__(self, _term_id, _name, _relations, _domain=""):
        self.term_id = _term_id
        self.name = _name
        self.relations = list(_relations)
        self.domain = _domain
        self.lipid_paths = {}


class OntologyResult:
    def __init__(
        self, _term, _number_background, _number_background_events, _targets, _pvalue
    ):
        self.term = _term
        self.number_background = _number_background
        self.number_background_events = _number_background_events
        self.targets = _targets
        self.pvalue = _pvalue


class EnrichmentOntology:
    def __init__(self, file_name, lipid_parser):
        self.ontology_terms = {}
        self.lipids = {}
        self.lipid_classes = {}
        self.carbon_chains = {}
        self.domains = set()
        self.search_terms = {}
        self.lipid_parser = lipid_parser
        self.num_background = 0

        try:
            with gzip.open(file_name) as input_stream:
                term_id, name, is_lipid_species, is_lipid_class, is_carbon_chain = (
                    "",
                    "",
                    False,
                    False,
                    False,
                )
                domain, relations, synonyms = "", set(), []

                for line in input_stream:
                    line = line.decode("utf8").strip().strip(" ")

                    if line == "[Term]":
                        if term_id != "" and name != "" and name.find("inding") == -1:
                            if term_id in self.ontology_terms:
                                raise Exception(
                                    "Term id '%s' in lipid ontology already defined."
                                    % term_id
                                )

                            term = OntologyTerm(term_id, name, relations, domain)

                            if is_lipid_species:
                                if name not in self.lipids:
                                    self.lipids[name] = term

                            if is_lipid_class:
                                for synonym in synonyms:
                                    if synonym not in self.lipid_classes:
                                        self.lipid_classes[synonym] = []
                                    self.lipid_classes[synonym].append(term)

                            if is_carbon_chain:
                                for synonym in synonyms:
                                    if synonym not in self.carbon_chains:
                                        self.carbon_chains[synonym] = term

                            if len(domain) > 0:
                                self.search_terms[term_id] = set()
                                self.domains.add(domain)

                            self.ontology_terms[term_id] = term

                        (
                            term_id,
                            name,
                            is_lipid_species,
                            is_lipid_class,
                            is_carbon_chain,
                        ) = ("", "", False, False, False)
                        domain, relations, synonyms = "", set(), []

                    elif line[:4] == "id: ":
                        term_id = line[4:]

                    elif line[:6] == "name: ":
                        name = line[6:]

                    elif line[:6] == "is_a: ":
                        relation = line[6:].split(" ! ")[0]
                        relations.add(relation)

                        if relation == "LS:0000001":
                            is_lipid_class = True
                        elif relation == "LS:0000002":
                            is_lipid_species = True
                        elif relation == "LS:0000003":
                            is_carbon_chain = True

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
                            "Term id '%s' in lipid ontology already defined." % term_id
                        )

                    term = OntologyTerm(term_id, name, relations, domain)

                    if is_lipid_species:
                        if name not in self.lipids:
                            self.lipids[name] = term

                    if is_lipid_class:
                        for synonym in synonyms:
                            if synonym not in self.lipid_classes:
                                self.lipid_classes[synonym] = []
                            self.lipid_classes[synonym].append(term)

                    if is_carbon_chain:
                        for synonym in synonyms:
                            if synonym not in carbon_chains:
                                self.carbon_chains[synonym] = term

                    if len(domain) > 0:
                        self.search_terms[term_id] = set()
                        self.domains.add(domain)

                    self.ontology_terms[term_id] = term

                if "external" in self.domains:
                    self.domains.remove("external")

        except Exception as e:
            logger.info("ERROR: %s" % e)



    def set_background_lipids(self, lipid_list):
        for term_id in self.search_terms:
            self.ontology_terms[term_id].lipid_paths.clear()

        self.compute_event_occurrance(lipid_list)

        for term_id, lipid_set in self.search_terms.items():
            lipid_set.clear()
            lipid_set |= set(
                path_lipid for path_lipid in self.ontology_terms[term_id].lipid_paths
            )

        self.num_background = len(lipid_list)



    def recursive_event_adding(
        self, lipid_input_name, term_id, visited_terms, recursion = 0, path = None
    ):
        if recursion > 0:
            return

        if path == None:
            path = []

        if term_id in self.ontology_terms:
            term = self.ontology_terms[term_id]
            path.append(term_id)

            if lipid_input_name not in term.lipid_paths:
                term.lipid_paths[lipid_input_name] = {}
            lipid_paths = term.lipid_paths[lipid_input_name]
            if path[0] not in lipid_paths:
                lipid_paths[path[0]] = []
                lipid_path = lipid_paths[path[0]]
                for path_id in path:
                    lipid_path.append(path_id)

            next_recursion = recursion + (1 if term.domain in self.domains else 0)
            for parent_term_id in term.relations:
                if parent_term_id not in visited_terms:
                    visited_terms.add(parent_term_id)

                    self.recursive_event_adding(
                        lipid_input_name,
                        parent_term_id,
                        visited_terms,
                        next_recursion,
                        path,
                    )

            path.pop()



    def compute_event_occurrance(self, lipid_list):
        for lipid_input_name in lipid_list:
            try:
                lipid = self.lipid_parser.parse(lipid_input_name)
            except Exception as e:
                logging.info("Error: %s" % e)
                continue

            if lipid.lipid.info.level.value > LipidLevel.MOLECULAR_SPECIES.value:
                lipid.lipid.info.level = LipidLevel.MOLECULAR_SPECIES
            lipid.sort_fatty_acyl_chains()
            lipid_name = lipid.get_lipid_string()
            lipid_name_species = lipid.get_lipid_string(LipidLevel.SPECIES)
            lipid_name_class = lipid.get_extended_class()

            visited_terms = set()
            if lipid_name in self.lipids:
                term = self.lipids[lipid_name]
                if term.term_id not in visited_terms:
                    visited_terms.add(term.term_id)
                    self.recursive_event_adding(
                        lipid_input_name, term.term_id, visited_terms
                    )

            if lipid_name_species in self.lipids:
                term = self.lipids[lipid_name_species]
                if term.term_id not in visited_terms:
                    visited_terms.add(term.term_id)
                    self.recursive_event_adding(
                        lipid_input_name, term.term_id, visited_terms
                    )

            if lipid_name_class in self.lipid_classes:
                for term in self.lipid_classes[lipid_name_class]:
                    if term.term_id not in visited_terms:
                        visited_terms.add(term.term_id)
                        self.recursive_event_adding(
                            lipid_input_name, term.term_id, visited_terms
                        )

            if lipid != None and lipid.lipid != None:
                for fa in lipid.lipid.fa_list:
                    if fa.num_carbon == 0:
                        continue

                    fa_string = (
                        "L"
                        if fa.lipid_FA_bond_type
                        in {LipidFaBondType.LCB_REGULAR, LipidFaBondType.LCB_EXCEPTION}
                        else "C"
                    ) + fa.to_string(LipidLevel.MOLECULAR_SPECIES)
                    if fa_string in self.carbon_chains:
                        visited_terms.add(self.carbon_chains[fa_string].term_id)
                        self.recursive_event_adding(
                            lipid_input_name,
                            self.carbon_chains[fa_string].term_id,
                            visited_terms,
                        )



    def enrichment_analysis(self, target_list, enrichment_domains):
        if self.num_background == 0 or len(enrichment_domains) == 0:
            return

        result_list, target_set = [], set()
        for target in target_list:
            target_set.add(target)

        for term_id, lipid_set in self.search_terms.items():
            if (
                len(lipid_set) == 0
                or term_id not in self.ontology_terms
                or self.ontology_terms[term_id].domain not in enrichment_domains
            ):
                continue

            target_number = len(lipid_set & target_set)

            a, b, c, d = (
                target_number,
                len(lipid_set) - target_number,
                len(target_list) - target_number,
                self.num_background - len(lipid_set) - len(target_list) + target_number,
            )
            p_hyp = stats.fisher_exact([[a, b], [c, d]], alternative = "two-sided")[1]
            result_list.append(
                OntologyResult(
                    self.ontology_terms[term_id],
                    self.num_background,
                    len(lipid_set),
                    len(target_list),
                    p_hyp,
                )
            )

        return result_list
