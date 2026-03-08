"""
MIT License

Copyright (c) 2026 Dominik Kopczynski  -  dominik.kopczynski {at} univie.ac.at
                   Cristina Coman  -  cristina.coman {at} univie.ac.at
                   Nils Hoffmann  -  n.hoffmann {at} fz-juelich.de
                   Robert Ahrends  -  robert.ahrends {at} univie.ac.at

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import gzip
from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.domain.LipidFaBondType import LipidFaBondType
import scipy.stats as stats
import logging
import json
import ctypes
import pathlib
from collections import defaultdict
import traceback
import pickle
import os
import csv
import sys
from enum import Enum
from statsmodels.stats.multitest import multipletests
from itertools import groupby


lipid_parser = LipidParser()
MOLECULE_HANDLING_ERROR = "molecule_handling_error"
MOLECULE_HANDLING_REMOVE = "molecule_handling_remove"
MOLECULE_HANDLING_IGNORE = "molecule_handling_ignore"


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
    with gzip.open(CHEBI_synonym_table_filename, "rt") as input_stream:
        CHEBI_synonym_table = {tokens[0].lower(): tokens[1] for line in input_stream.read().split("\n") if len(line) > 0 and (tokens := line.split("\t")) and len(tokens) >= 2}
else:
    logger.warning("No ChEBI synonyms table found")




def check_user_input(
    omics_included,
    omics_lists,
    ontology,
    ignore_unrecognizable_molecules,
    ignore_unknown,
):
    with_lipids, with_proteins, with_metabolites, with_transcripts = omics_included
    target_set = set()
    lipidome, regulated_lipids = {}, set()
    proteome, regulated_proteins = set(), set()
    metabolome, regulated_metabolites = set(), set()
    transcriptome, regulated_transcripts = set(), set()
    background_list = []

    (
        all_lipids_list,
        regulated_lipids_list,
        all_proteins_list,
        regulated_proteins_list,
        all_metabolites_list,
        regulated_metabolites_list,
        all_transcripts_list,
        regulated_transcripts_list,
    ) = omics_lists

    if with_lipids:
        if type(all_lipids_list) == str: all_lipids_list = all_lipids_list.split("\n")
        elif len(all_lipids_list) == 0:
            raise Exception("No background lipids are defined.")
        if type(regulated_lipids_list) == str: regulated_lipids_list = regulated_lipids_list.split("\n")
        elif len(regulated_lipids_list) == 0:
            raise Exception("No regulated lipids are defined.")

        for lipid_name in all_lipids_list:
            if len(lipid_name) == 0: continue
            try:
                lipid = lipid_parser.parse(lipid_name)
            except Exception as e:
                lower_lipid_name = lipid_name.lower()
                if lower_lipid_name in ontology.unspecific_lipids:
                    lipid = ontology.unspecific_lipids[lower_lipid_name]
                else:
                    if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                        lipid = None
                    elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                        continue
                    else:
                        raise Exception(f"Lipid name '{lipid_name}' unrecognizable! Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.")
            lipidome[lipid_name] = lipid

        for lipid_name in regulated_lipids_list:
            if len(lipid_name) == 0: continue
            if lipid_name not in lipidome:
                if ignore_unknown == MOLECULE_HANDLING_REMOVE: continue
                raise Exception(f"The regulated lipid '{lipid_name}' does not occur in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of non-background regulated molecules' setting.")
            regulated_lipids.add(lipid_name)

        if len(lipidome) == 0:
            raise Exception("No background lipid left after lipid recognition.")

        if len(regulated_lipids) == 0:
            raise Exception("No regulated lipid left after lipid recognition.")

        if len(regulated_lipids) > len(lipidome):
            raise Exception("Length of regulated lipid list must be smaller than background list.")

        left_lipids = regulated_lipids - lipidome.keys()
        if len(left_lipids) > 0:
            if ignore_unknown == MOLECULE_HANDLING_REMOVE:
                for lipid_name in left_lipids:
                    del lipidome[lipid_name]
            else:
                raise Exception("The regulated lipid" + (' ' if len(left_lipids) == 1 else 's ') + "'" + "', '".join(left_lipids) + ("' does" if len(left_lipids) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of non-background regulated molecules' setting.")

        target_set |= regulated_lipids
        background_list += [{"value": k, "label": k} for k in lipidome.keys()]

    if with_proteins:
        if type(all_proteins_list) == str: all_proteins_list = all_proteins_list.split("\n")
        elif len(all_proteins_list) == 0:
            raise Exception("No background proteins are defined.")

        if type(regulated_proteins_list) == str: regulated_proteins_list = regulated_proteins_list.split("\n")
        elif len(regulated_proteins_list) == 0:
            raise Exception("No regulated proteins are defined.")

        proteome = set(protein for protein in all_proteins_list if len(protein) > 0)
        regulated_proteins = set(protein for protein in regulated_proteins_list if len(protein) > 0)

        background_list += [{"value": pp, "label": p + (' (' + ontology.proteins[pp].name + ')' if (pp in ontology.proteins) else '')} for p in proteome if (pp := "UNIPROT:" + p)]
        proteome = set(p.split("-")[0] for p in proteome)
        regulated_proteins = set(rp.split("-")[0] for rp in regulated_proteins)
        left_proteins = proteome - ontology.clean_protein_ids
        if len(left_proteins) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                proteome -= left_proteins
            else:
                raise Exception("The protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' is" if len(left_proteins) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.")

        left_proteins = regulated_proteins - ontology.clean_protein_ids
        if len(left_proteins) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                regulated_proteins -= left_proteins
            else:
                raise Exception("The protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' is" if len(left_proteins) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.")

        left_proteins = regulated_proteins - proteome
        if len(left_proteins) > 0:
            if ignore_unknown == MOLECULE_HANDLING_REMOVE:
                regulated_proteins -= left_proteins
            else:
                raise Exception("The regulated protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' does" if len(left_proteins) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of non-background regulated molecules' setting.")

        if len(proteome) == 0:
            raise Exception("No background protein left after protein recognition.")

        if len(regulated_proteins) == 0:
            raise Exception("No regulated protein left after protein recognition.")

        if len(regulated_proteins) > len(proteome):
            raise Exception("Length of regulated protein list must be smaller than background list.")

        proteome = set([f"UNIPROT:{protein}" for protein in proteome])
        regulated_proteins = set([f"UNIPROT:{protein}" for protein in regulated_proteins])
        target_set |= regulated_proteins

    if with_metabolites:
        if type(all_metabolites_list) == str: all_metabolites_list = all_metabolites_list.split("\n")
        elif len(all_metabolites_list) == 0:
            raise Exception("No background metabolites are defined.")

        if type(regulated_metabolites_list) == str: regulated_metabolites_list = regulated_metabolites_list.split("\n")
        elif len(regulated_metabolites_list) == 0:
            raise Exception("No regulated metabolites are defined.")

        metabolome = set(metabolite for metabolite in all_metabolites_list if len(metabolite) > 0)
        regulated_metabolites = set(metabolite for metabolite in regulated_metabolites_list if len(metabolite) > 0)

        background_list += [{"value": m, "label": m} for m in metabolome]
        left_metabolites = metabolome - ontology.clean_metabolite_ids - ontology.metabolites.keys()
        left_metabolites -= set([m for m in left_metabolites if m.lower() in ontology.metabolite_names.keys()])
        if len(left_metabolites) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                metabolome -= left_metabolites
            else:
                raise Exception("The metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' is" if len(left_metabolites) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.")

        left_metabolites = regulated_metabolites - ontology.clean_metabolite_ids - ontology.metabolites.keys()
        left_metabolites -= set([m for m in left_metabolites if m.lower() in ontology.metabolite_names.keys()])
        if len(left_metabolites) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                regulated_metabolites -= left_metabolites
            else:
                raise Exception("The metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' is" if len(left_metabolites) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.")

        left_metabolites = regulated_metabolites - metabolome
        if len(left_metabolites) > 0:
            if ignore_unknown == MOLECULE_HANDLING_REMOVE:
                regulated_metabolites -= left_metabolites
            else:
                raise Exception("The regulated metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' does" if len(left_metabolites) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.")

        if len(metabolome) == 0:
            raise Exception("No background metabolite left after metabolite recognition.")

        if len(regulated_metabolites) == 0:
            raise Exception("No regulated metabolite left after metabolite recognition.")

        if len(regulated_metabolites) > len(metabolome):
            raise Exception("Length of regulated metabolite list must be smaller than background list.")

        metabolome = set([f"CHEBI:{metabolite}" if (type(metabolite) == int or not metabolite.startswith("CHEBI:")) and metabolite.lower() not in ontology.metabolite_names else metabolite for metabolite in metabolome])
        regulated_metabolites = set([f"CHEBI:{metabolite}" if (type(metabolite) == int or not metabolite.startswith("CHEBI:")) and metabolite.lower() not in ontology.metabolite_names else metabolite for metabolite in regulated_metabolites])
        target_set |= regulated_metabolites

    if with_transcripts:
        if type(all_transcripts_list) == str: all_transcripts_list = all_transcripts_list.split("\n")
        elif len(all_transcripts_list) == 0:
            raise Exception("No background transcripts are defined.")

        if type(regulated_transcripts_list) == str: regulated_transcripts_list = regulated_transcripts_list.split("\n")
        elif len(regulated_transcripts_list) == 0:
            raise Exception("No regulated transcripts are defined.")

        transcriptome = set(transcript for transcript in all_transcripts_list if len(transcript) > 0)
        regulated_transcripts = set(transcript for transcript in regulated_transcripts_list if len(transcript) > 0)

        background_list += [{"value": t, "label": t + (' (' + ontology.transcripts[tt].name + ')' if (tt in ontology.transcripts) else '')} for t in transcriptome if (tt := t.split(".")[0])]
        transcript_keys = set(ontology.transcripts.keys())
        left_transcripts = set(t for t in transcriptome if t.split(".")[0] not in transcript_keys)
        if len(left_transcripts) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                transcriptome -= left_transcripts
            else:
                raise Exception("The transcript" + (' ' if len(left_transcripts) == 1 else 's ') + "'" + "', '".join(left_transcripts) + ("' is" if len(left_transcripts) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.")

        left_transcripts = set(t for t in regulated_transcripts if t.split(".")[0] not in transcript_keys)
        if len(left_transcripts) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                regulated_transcripts -= left_transcripts
            else:
                raise Exception("The transcript" + (' ' if len(left_transcripts) == 1 else 's ') + "'" + "', '".join(left_transcripts) + ("' is" if len(left_transcripts) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.")

        left_transcripts = regulated_transcripts - transcriptome
        if len(left_transcripts) > 0:
            if ignore_unknown == MOLECULE_HANDLING_REMOVE:
                regulated_transcripts -= left_transcripts
            else:
                raise Exception("The regulated transcript" + (' ' if len(left_transcripts) == 1 else 's ') + "'" + "', '".join(left_transcripts) + ("' does" if len(left_transcripts) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.")
        if len(transcriptome) == 0:
            raise Exception("No background transcript left after transcript recognition.")

        if len(regulated_transcripts) == 0:
            raise Exception("No regulated transcript left after transcript recognition.")

        if len(regulated_transcripts) > len(transcriptome):
            raise Exception("Length of regulated transcript list must be smaller than background list.")

        target_set |= regulated_transcripts

    return (
        target_set,
        lipidome,
        regulated_lipids,
        proteome,
        regulated_proteins,
        metabolome,
        regulated_metabolites,
        transcriptome,
        regulated_transcripts,
        background_list,
    )


EXCLUDE_TERMS = {"", "external"}
class OntologyTerm:
    __slots__ = ("term_id", "term_id_str", "term_type", "name", "domain", "relations", "categories", "synonyms")

    def __init__(self, _term_id, _name, _relations, _synonyms, _domain, _categories):
        self.term_type = TermType(int(_relations[moea_pos + 5 : moea_pos + 12])) if (moea_pos := _relations.find("MOEA:00000")) > -1 else TermType.UNCLASSIFIED_TERM
        self.term_id = sorted(_term_id.split("|"))
        self.term_id_str = "|".join(self.term_id)
        self.name = _name
        self.domain = set(_domain.split("|")) - EXCLUDE_TERMS
        self.relations = _relations.split("|")
        self.categories = set(_categories.split("|")) - EXCLUDE_TERMS
        self.synonyms = _synonyms.split("|")

    def get_term_id(self, space = False):
        return " | ".join(sorted(self.term_id))



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
    def __init__(self, file_name, ontology_name = "undefined"):
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

        try:
            csv.field_size_limit(sys.maxsize)
            with gzip.open(file_name, mode="rt", encoding="utf-8", newline = "") as f:
                for term in [OntologyTerm(*row) for row in csv.reader(f, delimiter = "\t") if len(row) == 6]:
                    for d in term.domain: self.domains.add(d)
                    for t_id in term.term_id: self.ontology_terms[t_id] = term

                    str_term_id, name = term.term_id_str, term.name
                    match term.term_type:
                        case TermType.LIPID_CLASS: # lipid class
                            if name not in self.lipid_classes: self.lipid_classes[name] = []
                            self.lipid_classes[name].append(term)
                            for synonym in term.synonyms:
                                if synonym not in self.lipid_classes: self.lipid_classes[synonym] = []
                                self.lipid_classes[synonym].append(term)

                        case TermType.LIPID_SPECIES: # lipid species
                            if name not in self.lipids: self.lipids[name] = term

                        case TermType.CARBON_CHAIN: # carbon chain
                            if name not in self.carbon_chains: self.carbon_chains[name] = term
                            for synonym in term.synonyms:
                                self.carbon_chains[synonym] = term

                        case TermType.UNSPECIFIC_LIPID:
                            if name.lower() not in self.unspecific_lipids: self.unspecific_lipids[name.lower()] = term
                            for synonym in term.synonyms:
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


        except Exception as e:
            logger.error("".join(traceback.format_tb(e.__traceback__)))
            logger.error(e)

        self.clean_protein_ids = set([key.replace("UNIPROT:", "") for key in self.proteins.keys()])
        self.clean_metabolite_ids = set([key.replace("CHEBI:", "") for key in self.metabolites.keys()])
        self.metabolite_names = {term.name.lower(): term for term in self.metabolites.values()}
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
        try:
            ontology_keys, ontology_terms = set(self.ontology_terms), self.ontology_terms
            for term in self.ontology_terms.values():
                if not term.relations or type(term.relations[0]) == OntologyTerm: continue
                term.relations = [
                    k for k, _ in groupby(
                        sorted(
                            [ontology_terms[t] for t in term.relations if t in ontology_keys],
                            key = lambda term: term.term_id_str
                        )
                    )
                ]
        except Exception as e:
            logger.error("".join(traceback.format_tb(e.__traceback__)))
            logger.error(e)


    def set_background(
        self,
        lipid_dict = {},
        protein_set = set(),
        metabolite_set = set(),
        transcript_set = {},
        use_bounded_fatty_acyls = False,
    ):
        search_terms = defaultdict(list)
        all_parent_nodes = {}
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

            if use_bounded_fatty_acyls:
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

        return search_terms, all_parent_nodes


    def enrichment_analysis(self, search_terms, num_background, target_set, enrichment_domains, term_regulation = "greater", multiple_test_correction = "no"):
        if len(target_set) == 0 or num_background < 2 or len(enrichment_domains) == 0: return []

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

        results = [result for result in result_list if result != None]
        if multiple_test_correction != "no" and len(results) > 1:
            pvalues = [r.pvalue for r in results]
            pvalues = multipletests(pvalues, method = multiple_test_correction)[1]
            for pvalue, r in zip(pvalues, results): r.pvalue_corrected = pvalue
        results.sort(key = lambda row: (row.pvalue_corrected, row.term.name))
        return results
