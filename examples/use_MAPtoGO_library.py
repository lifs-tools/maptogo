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

import sys
import os
sys.path.append(os.path.abspath(".."))
import pandas as pd
from EnrichmentDataStructure import *



## Load example dataset
print("Load example dataset")
df = pd.read_excel("../Data/examples.xlsx", "Simplex")
background_lipids_list = list(df["BackgroundLipids"].dropna())
regulated_lipids_list = list(df["RegulatedLipids"].dropna())
background_proteins_list = list(df["BackgroundProteins"].dropna())
regulated_proteins_list = list(df["RegulatedProteins"].dropna())
background_metabolites_list = list(df["BackgroundMetabolites"].dropna())
regulated_metabolites_list = list(df["RegulatedMetabolites"].dropna())
organism_taxonomy = df["Organism"][0]



## load correct ontology
print(f"Load {organism_taxonomy} ontology")
ontology = EnrichmentOntology(f"../Data/ontology_{organism_taxonomy.split(":")[1]}.gz")



## define which omics layer is included in analysis (lipidomics, proteomics, metabolomics, transcriptomics)
omics_included = [True, True, True, False]
omics_lists = [
    background_lipids_list,
    regulated_lipids_list,
    background_proteins_list,
    regulated_proteins_list,
    background_metabolites_list,
    regulated_metabolites_list,
    [],
    [],
]
ignore_unrecognizable_molecules = MOLECULE_HANDLING_REMOVE
ignore_unknown = MOLECULE_HANDLING_REMOVE
use_bounded_fatty_acyls = False



## parse all biomolocule entries
print("Parse all biomolocule entries")
try:
    (
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
    ) = check_user_input(
        omics_included,
        omics_lists,
        ontology,
        ignore_unrecognizable_molecules,
        ignore_unknown,
    )
    num_background = len(lipidome) + len(proteome) + len(metabolome) + len(transcriptome)
except Exception as error_message:
    print(error_message)
    exit(-1)



## set enrichment background
print("Set enrichment background")
(
    search_terms,
    all_parent_nodes,
) = ontology.set_background(
    lipid_dict = lipidome,
    protein_set = proteome,
    metabolite_set = metabolome,
    transcript_set = transcriptome,
    use_bounded_fatty_acyls = use_bounded_fatty_acyls,
)



## run enrichment analysis
print("Run enrichment analysis")
# possible domains: "Biological process", "Cellular component", "Disease", "Metabolic and signalling pathway", "Molecular function", "Phenotype", "Physical or chemical properties"
domains = {"Biological process"}

# possible term_regulation: "greater", "less", "two-sided"
term_regulation = "greater"

# possible multiple_test_correction:
# "no" -> No correction
# "bonferroni" -> Bonferroni
# "fdr_bh" -> Benjamini/Hochberg
# "fdr_by" -> Benjamini/Yekutieli
# "holm-sidak" -> Step down method using Sidak adjustment
# "holm" -> Step-down method using Bonferroni adjustments
# "simes-hochberg" -> Step-up method (independent)
# "fdr_tsbh" -> Two stage fdr correction
multiple_test_correction = "fdr_bh"

results = ontology.enrichment_analysis(
    search_terms,
    num_background,
    target_set,
    domains,
    term_regulation,
    multiple_test_correction,
)


### printing the results in p-value reverse order
for result in results[::-1]:
    expected = round((result.fisher_data[0] + result.fisher_data[1]) * (result.fisher_data[0] + result.fisher_data[2]) / sum(result.fisher_data))
    print("domain", " | ".join(result.term.domain))
    print("term", result.term.name)
    print("termid", result.term.term_id_str)
    print("count", f"{result.fisher_data[0]} ({expected}) / {len(result.source_terms)}")
    print("pvalue", result.pvalue_corrected)
    print("log_odds_ratio", result.lor)
    print()
