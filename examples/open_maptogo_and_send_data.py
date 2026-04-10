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


import requests
import json


# The URL you want to send the POST request to
url = "https://maptogo.lifs-tools.org"

# Data to send in the POST request (as a dictionary)
data = {
    "separate_updown_switch": True,

    "all_lipids": ["PS 34:2", "PS 34:1", "PS 36:4", "PS 36:3", "PS 36:2", "PS 36:1", "PS 38:5", "PS 38:4", "PS 38:3", "PS 38:2", "PS 38:1", "PS 39:3", "PS 40:7", "PS 40:6", "PS 40:5", "PS 40:4", "PS 40:3", "PS 40:2", "PS 40:1", "PI 32:1", "PI 34:2", "PI 34:1", "PI 36:5", "PI 36:4", "PI 36:3", "PI 36:2", "PI 37:4", "PI 37:3", "PI 37:2", "PI 38:6", "PI 38:5", "PI 38:4", "PI 38:3", "PI 39:5", "PI 40:7", "PI 40:6", "PI 40:5", "TAG 44:1", "TAG 46:2", "TAG 47:2", "TAG 48:1", "TAG 48:2", "TAG 48:3", "TAG 49:2", "TAG 49:3", "TAG 50:2", "TAG 50:3", "TAG 50:4", "TAG 51:2", "TAG 51:3", "TAG 52:3", "TAG 52:4", "TAG 52:5", "TAG 53:3", "TAG 53:4", "TAG 54:3", "TAG 54:4", "TAG 54:5", "TAG 54:6", "TAG 54:7", "TAG 56:2", "TAG 56:3", "TAG 56:4", "TAG 56:5", "TAG 56:6", "TAG 56:7", "TAG 56:8", "TAG 58:2", "TAG 58:3", "TAG 58:4", "TAG 58:5", "TAG 58:6", "TAG 58:7", "TAG 58:8", "TAG 58:9", "TAG 60:3", "TAG 60:4", "PC-O 32:0", "PC-O 32:1", "PC-O 34:0", "PC-O 34:1", "PC-O 36:3", "PC-O 38:3", "PC-O 38:4", "PC 32:0", "PC 32:1", "PC 32:2", "PC 33:1", "PC 33:2", "PC 34:1", "PC 34:2", "PC 34:3", "PC 35:1", "PC 35:2", "PC 36:1", "PC 36:2", "PC 36:3", "PC 36:4", "PC 38:4", "PC 38:5", "PC 40:5", "PC 40:6", "LPC 14:0", "LPC 15:0", "LPC 16:1", "LPC 16:0", "LPC 17:0", "LPC 18:2", "LPC 18:1", "LPC 18:0", "LPC 20:4", "LPC 20:3", "LPC 20:1", "LPC 20:0", "LPC 22:6", "LPC 22:5", "LPC 22:4", "DG 30:0", "DG 32:0", "DG 32:1", "DG 34:0", "DG 34:1", "DG 34:2", "DG 36:1", "DG 36:2", "DG 36:4", "DG 38:5", "DG 43:6", "DG 46:1", "DG 48:1", "DG 50:1", "Cer 34:2", "Cer 34:1", "Cer 36:1", "Cer 40:2", "Cer 40:1", "Cer 42:3", "Cer 42:2", "Cer 42:1", "SM 32:1", "SM 33:1", "SM 34:1", "SM 34:0", "SM 35:1", "SM 36:1", "SM 38:2", "SM 38:1", "SM 39:1", "SM 40:2", "SM 40:1", "SM 41:1", "SM 42:2", "SM 42:1", "PG 32:0", "PG 32:1", "PG 34:1", "PG 34:2", "PG 34:3", "PG 34:4", "PG 36:1", "PG 36:3", "PG 36:4", "PG 38:2", "PG 38:3", "PG 38:4", "PG 38:5", "PG 38:6", "PG 40:3", "PG 40:5", "PG 40:6", "PG 40:7", "PG 42:7", "PE-O 32:2", "PE-O 34:1", "PE-O 34:2", "PE-O 34:3", "PE-O 35:2", "PE-O 35:4", "PE-O 36:1", "PE-O 36:2", "PE-O 36:3", "PE-O 36:4", "PE-O 36:5", "PE-O 36:6", "PE-O 37:4", "PE-O 37:5", "PE-O 37:6", "PE-O 38:1", "PE-O 38:2", "PE-O 38:3", "PE-O 38:4", "PE-O 38:5", "PE-O 38:6", "PE-O 38:7", "PE-O 40:2", "PE-O 40:3", "PE-O 40:5", "PE 32:1", "PE 32:2", "PE 33:1", "PE 34:1", "PE 34:2", "PE 35:1", "PE 35:2", "PE 36:1", "PE 36:2", "PE 36:3", "PE 36:4", "PE 36:6", "PE 37:2", "PE 37:3", "PE 37:4", "PE 38:2", "PE 38:3", "PE 38:4", "PE 38:5", "PE 38:6", "PE 38:7", "PE 39:5", "PE 39:6", "PE 40:1", "PE 40:2", "PE 40:4", "PE 40:5", "PE 40:6", "PE 40:7", "PE 42:7", "PA 32:0", "PA 32:1", "PA 34:0", "PA 34:1", "PA 34:2", "PA 35:2", "PA 36:2", "PA 36:1", "PA 38:4", "LPE 16:0", "LPE 16:1", "LPE 17:0", "LPE 18:0", "LPE 18:1", "LPE 20:4", "LPE 20:5", "LPE 22:5", "LPE 22:6", "LPA 16:0", "LPA 18:3", "LPA 18:0"],

    "upregulated_lipids": ["PS 34:2", "PS 34:1", "PS 36:4", "PS 36:3", "PS 36:2", "PS 36:1", "PS 38:5", "PS 38:4", "PS 38:2", "PS 38:1", "PS 39:3", "PS 40:7", "PS 40:6", "PS 40:5", "PS 40:4", "PS 40:3", "PS 40:2", "PS 40:1", "PI 32:1", "PI 34:2", "PI 34:1", "PI 36:5", "PI 36:4", "PI 36:3", "PI 36:2", "PI 37:4", "PI 37:3", "PI 37:2", "PI 38:6", "PI 38:5", "PI 38:4", "PI 38:3", "PI 39:5", "PI 40:7", "PI 40:6", "PI 40:5", "TAG 44:1", "TAG 46:2", "TAG 47:2", "TAG 48:1", "TAG 48:2", "TAG 49:2", "TAG 49:3", "TAG 50:2", "TAG 52:5", "TAG 54:7", "TAG 58:2", "PC 32:1", "PC 32:2", "PC 33:1", "PC 33:2", "PC 34:1", "PC 34:2", "PC 35:1", "PC 36:1", "LPC 14:0", "LPC 15:0", "LPC 16:1", "LPC 16:0", "LPC 17:0", "LPC 18:2", "LPC 18:1", "LPC 18:0", "LPC 20:3", "LPC 22:6", "DG 30:0", "DG 32:0", "DG 32:1", "DG 34:0", "DG 34:1", "DG 46:1", "Cer 34:2", "Cer 40:2", "Cer 42:3", "SM 33:1", "PG 32:0", "PG 32:1", "PG 34:2", "PG 34:3", "PE-O 32:2", "PE-O 36:2", "PE 32:1", "PE 32:2", "PE 33:1", "PE 34:1", "PE 34:2", "PE 36:4", "PE 36:6", "PE 38:7", "PE 40:6", "PA 32:0", "PA 32:1", "PA 34:0", "PA 34:1", "PA 34:2", "PA 35:2", "LPE 16:0", "LPE 17:0", "LPE 18:0", "LPE 18:1", "LPA 16:0", "LPA 18:0"],

    "downregulated_lipids": ["TAG 58:5", "TAG 58:6", "PC-O 34:1", "PC-O 36:3", "PC-O 38:3", "PC-O 38:4", "DG 36:2", "SM 34:0", "PG 34:4", "PG 36:1", "PG 36:3", "PG 36:4", "PG 38:2", "PG 38:3", "PG 38:4", "PG 38:5", "PG 40:3", "PG 40:5", "PG 40:6", "PG 40:7", "PG 42:7", "PE-O 34:2", "PE-O 35:4", "PE-O 36:3", "PE-O 36:4", "PE-O 36:5", "PE-O 37:5", "PE-O 38:1", "PE-O 38:2", "PE-O 38:3", "PE-O 38:4", "PE-O 38:5", "PE-O 40:3", "PE-O 40:5", "PE 37:2", "PE 38:2", "PE 38:3", "PE 39:5", "PE 40:1", "PE 40:2", "PE 40:4", "PE 40:7", "PE 42:7", "PA 36:1", "PA 38:4", "LPE 20:5", "LPE 22:5", "LPE 22:6", "LPA 18:3"],

    # "all_proteins": ...,
    # "regulated_proteins": ...,
    # "upregulated_proteins": ...,
    # "downregulated_proteins": ...,
    # "all_metabolites": ...,
    # "regulated_metabolites": ...,
    # "upregulated_metabolites": ...,
    # "downregulated_metabolites": ...,
    # "all_transcripts": ...,
    # "regulated_transcripts": ...,
    # "upregulated_transcripts": ...,
    # "downregulated_transcripts": ...,
    # "use_bounded_fatty_acyls": ...,
    "use_lipids": True,
    # "use_proteins": ...,
    # "use_metabolites": ...,
    # "use_transcripts": ...,

    #   Mus musculus -> "NCBITaxon:10090"
    #   Homo sapiens -> "NCBITaxon:9606"
    #   Bacillus cereus -> "NCBITaxon:405534"
    #   Saccharomyces cerevisiae -> "NCBITaxon:4932"
    #   Escherichia coli -> "NCBITaxon:562"
    #   Drosophila melanogaster -> "NCBITaxon:7227"
    #   Rattus norvegicus -> "NCBITaxon:10116"
    #   Bos taurus -> "NCBITaxon:9913"
    #   Caenorhabditis elegans -> "NCBITaxon:6239"
    #   Pseudomonas aeruginosa -> "NCBITaxon:287"
    #   Arabidopsis thaliana -> "NCBITaxon:3702"
    "organism": "NCBITaxon:9606",

    #   regulated_molecule_handling can be "molecule_handling_error", "molecule_handling_remove"
    # "molecule_handling": ...,

    #   regulated_molecule_handling can be "molecule_handling_error", "molecule_handling_remove",
    #   "molecule_handling_ignore"
    # "regulated_molecule_handling": ...,

    #   possible term_regulation:
    #   "greater" -> significant biomolocules are overrepresented in domain term (default)
    #   "less" -> significant biomolocules are underrepresented in domain term, means no effect on the term
    #   "two-sided" -> either way
    # "term_representation": ...,

    #   possible multiple_test_correction:
    #   "no" -> No correction
    #   "bonferroni" -> Bonferroni
    #   "fdr_bh" -> Benjamini/Hochberg (default)
    #   "fdr_by" -> Benjamini/Yekutieli
    #   "holm-sidak" -> Step down method using Sidak adjustment
    #   "holm" -> Step-down method using Bonferroni adjustments
    #   "simes-hochberg" -> Step-up method (independent)
    #   "fdr_tsbh" -> Two stage fdr correction
    # "test_method" = "fdr_bh"

    #   possible domains: "Biological process", "Cellular component", "Disease",
    #   "Metabolic and signalling pathway", "Molecular function", "Phenotype", "Physical or chemical properties", or any combination
    "domains": ["Molecular function"],

    # "use_upregulated_lipids": ...,
    # "use_downregulated_lipids": ...,
    # "use_upregulated_proteins": ...,
    # "use_downregulated_proteins": ...,
    # "use_upregulated_metabolites": ...,
    # "use_downregulated_metabolites": ...,
    # "use_upregulated_transcripts": ...,
    # "use_downregulated_transcripts": ...,
}
payload = json.dumps(data)

headers = {"Content-Type": "application/json"}
response = requests.post(f"{url}/submit", data=payload, headers=headers)
if response and "uid" in (response := json.loads(response.text)):
    print(response['uid'])
    import webbrowser
    webbrowser.open(f"{url}/?uid={response['uid']}")
else:
    print("An unexpected error occurred")
