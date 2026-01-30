from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidLevel import LipidLevel
import gzip
import sys
import pickle
import os
import pandas as pd
import zipfile

with zipfile.ZipFile("Data/BioDolphin_vr1.1.zip") as z:
    with z.open("BioDolphin_vr1.1/BioDolphin_vr1.1.txt") as f:
        df_dolphin = pd.read_csv(f, sep = "\t")
        df_dolphin = df_dolphin[["BioDolphinID", "protein_UniProt_ID", "lipid_ChEBI"]]
        df_dolphin = df_dolphin[df_dolphin["protein_UniProt_ID"].str.len() > 0]
        df_dolphin = df_dolphin[df_dolphin["lipid_ChEBI"].str.len() > 0]


all_species = len(sys.argv) > 1 and sys.argv[1] == "all"
parser = LipidParser()

PROTEIN_CATEGORY_RULES = {
    "Enzyme": [
        "catalytic activity", "EC=", "oxidoreductase", "transferase",
        "hydrolase", "lyase", "isomerase", "ligase"
    ],
    "Structural": [
        "structural molecule activity", "cytoskeleton", "extracellular matrix",
        "filament", "scaffold", "collagen"
    ],
    "Transport": [
        "transporter activity", "channel activity", "carrier activity",
        "oxygen transport", "ion transport", "permease"
    ],
    "Storage": [
        "nutrient reservoir activity", "storage", "ferritin", "casein", "ovalbumin"
    ],
    "Signaling": [
        "signal transduction", "receptor", "hormone", "ligand",
        "cytokine", "transcription factor", "gpcr"
    ],
    "Defense": [
        "immune response", "antibody", "defensin", "complement",
        "toxin", "antimicrobial", "venom"
    ],
    "Motor": [
        "motor activity", "myosin", "kinesin", "dynein", "movement"
    ]
}


# Map EC class (first digit) to category
EC_CATEGORY = {
    "1": "Oxidation–Reduction (Oxidoreductases)",
    "2": "Transfer Reactions (Transferases)",
    "3": "Hydrolysis (Hydrolases)",
    "4": "Group Elimination (Lyases)",
    "5": "Isomerization (Isomerases)",
    "6": "Ligation (Ligases)"
}


gene_biotype_dict = {
    "protein_coding": "Standard protein-coding gene",
    "IG_V_gene": "Immunoglobulin Variable gene segment",
    "IG_C_gene": "Immunoglobulin Constant gene segment",
    "IG_D_gene": "Immunoglobulin Diversity gene segment",
    "TR_V_gene": "T-cell Receptor Variable gene segment",
    "TR_C_gene": "T-cell Receptor Constant gene segment",
    "TR_D_gene": "T-cell Receptor Diversity gene segment",
    "TR_J_gene": "T-cell Receptor Joining gene segment",
    "miRNA": "MicroRNA",
    "rRNA": "Ribosomal RNA",
    "snRNA": "Small Nuclear RNA",
    "snoRNA": "Small Nucleolar RNA",
    "scRNA": "Small Cajal body-specific RNA",
    "sRNA": "Small RNA",
    "misc_RNA": "Miscellaneous RNA",
    "lincRNA": "Long Intergenic Non-Coding RNA",
    "lncRNA": "Long Non-Coding RNA",
    "antisense": "Antisense RNA",
    "sense_intronic": "Sense Intronic RNA",
    "sense_overlapping": "Sense Overlapping RNA",
    "bidirectional_promoter_lncRNA": "Bidirectional Promoter Long Non-Coding RNA",
    "processed_transcript": "Processed Transcript",
    "pseudogene": "Generic pseudogene",
    "IG_pseudogene": "Immunoglobulin Pseudogene",
    "TR_pseudogene": "T-cell Receptor Pseudogene",
    "miRNA_pseudogene": "MicroRNA Pseudogene",
    "rRNA_pseudogene": "Ribosomal RNA Pseudogene",
    "snRNA_pseudogene": "Small Nuclear RNA Pseudogene",
    "snoRNA_pseudogene": "Small Nucleolar RNA Pseudogene",
    "scRNA_pseudogene": "Small Cajal body-specific RNA Pseudogene",
    "misc_RNA_pseudogene": "Miscellaneous RNA Pseudogene",
    "tRNA_pseudogene": "Transfer RNA Pseudogene",
    "Mt_tRNA_pseudogene": "Mitochondrial Transfer RNA Pseudogene",
    "TEC": "To be Experimentally Confirmed",
    "nonsense_mediated_decay": "Transcript subject to Nonsense-Mediated Decay",
    "ribozyme": "Ribozyme",
    "scaRNA": "Small Cajal body-specific RNA",
    "scRNA": "Small Cajal body-specific RNA",
    "Mt_rRNA": "Mitochondrial Ribosomal RNA",
    "Mt_tRNA": "Mitochondrial Transfer RNA",
    "misc_RNA": "Miscellaneous RNA",
    "sRNA": "Small RNA",
    "scaRNA": "Small Cajal body-specific RNA",
}

CHEBI_CATEGORY_MAPPING = {
    "CHEBI:18059": "Lipid",
    "CHEBI:33708": "Amino acid",
    "CHEBI:16670": "Peptide",
    "CHEBI:16646": "Carbohydrate",
    "CHEBI:33838": "Nucleoside",
    "CHEBI:22315": "Alkaloid",
    "CHEBI:26873": "Terpenoid",
    "CHEBI:35341": "Steroid",
    "CHEBI:134179": "Volatile organic compound",
    "CHEBI:64709": "Organic acid",
    "CHEBI:36976": "Nucleotide",
    "CHEBI:47916": "Flavonoids",
    "CHEBI:166890": "Phenolic acid",
    "CHEBI:22580": "Anthraquinone",
    "CHEBI:28794": "Coumarins",
}




if all_species:
    print("Creating all ontologies")
    species = {
        'Homo sapiens': "9606",
        'Mus musculus': "10090",
        'Bacillus cereus': "405534",
        'Saccharomyces cerevisiae': "4932",
        'Escherichia coli': "562",
        'Drosophila melanogaster': "7227",
        'Rattus norvegicus': "10116",
        'Bos taurus': "9913",
        'Caenorhabditis elegans': "6239",
        'Pseudomonas aeruginosa': "287",
        'Arabidopsis thaliana': "3702",
    }
    ensembl_files = [
        ("Data/Homo_sapiens.uniprot.tsv.gz", "9606", "Data/Homo_sapiens.all.tsv.gz"),
        ("Data/Mus_musculus.uniprot.tsv.gz", "10090", "Data/Mus_musculus.all.tsv.gz"),
        ("Data/Saccharomyces_cerevisiae.uniprot.tsv.gz", "4932", "Data/Saccharomyces_cerevisiae.all.tsv.gz"),
        ("Data/Drosophila_melanogaster.uniprot.tsv.gz", "7227", "Data/Drosophila_melanogaster.all.tsv.gz"),
        ("Data/Rattus_norvegicus.uniprot.tsv.gz", "10116", "Data/Rattus_norvegicus.all.tsv.gz"),
        ("Data/Bos_taurus.uniprot.tsv.gz", "9913", "Data/Bos_taurus.all.tsv.gz"),
        ("Data/Caenorhabditis_elegans.uniprot.tsv.gz", "6239", "Data/Caenorhabditis_elegans.all.tsv.gz"),
    ]
    reactomes = {
        "9606": "HSA",
        "10090": "MMU",
        "4932": "SCE",
        "562": "ECO",
        "7227": "DME",
        "10116": "RNO",
        "9913": "BTA",
        "6239": "CEL",
        "287": None,
        "3702": None,
        "405534": None,
    }

else:
    print("Creating Human & Mouse ontology")
    species = {
        'Homo sapiens': "9606",
        'Mus musculus': "10090",
        'Rattus norvegicus': "10116",
    }
    ensembl_files = [
        ("Data/Homo_sapiens.uniprot.tsv.gz", "9606", "Data/Homo_sapiens.all.tsv.gz"),
        ("Data/Mus_musculus.uniprot.tsv.gz", "10090", "Data/Mus_musculus.all.tsv.gz"),
        ("Data/Rattus_norvegicus.uniprot.tsv.gz", "10116", "Data/Rattus_norvegicus.all.tsv.gz"),
    ]
    reactomes = {
        "9606": "HSA",
        "10090": "MMU",
        "10116": "RNO",
    }

species_set = set(species.values())


class Term:
    def __init__(self, _id, _name, _relations = None, _synonyms = None, _namespace = None, _xref = None, upper = False, _categories = None):
        if _relations == None: _relations = set()
        if _synonyms == None: _synonyms = set()
        if _namespace == None: _namespace = set()
        if _xref == None: _xref = set()
        if _categories == None: _categories = set()

        self.id = set(_id) if type(_id) in {list, set} else {_id}
        self.name = _name[0].upper() + _name[1:] if upper else _name
        self.relations = set(_relations) if type(_relations) in {list, set} else {_relations}
        self.synonyms = set(_synonyms) if type(_synonyms) in {list, set} else {_synonyms}
        self.namespace = set(_namespace) if type(_namespace) in {list, set} else {_namespace}
        self.xref = set(_xref) if type(_xref) in {list, set} else {_xref}
        self.categories = set(_categories) if type(_categories) in {list, set} else {_categories}


    def get_first_id(self):
        return list(self.id)[0]
        
    def joining(self, a_set):
        return "|".join(sorted(list(a_set)))

    def to_string(self):
        return f"{self.joining(self.id)}\t{self.name}\t{self.joining(self.relations)}\t{self.joining(self.synonyms)}\t{self.joining(self.namespace)}\t{self.joining(self.categories)}\n"


    def copy(self):
        return Term(set(self.id), self.name, set(self.relations), set(self.synonyms), set(self.namespace), set(self.xref), _categories = set(self.categories))


    def merge(self, copy):
        self.id |= copy.id
        self.relations |= copy.relations
        self.synonyms |= copy.synonyms
        self.namespace |= copy.namespace
        self.xref |= copy.xref
        self.categories |= copy.categories


lipid_maps_terms = {}
with open("Data/LM.csv") as fin:
    print("Readin LM")
    for line in fin:
        line = line.strip()
        if len(line) == 0: continue
        tokens = line.split("\t")
        if len(tokens) < 3 or len(tokens[0]) == 0 or len(tokens[1]) == 0: continue
        lipid_maps_terms[tokens[0]] = Term(tokens[0], tokens[2], {"MOEA:0000004", "CHEBI:" + tokens[1]}, tokens[3:] if len(tokens) > 3 else set())


def get_terms(term_file, prefix, upper = False, with_relationship = False):
    ontology_terms = {}

    opener = open if term_file.lower()[-3:] != ".gz" else gzip.open

    with opener(term_file, "rt") as obo:
        term_id, name, relations, synonyms, namespace, is_obsolete, xref = "", "", set(), [], "", False, set()
        for i, line in enumerate(obo):
            line = line.strip()
            if len(line) == 0: continue
            if i > 0 and i % 100000 == 0: print(term_file, i)

            if line == "[Term]":
                if len(term_id) > 0 and len(name) > 0 and not is_obsolete and term_id.startswith(prefix):
                    ontology_terms[term_id] = Term(term_id, name, relations, synonyms, namespace, xref, upper)
                term_id, name, relations, synonyms, is_obsolete, namespace, xref = "", "", set(), [], False, "", set()

            elif line.startswith("id:"):
                term_id = line[3:].strip(" ")

            elif line.startswith("name:"):
                name = line[5:].strip(" ")

            elif line.startswith("is_a:"):
                relation = line[5:].split(" ! ")[0].strip(" ").split(" ")[0]
                relations.add(relation)

            elif with_relationship and line.startswith("relationship: "):
                relation = "GO:" + line[14:].split("GO:")[1].split(" ! ")[0]
                relations.add(relation)

            elif line.startswith( "synonym:"):
                synonyms.append(line.split("\"")[1])

            elif line.startswith("namespace: "):
                namespace = line[11:]

            elif line.startswith("xref: "):
                xref.add(line[6:].split(" ")[0])

            elif line.startswith("is_obsolete:"):
                is_obsolete = True

        if len(term_id) > 0 and len(name) > 0 and not is_obsolete and term_id.startswith(prefix):
            ontology_terms[term_id] = Term(term_id, name, relations, synonyms, namespace, xref, upper)

    return ontology_terms



def get_lipid_terms():
    # import LION ontology
    term_id, name, relations, synonyms, is_class, is_species, namespace, category = "", "", set(), [], False, False, "", None
    ontology_terms = {}
    lipid_class_terms = {}
    lipid_species_terms = {}
    if True:
        with open("Data/LION.obo", "rt") as obo:
            for i, line in enumerate(obo):
                line = line.strip()
                if len(line) == 0: continue
                if i > 0 and i % 10000 == 0: print("LION:", i)

                if line == "[Term]":
                    if len(term_id) > 0 and len(name) > 0:
                        term = Term(term_id, name, relations, synonyms, namespace, _categories = category)
                        ontology_terms[term_id] = term

                        if is_class:
                            if name not in lipid_class_terms: lipid_class_terms[name] = []
                            lipid_class_terms[name].append(term)
                            for synonym in synonyms:
                                if synonym not in lipid_class_terms: lipid_class_terms[synonym] = []
                                lipid_class_terms[synonym].append(term)

                        if is_species:
                            if name not in lipid_species_terms: lipid_species_terms[name] = []
                            lipid_species_terms[name].append(term)

                    term_id, name, relations, synonyms, is_class, is_species, namespace, category = "", "", set(), [], False, False, "", None

                elif line.startswith("id:"):
                    term_id = line[3:].strip(" ")

                elif line.startswith("name:"):
                    name = line[5:].strip(" ")
                    try:
                        parse_name = name
                        if parse_name[0] != "C" or parse_name[1] not in {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"}:
                            if parse_name[:3] == "FFA": parse_name = parse_name[1:]
                            lipid = parser.parse(parse_name)
                            # if lipid.lipid.info.level.value > LipidLevel.MOLECULAR_SPECIES.value:
                            #     lipid.lipid.info.level = LipidLevel.MOLECULAR_SPECIES
                            # lipid.sort_fatty_acyl_chains()
                            name = lipid.get_lipid_string()
                            category = lipid.get_lipid_string(LipidLevel.CATEGORY)
                            is_species = True
                            relations.add("MOEA:0000002")
                    except Exception as e:
                        pass


                elif line.startswith("is_a:"):
                    relation = line[5:].split(" ! ")[0].strip(" ")
                    relations.add(relation)
                    if relation == "MOEA:0000001": is_class = True
                    elif relation == "MOEA:0000002": is_species = True

                elif line.startswith("synonym:"):
                    synonyms.append(line.split("\"")[1])

                elif line.startswith("namespace: "):
                    namespace = line[11:]


    if len(term_id) > 0 and len(name) > 0:
        term = Term(term_id, name, relations, synonyms, namespace, _categories = category)
        ontology_terms[term_id] = term

        if is_class:
            for synonym in synonyms:
                if synonym not in lipid_class_terms: lipid_class_terms[synonym] = []
                lipid_class_terms[synonym].append(term)

        if is_species:
            if name not in lipid_species_terms: lipid_species_terms[name] = []
            lipid_species_terms[name].append(term)




    # import chebi to lipid dictionary
    chebi_to_lipid = {}
    UNDEFINED = "UNDEFINED"
    for file_prefix in ["", "ChEBI_", "LipidMaps_", "SwissLipids_"]:
        for line in open(f"Data/{file_prefix}chebi_to_lipid.csv", "rt").read().split("\n"):
            if len(line) < 1: continue
            tokens = line.split("\t")
            if len(tokens) < 2: continue

            the_lipid = None
            for lipid_name in tokens[1:]:
                if lipid_name == "": continue
                try:
                    lipid = parser.parse(lipid_name)
                    if not the_lipid or the_lipid.lipid.info.level.value < lipid.lipid.info.level.value:
                        the_lipid = lipid
                except:
                    if lipid_name.startswith(UNDEFINED): continue
                    if lipid_name in lipid_class_terms:
                        for lipid_class_term in lipid_class_terms[lipid_name]:
                            for chebi_id in tokens[0].split(" | "):
                                lipid_class_term.relations.add("CHEBI:" + chebi_id)
            if not the_lipid: continue

            for t in tokens[0].split(" | "):
                if t not in chebi_to_lipid: chebi_to_lipid[t] = the_lipid
                elif chebi_to_lipid[t].lipid.info.level.value < the_lipid.lipid.info.level.value:
                    chebi_to_lipid[t] = the_lipid


    chebi_lipid_terms = {}
    ls_id = 1000
    for chebi_id, lipid in chebi_to_lipid.items():
        lipid_name = lipid.get_lipid_string()
        if lipid_name.startswith(UNDEFINED) or lipid_name == "": continue

        lipid_category = lipid.get_lipid_string(LipidLevel.CATEGORY)
        lipid_translates = [lipid.get_lipid_string()]

        for lipid_level in [LipidLevel.COMPLETE_STRUCTURE, LipidLevel.FULL_STRUCTURE, LipidLevel.STRUCTURE_DEFINED, LipidLevel.SN_POSITION, LipidLevel.MOLECULAR_SPECIES, LipidLevel.SPECIES, LipidLevel.CLASS]:
            if lipid.lipid.info.level.value > lipid_level.value:
                lipid_level_name = lipid.get_lipid_string(lipid_level)
                if lipid_level_name != lipid_translates[-1]:
                    lipid_translates.append(lipid_level_name)

        term_list = []
        for lipid_translate in lipid_translates:
            if lipid_translate not in lipid_species_terms:
                if lipid_translate not in chebi_lipid_terms:
                    term_id = f"MOEA:{ls_id:07}"
                    term = Term(term_id, lipid_translate, {"MOEA:0000002"}, _categories = {lipid_category})
                    ls_id += 1
                    chebi_lipid_terms[lipid_translate] = term
                    term_list.append(term)
                else:
                    term_list.append(chebi_lipid_terms[lipid_translate])

            else:
                term_list.append(lipid_species_terms[lipid_translate][-1])

        term_list[0].relations.add("CHEBI:" + chebi_id)
        for i in range(len(term_list) - 1):
            term_list[i].relations.add(term_list[i + 1].get_first_id())

    return ontology_terms, chebi_lipid_terms




if not os.path.exists("Data/LION.pkl.gz"):
    ontology_terms, chebi_lipid_terms = get_lipid_terms()
    with gzip.open("Data/LION.pkl.gz", 'wb') as output_stream:
        pickle.dump([ontology_terms, chebi_lipid_terms], output_stream)
else:
    print("Reading LION pickle")
    with gzip.open("Data/LION.pkl.gz") as input_stream:
        ontology_terms, chebi_lipid_terms = pickle.load(input_stream)



# adding chebi terms with names and relations
chebi_terms_all = {}
with open("Data/chebi.csv", "rt") as infile:
    print("reading chebi")
    for line in infile:
        tokens = line.strip().split("\t")
        chebi_terms_all[tokens[0]] = Term("CHEBI:" + tokens[0], tokens[1], {"MOEA:0000008"} | {"CHEBI:" + t for t in tokens[2:]} if len(tokens) > 2 else {"MOEA:0000008"})


# max_depth = 2
# for chebi_id, c_term in chebi_obo.items():
#     if not chebi_id.startswith("CHEBI:"): continue
#     chebi_id = chebi_id.replace("CHEBI:", "")
#     chebi_term = chebi_terms_all[chebi_id]
#     queue = [(rel, 1) for rel in c_term.relations]
#     while queue:
#         current_id, depth = queue.pop()
#         current_term = chebi_obo[current_id]
#         chebi_term.categories.add(current_term.name.replace("|", "/"))
#         if depth + 1 > max_depth: continue
#         for relation_id in current_term.relations:
#             queue.append((relation_id, depth + 1))

print("set metabolite categories")
chebi_tree = {}
chebi_obo = get_terms("Data/chebi.obo.gz", "CHEBI:")
for chebi_id, c_term in chebi_obo.items():
    chebi_pure_id = chebi_id.split(":")[1]
    if chebi_pure_id not in chebi_terms_all: chebi_terms_all[chebi_pure_id] = Term(chebi_id, c_term.name, {"MOEA:0000008"})
    for relation in c_term.relations:
        if relation not in chebi_tree: chebi_tree[relation] = [chebi_id]
        else: chebi_tree[relation].append(chebi_id)

# for chebi_id, chebi_relations in chebi_tree.items():
#     chebi_pure_id = chebi_id.split(":")[1]
#     if chebi_pure_id not in chebi_terms_all: continue
#     chebi_terms_all[chebi_pure_id].relations |= set(chebi_relations)


for chebi_id, category in CHEBI_CATEGORY_MAPPING.items():
    if chebi_id not in chebi_tree: continue
    queue = [chebi_id]
    visited_terms = set()
    while queue:
        current_chebi_id = queue.pop()
        chebi_pure_id = current_chebi_id.split(":")[1]
        if current_chebi_id in visited_terms: continue
        visited_terms.add(current_chebi_id)

        if chebi_pure_id in chebi_terms_all:
            chebi_terms_all[chebi_pure_id].categories.add(category)
        if current_chebi_id in chebi_tree:
            queue += chebi_tree[current_chebi_id]


for _, chebi_term in chebi_terms_all.items():
    if not chebi_term.categories:
        chebi_term.categories.add("Unclassified metabolite")


with open("Data/rhea_to_chebi.csv", "rt") as infile:
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 2: continue
        rhea_id = "RHEA:" + tokens[0]
        for chebi_term in tokens[1:]:
            if chebi_term in chebi_terms_all:
                chebi_terms_all[chebi_term].relations.add(rhea_id)


chebi_terms_organisms = {org: {} for _, org in species.items()}
pathway_terms = {}
pathway_terms_organisms = {org: {} for _, org in species.items()}
## PATHBANK
# with open("Data/chebi_to_pathway.csv", "rt") as infile:
#     print("Readin chebi_to_pathway")
#     for line in infile:
#         tokens = line.strip().split("\t")
#         if len(tokens) < 3: continue
#
#         if tokens[0] not in chebi_terms_all: continue
#         if tokens[1] not in chebi_terms_organisms: chebi_terms_organisms[tokens[1]] = {}
#         if tokens[0] not in chebi_terms_organisms[tokens[1]]: chebi_terms_organisms[tokens[1]][tokens[0]] = set()
#         chebi_terms_organisms[tokens[1]][tokens[0]] |= set(tokens[2:])

## REACTOME
with gzip.open("Data/ChEBI2Reactome_All_Levels.txt.gz", "rt") as infile:
    print("Readin ChEBI2Reactome")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 6: continue
        organism = tokens[5]
        if organism not in species: continue
        organism = species[organism]
        chebi_id = tokens[0]
        reactome_id = tokens[1]
        reactome_name = tokens[3]
        if chebi_id not in chebi_terms_organisms[organism]: chebi_terms_organisms[organism][chebi_id] = set()
        chebi_terms_organisms[organism][chebi_id].add(reactome_id)
        if reactome_id not in pathway_terms_organisms[organism]:
            term = Term(reactome_id, reactome_name, _namespace = "Metabolic and signalling pathway")
            pathway_terms_organisms[organism][reactome_id] = term
            pathway_terms[reactome_id] = term


# import uniprot accession to gene name
uniprot_data = {}
uniprot_terms_organisms = {organism: {} for organism in species.values()}
uniprot_to_organism = {}
uniprot_to_ec = {}
ensembl_terms_organisms = {organism: {} for ens, organism, ens_cds in ensembl_files}
ensembl_terms = {}
gene_terms = {}
gene_term_organisms = {organism: {} for _, organism, _ in ensembl_files}
with gzip.open("Data/uniprot.csv.gz", "rt") as infile:
    print("Readin uniprot")
    for i, line in enumerate(infile):
        if i == 0: continue
        tokens = line.strip().split("\t")

        if len(tokens) < 2: continue

        uniprot = tokens[0]
        organisms_id = tokens[1]
        if organisms_id not in species_set: continue

        unreviewed = len(tokens) > 10 and tokens[10] == "unreviewed"
        gene_name = tokens[2] if len(tokens) > 2 and len(tokens[2]) > 0 else ""
        protein_name = gene_name if len(gene_name) > 0 else uniprot

        go_terms = {"MOEA:0000006" if unreviewed else "MOEA:0000005"}
        #categories = set(tokens[5].split(", ")) if len(tokens) > 5 and len(tokens[5]) > 0 else set()
        uniprot_term = "UNIPROT:" + uniprot

        function_terms = []

        if len(tokens) > 3 and len(tokens[3]) > 0 and tokens[3].find("[") > -1:
            function_terms.append(tokens[3].lower())
            for entry in tokens[3].split("["):
                go_term = entry.split("]")[0]
                if len(go_term) != 10 or go_term[:3] != "GO:": continue
                go_terms.add(go_term)


        if len(tokens) > 4 and len(tokens[4]) > 0 and tokens[4][0] in EC_CATEGORY:
            uniprot_to_ec[uniprot_term] = EC_CATEGORY[tokens[4][0]]

        gene_ids = set()
        if len(tokens) > 7 and len(tokens[7]) > 1:
            gene_ids.add(tokens[7].strip('";').split(";")[0])

        if len(tokens) > 9 and len(tokens[9]) > 1:
            for ncbi_id in tokens[9].strip(';').split(";"):
                gene_ids.add(f"NCBI:{ncbi_id}")

        if len(tokens) > 11 and len(tokens[11]) > 1: # Keywords
            function_terms.append(tokens[11].lower())

        if len(tokens) > 12 and len(tokens[12]) > 1: # Functions [cc]
            function_terms.append(tokens[12].lower())

        categories = set()
        for cat, patterns in PROTEIN_CATEGORY_RULES.items():
            for pattern in patterns:
                if any(pattern in txt for txt in function_terms):
                    categories.add(cat)

        if not categories: categories = {"Unclassified protein"}

        term = Term(uniprot_term, protein_name, _relations = go_terms | gene_ids, _categories = categories)
        uniprot_terms_organisms[organisms_id][uniprot] = term
        uniprot_data[uniprot] = term
        uniprot_to_organism[uniprot] = tokens[1]

        if len(gene_ids) > 0 and organisms_id in gene_term_organisms:
            # check if at least on gene id is already registered
            common_gene_ids = gene_ids & gene_terms.keys()
            if len(common_gene_ids) > 0:
                common_gene = gene_terms[list(common_gene_ids)[0]]
                common_gene.id |= gene_ids
                for cgi in common_gene_ids | gene_ids:
                    gene_terms[cgi] = common_gene

            else:
                gene_term = Term(gene_ids, f"{tokens[2] if len(tokens) > 2 and len(tokens[2]) > 0 else uniprot} (Gene)", {"MOEA:0000012"})
                for gene_id in gene_ids:
                    gene_terms[gene_id] = gene_term
                    gene_term_organisms[organisms_id][gene_id] = gene_term

        if len(tokens) > 8 and len(tokens[8]) > 0 and organisms_id in gene_term_organisms:
            ensembl_organism = ensembl_terms_organisms[organisms_id]
            for ensembl_sets in tokens[8].strip(";").split("\";\""):
                ensembls = ensembl_sets.split(" ")
                if len(ensembls) < 3: continue
                ENT = ensembls[0].strip('" ;.').split(".")[0]
                ENP = ensembls[1].strip('" ;.').split(".")[0]
                ENG = ensembls[2].strip('" ;.').split(".")[0]


                if ENG in ensembl_organism:
                    ensembl_organism[ENG].relations |= {uniprot_term} | gene_ids
                else:
                    ensembl_term = Term(ENG, gene_name if len(gene_name) > 0 else ENG, {"MOEA:0000011", uniprot_term} | gene_ids)
                    ensembl_organism[ENG] = ensembl_term
                    ensembl_terms[ENG] = ensembl_term

                if ENP in ensembl_organism:
                    ensembl_organism[ENP].relations |= {uniprot_term, ENG}
                else:
                    ensembl_term = Term(ENP, gene_name if len(gene_name) > 0 else ENP, {"MOEA:0000007", uniprot_term, ENG})
                    ensembl_organism[ENP] = ensembl_term
                    ensembl_terms[ENP] = ensembl_term

                if ENT in ensembl_organism:
                    ensembl_organism[ENT].relations |= {ENP, ENG}
                else:
                    ensembl_term = Term(ENT, gene_name if len(gene_name) > 0 else ENT, {"MOEA:0000010", ENP, ENG})
                    ensembl_organism[ENT] = ensembl_term
                    ensembl_terms[ENT] = ensembl_term



for ensembl_file, organism, ensembl_cds in ensembl_files:
    ensembl_organism = ensembl_terms_organisms[organism]
    print(f"Readin Ensembl {organism}")
    try:
        with gzip.open(ensembl_file, "rt") as infile:
            for i, line in enumerate(infile):
                if i == 0: continue
                tokens = line.split("\t")
                if len(tokens) < 4: continue
                uniprot = tokens[3].split("-")[0]
                ENG, ENT, ENP = tokens[:3]

                uniprot_term = "UNIPROT:" + uniprot
                uniprot_known = uniprot in uniprot_data

                gene_name = uniprot_data[uniprot].name if (uniprot in uniprot_data and list(uniprot_data[uniprot].id)[0].find(uniprot_data[uniprot].name) < 0) else ""

                if ENG not in ensembl_organism:
                    ensembl_organism[ENG] = Term(ENG, gene_name if len(gene_name) > 0 else ENG, {"MOEA:00000011"} | ({uniprot_term} if uniprot_known else set()))
                elif uniprot_known:
                    ensembl_organism[ENG].relations.add(uniprot_term)

                if ENP != ENG:
                    if ENP not in ensembl_organism:
                        ensembl_organism[ENP] = Term(ENP, gene_name if len(gene_name) > 0 else ENP, {"MOEA:0000007", ENG} | ({uniprot_term} if uniprot_known else set()))
                    else:
                        ensembl_organism[ENP].relations |= {ENG} | ({uniprot_term} if uniprot_known else set())

                if ENT not in ensembl_organism:
                    ensembl_organism[ENT] = Term(ENT, gene_name if len(gene_name) > 0 else ENT, {"MOEA:0000010", ENG, ENP})
                else:
                    ensembl_organism[ENT].relations |= {ENG, ENP}



        with gzip.open(ensembl_file.replace(".uniprot.", ".all.")) as infile:
            for i, line in enumerate(infile):
                if i == 0: continue
                tokens = line.decode("utf8").split("\t")
                if len(tokens) < 4: continue

                ENG, ENT = tokens[2], tokens[3]
                if ENG not in ensembl_organism: continue

                ENT_term = None
                if ENT not in ensembl_organism:
                    ENT_term = Term(ENT, ENT, {"MOEA:0000010", ENG})
                    ensembl_organism[ENT] = ENT_term

                if len(tokens) > 4 and len(tokens[4]) > 0:
                    ENP = tokens[4]
                    if ENP != ENG and ENP not in ensembl_organism:
                        ensembl_organism[ENP] = Term(ENP, ENP, {"MOEA:0000007", ENG})
                        if ENT_term != None:
                            ENT_term.relations.add(ENP)

    except Exception as e:
        print(e)
        exit(-1)

    try:
        with gzip.open(ensembl_cds, "rt") as infile:
            for inline in infile:
                if inline[0] != ">": continue
                tokens = inline[1:].strip().split(" ")
                ensembl_prot = tokens[0].split(".")[0]
                ensembl_gene = tokens[3][5:].split(".")[0]
                biotype = gene_biotype_dict[bt] if (bt := tokens[4][13:]) in gene_biotype_dict else "Unclassified gene"
                if ensembl_prot in ensembl_terms: ensembl_terms[ensembl_prot].categories.add(biotype)
                if ensembl_gene in ensembl_terms: ensembl_terms[ensembl_gene].categories.add(biotype)

                hgnc_pos = inline.find("Acc:HGNC:")
                if hgnc_pos  == -1: continue
                hgnc_id = inline[hgnc_pos + 4 : ].split("]")[0]
                if ensembl_prot in ensembl_terms: ensembl_terms[ensembl_prot].relations.add(hgnc_id)
                if ensembl_gene in ensembl_terms: ensembl_terms[ensembl_gene].relations.add(hgnc_id)

    except Exception as e:
        print(e)
        exit(-1)


with gzip.open("Data/goa_uniprot.csv.gz", "rt") as input_stream:
    for line in input_stream:
        tokens = line.strip().split("\t")
        if len(tokens) < 2: continue
        if tokens[0] in uniprot_data: uniprot_data[tokens[0]].relations |= set(tokens[1:]) - {""}


## PATHBANK
# with open("Data/pathbank_to_uniprot.csv", "rt") as infile:
#     print("Readin pathbank_to_uniprot")
#     for line in infile:
#         tokens = line.strip().split("\t")
#         if len(tokens) < 3: continue
#
#         i = 2
#         organism = None
#         while i < len(tokens):
#             if tokens[i] in uniprot_to_organism:
#                 organism = uniprot_to_organism[tokens[i]]
#                 break
#             i += 1
#
#         if organism == None: continue
#         if organism not in pathway_terms_organisms: pathway_terms_organisms[organism] = {}
#         pathway_terms_organisms[organism][tokens[0]] = Term(tokens[0], tokens[1], _namespace = "Metabolic and signalling pathway")
#         for term_id in tokens[2:]:
#             if term_id in uniprot_data:
#                 uniprot_data[term_id].relations.add(tokens[0])

## REACTOME
with gzip.open("Data/UniProt2Reactome_All_Levels.txt.gz", "rt") as infile:
    print("Readin UniProt2Reactome")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 6: continue
        organism = tokens[5]
        uniprot = tokens[0]
        if organism not in species or uniprot not in uniprot_data: continue
        organism = species[organism]
        reactome_id = tokens[1]
        reactome_name = tokens[3]
        uniprot_data[uniprot].relations.add(reactome_id)
        if reactome_id not in pathway_terms_organisms[organism]:
            term = Term(reactome_id, reactome_name, _namespace = "Metabolic and signalling pathway")
            pathway_terms_organisms[organism][reactome_id] = term
            pathway_terms[reactome_id] = term


with gzip.open("Data/Ensembl2Reactome_All_Levels.txt.gz", "rt") as infile:
    print("Readin Ensembl2Reactome")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 6: continue
        organism = tokens[5]
        if organism not in species: continue
        organism = species[organism]
        ensembl = tokens[0]
        if ensembl not in ensembl_terms_organisms[organism]: continue
        reactome_id = tokens[1]
        reactome_name = tokens[3]
        ensembl_terms_organisms[organism][ensembl].relations.add(reactome_id)
        if reactome_id not in pathway_terms_organisms[organism]:
            term = Term(reactome_id, reactome_name, _namespace = "Metabolic and signalling pathway")
            pathway_terms_organisms[organism][reactome_id] = term
            pathway_terms[reactome_id] = term



with open("Data/ReactomePathwaysRelation.txt") as infile:
    print("Readin ReactomePathwaysRelation")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 2 or tokens[0] not in pathway_terms or tokens[1] not in pathway_terms: continue
        pathway_terms[tokens[1]].relations.add(tokens[0])

    #
    # for rhea_term_id, uniprot_terms in rhea_terms_organisms[tax_id].items():
    #     rhea_terms[rhea_term_id] = Term(
    #         f"RHEA:{rhea_term_id}",
    #         f"Rhea reaction {rhea_term_id}",
    #         set(["UNIPROT:" + uniprot for uniprot in uniprot_terms])
    #     )


rhea_terms_organisms = {tax_id: {} for tax_id in uniprot_terms_organisms}
with open("Data/rhea2uniprot_sprot.tsv", "rt") as infile:
    print("Readin rhea2 sprot")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) != 4: continue

        if tokens[3] not in uniprot_to_organism: continue
        organism = uniprot_to_organism[tokens[3]]

        if tokens[0] not in rhea_terms_organisms[organism]:
            rhea_terms_organisms[organism][tokens[0]] = Term(
                f"RHEA:{tokens[0]}",
                f"Rhea reaction {tokens[0]}",
                {"UNIPROT:" + tokens[3], "MOEA:0000009"}
            )
        else:
            rhea_terms_organisms[organism][tokens[0]].relations.add("UNIPROT:" + tokens[3])

        if tokens[2] not in rhea_terms_organisms[organism]:
            rhea_terms_organisms[organism][tokens[2]] = Term(
                f"RHEA:{tokens[2]}",
                f"Rhea reaction {tokens[2]}",
                {"UNIPROT:" + tokens[3], "MOEA:0000009"}
            )
        else:
            rhea_terms_organisms[organism][tokens[2]].relations.add("UNIPROT:" + tokens[3])



for _, rhea_terms in rhea_terms_organisms.items():
    for _, rhea_term in rhea_terms.items():
        rhea_term_categories = {uniprot_to_ec[relation] for relation in rhea_term.relations if relation in uniprot_to_ec}
        if len(rhea_term_categories) == 0: rhea_term_categories = {"Unclassified reaction"}
        rhea_term.categories |= rhea_term_categories



with gzip.open("Data/rhea2uniprot_trembl.tsv.gz", "rt") as infile:
    print("Readin rhea2 trembl")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) != 4: continue

        if tokens[3] not in uniprot_to_organism: continue
        organism = uniprot_to_organism[tokens[3]]

        if tokens[0] not in rhea_terms_organisms[organism]:
            rhea_terms_organisms[organism][tokens[0]] = Term(
                f"RHEA:{tokens[0]}",
                f"Rhea reaction {tokens[0]}",
                {"UNIPROT:" + tokens[3], "MOEA:0000009"}
            )
        else:
            rhea_terms_organisms[organism][tokens[0]].relations.add("UNIPROT:" + tokens[3])

        if tokens[2] not in rhea_terms_organisms[organism]:
            rhea_terms_organisms[organism][tokens[2]] = Term(
                f"RHEA:{tokens[2]}",
                f"Rhea reaction {tokens[2]}",
                {"UNIPROT:" + tokens[3], "MOEA:0000009"}
            )
        else:
            rhea_terms_organisms[organism][tokens[2]].relations.add("UNIPROT:" + tokens[3])








print("Readin reactome reactions")
with open("Data/ChEBI2ReactomeReactions.txt", "rt") as input_stream:
    for line in input_stream:
        tokens = line.strip().split("\t")

        if tokens[0] in chebi_terms_all: chebi_terms_all[tokens[0]].relations.add(tokens[1])


reactome_reactions = {}
with open("Data/UniProt2ReactomeReactions.txt", "rt") as input_stream:
    for line in input_stream:
        tokens = line.strip().split("\t")

        uniprot_id = "UNIPROT:" + tokens[0]
        if tokens[1] not in reactome_reactions: reactome_reactions[tokens[1]] = Term(tokens[1], f"Reactome reaction {tokens[1]}", {uniprot_id, "MOEA:0000009"})
        else: reactome_reactions[tokens[1]].relations.add(uniprot_id)









disease_and_phenotype_terms = get_terms("Data/doid-base.obo", "DOID:", True)
for _, doid_term in disease_and_phenotype_terms.items():
    doid_term.namespace.add("Disease")
for k, v in {tokens[0]: set(tokens[1:]) for line in open("Data/uniprot_to_doid.csv").read().split("\n") if len(line) > 0 and (tokens := line.split("\t")) and len(tokens[1]) > 0}.items():
    if k in uniprot_data:
        uniprot_data[k].relations |= v



mondo_terms = get_terms("Data/mondo.obo", "MONDO:", True)
omim_to_mondo = {}
orphanet_to_mondo = {}
hpo_to_mondo = {}
for mondo_id, mondo_term in mondo_terms.items():
    mondo_term.namespace.add("Disease")
    found_xref = False
    for xref_id in mondo_term.xref:
        if xref_id.startswith("OMIM:"): omim_to_mondo[xref_id] = mondo_id
        if xref_id.startswith("Orphanet:"): orphanet_to_mondo[xref_id] = mondo_id
        if xref_id.startswith("HP:"): hpo_to_mondo[xref_id] = mondo_id
        if xref_id in disease_and_phenotype_terms:
            found_xref = True
            disease_and_phenotype_terms[xref_id].merge(mondo_term)
            disease_and_phenotype_terms[mondo_id] = disease_and_phenotype_terms[xref_id]

    if not found_xref:
        disease_and_phenotype_terms[mondo_id] = mondo_term



hpo_terms = get_terms("Data/hp.obo", "HP:", True)
for hpo_term_id, hpo_term in hpo_terms.items():
    hpo_term.namespace.add("Phenotype")
    if hpo_term_id in hpo_to_mondo and hpo_to_mondo[hpo_term_id] in disease_and_phenotype_terms:
        disease_and_phenotype_terms[hpo_to_mondo[hpo_term_id]].merge(hpo_term)
        disease_and_phenotype_terms[hpo_term_id] = disease_and_phenotype_terms[hpo_to_mondo[hpo_term_id]]
    else:
        disease_and_phenotype_terms[hpo_term_id] = hpo_term



print("Readin genes_to_phenotype")
for line in open("Data/genes_to_phenotype.txt").read().split("\n"):
    if len(line) < 2: continue
    tokens = line.split("\t")
    if len(tokens) < 3 or tokens[1] == "-": continue
    ncbi_id, ncbi_name, hpo_id, orphanet_id = tokens[0], tokens[1], tokens[2], tokens[5].replace("ORPHA:", "Orphanet:")
    ncbi_id = "NCBI:" + ncbi_id
    gene_name = ncbi_name + " (Gene)"
    if ncbi_id not in gene_terms:
        term = Term(ncbi_id, gene_name, {hpo_id} | {"MOEA:0000012"})
        gene_terms[ncbi_id] = term
        gene_term_organisms["9606"][ncbi_id] = term
    else:
        gene_terms[ncbi_id].name = gene_name
        gene_terms[ncbi_id].relations.add(hpo_id)
    if orphanet_id in orphanet_to_mondo: gene_terms[ncbi_id].relations.add(orphanet_to_mondo[orphanet_id])




print("Readin genes_to_disease")
for line in open("Data/genes_to_disease.txt").read().split("\n"):
    if len(line) < 2: continue
    tokens = line.split("\t")
    if len(tokens) < 3 or tokens[1] == "-": continue
    ncbi_id, ncbi_name, omim_id = tokens[0], tokens[1], tokens[3]
    if omim_id not in omim_to_mondo: continue
    ncbi_id = ncbi_id.replace("NCBIGene:", "NCBI:")
    gene_name = ncbi_name + " (Gene)"
    if ncbi_id not in gene_terms:
        term = Term(ncbi_id, gene_name, {omim_to_mondo[omim_id]} | {"MOEA:0000012"})
        gene_terms[ncbi_id] = term
        gene_term_organisms["9606"][ncbi_id] = term

    else:
        gene_terms[ncbi_id].name = gene_name
        gene_terms[ncbi_id].relations.add(omim_to_mondo[omim_id])



with open("Data/HGNC.csv") as infile:
    ensembl_organism = ensembl_terms_organisms["9606"]
    print("Readin HGNC")
    for i, line in enumerate(infile):
        if i == 0: continue
        tokens = line.strip().split("\t")
        if len(tokens) < 2: continue
        hgnc_id = tokens[0]
        hgnc_name = tokens[1] + " (Gene)"
        hgnc_term = Term(hgnc_id, hgnc_name, {"MOEA:0000012"}) if hgnc_id not in gene_terms else gene_terms[hgnc_id]
        gene_ids = {hgnc_id}

        uniprot = tokens[2] if len(tokens) > 2 and len(tokens[2]) > 0 else None
        ncbi_id = "NCBI:" + tokens[3] if len(tokens) > 3 and len(tokens[3]) > 0 else None
        ENG = tokens[4] if len(tokens) > 4 and len(tokens[4]) > 0 else None

        if uniprot != None and uniprot in uniprot_data:
            uniprot_data[uniprot].relations.add(hgnc_id)

        if ncbi_id != None:
            gene_ids.add(ncbi_id)
            if ncbi_id in gene_terms:
                gene_terms[ncbi_id].merge(hgnc_term)
                if hgnc_id not in gene_terms:
                    gene_terms[hgnc_id] = gene_terms[ncbi_id]
                    gene_term_organisms["9606"][hgnc_id] = gene_terms[ncbi_id]
                else:
                    gene_terms[hgnc_id].merge(gene_terms[ncbi_id])

            else:
                gene_terms[hgnc_id] = hgnc_term
                gene_term_organisms["9606"][hgnc_id] = hgnc_term
                hgnc_term.id.add(ncbi_id)
                gene_terms[ncbi_id] = hgnc_term
                gene_term_organisms["9606"][ncbi_id] = hgnc_term

        if ENG != None:
            if ENG in ensembl_organism:
                ensembl_organism[ENG].relations |= gene_ids

            else:
                ensembl_organism[ENG] = Term(ENG, ENG, {"MOEA:0000011"} | gene_ids)


print("Readin hgnc_to_mondo")
for line in open("Data/hgnc_to_mondo.csv").read().split("\n"):
    if len(line) < 1: continue
    tokens = line.split("\t")
    if len(tokens) < 2: continue
    hgnc_id, mondo_ids = tokens[0], set(tokens[1:])
    if hgnc_id in gene_terms: gene_terms[hgnc_id].relations |= mondo_ids



for ensembl_id, ensembl_term in ensembl_terms.items():
    for relation_id in ensembl_term.relations:
        if relation_id in gene_terms:
            gene_terms[relation_id].categories |= ensembl_term.categories


for _, gene_term in gene_terms.items():
    if not gene_term.categories:
        gene_term.categories.add("Unclassified gene")



go_terms = get_terms("Data/go-basic.obo", "GO:", upper = True, with_relationship = True)
namespaces = {
    "biological_process": "Biological process",
    "cellular_component": "Cellular component",
    "molecular_function": "Molecular function"
}
for go_term_id, go_term in go_terms.items():
    go_term.namespace = {namespaces[n] if n in namespaces else n for n in go_term.namespace}


for uniprot_term in uniprot_data.values():
    if "Enzyme" in uniprot_term.categories:
        if "MOEA:0000005" in uniprot_term.relations:
            uniprot_term.namespace.add("Enzymatic activity (Swiss-Prot)")
        uniprot_term.namespace.add("Enzymatic activity (Swiss-Prot + TrEMBL)")



# do several organisms
for tax_name, tax_id in species.items():
    print(tax_id)
    output = []
    for go_id, term in go_terms.items():
        output.append(term.to_string())

    # adding all remaining categories for multiomics ontology
    output.append(Term("MOEA:0000004", "Unspecific lipid").to_string())
    output.append(Term("MOEA:0000005", "Uniprot protein swissprot").to_string())
    output.append(Term("MOEA:0000006", "Uniprot protein trembl").to_string())
    output.append(Term("MOEA:0000007", "Ensembl protein").to_string())
    output.append(Term("MOEA:0000008", "ChEBI metabolite").to_string())
    output.append(Term("MOEA:0000009", "Generic reaction").to_string())
    output.append(Term("MOEA:0000010", "Ensembl transcript").to_string())
    output.append(Term("MOEA:0000011", "Ensembl gene").to_string())
    output.append(Term("MOEA:0000012", "Gene").to_string())

    reactome_tag = reactomes[tax_id]
    chebi_terms_organism = chebi_terms_organisms[tax_id]

    # # KEGG
    # try:
    #     kegg_to_uniprot = {"KEGG:" + t[0]: (t[1], set(t[2:])) for line in open(f"Data/kegg_to_uniprot_{reactome_tag.lower()}.csv").read().split("\n") if (t := line.split("\t")) and len(t) > 1}
    #     kegg_to_chebi = {"KEGG:" + t[0]: (t[1], set(t[2:])) for line in open(f"Data/kegg_to_chebi_{reactome_tag.lower()}.csv").read().split("\n") if (t := line.split("\t")) and len(t) > 1}
    #
    #
    #     for kegg_id, (kegg_pathway_name, chebi_ids) in kegg_to_chebi.items():
    #         for chebi_id in chebi_ids:
    #             chebi_id_clean = chebi_id[6:]
    #             if chebi_id_clean not in chebi_terms_organism: chebi_terms_organism[chebi_id_clean] = {kegg_id}
    #             else: chebi_terms_organism[chebi_id_clean].add(kegg_id)
    #
    #     for kegg_id, (kegg_pathway_name, uniprot_ids) in kegg_to_uniprot.items():
    #         output.append(Term(kegg_id, kegg_pathway_name, _namespace = "KEGG").to_string())
    #         for uniprot_id in uniprot_ids:
    #             if uniprot_id in uniprot_data:
    #                 uniprot_data[uniprot_id].relations.add(kegg_id)
    # except Exception as e:
    #     print(e)


    for i, row in df_dolphin[df_dolphin["protein_UniProt_ID"].isin(list(uniprot_terms_organisms[tax_id].keys()))].iterrows():
        output.append(Term("BD:" + row["BioDolphinID"], f"Interaction {row["BioDolphinID"]}", {"UNIPROT:" + row["protein_UniProt_ID"], "MOEA:0000009"}).to_string())
        chebi_id = row["lipid_ChEBI"][6:]
        if chebi_id not in chebi_terms_organism: chebi_terms_organism[chebi_id] = {"BD:" + row["BioDolphinID"]}
        else: chebi_terms_organism[chebi_id].add("BD:" + row["BioDolphinID"])


    if reactome_tag:
        for reactome_id, reactome_term in reactome_reactions.items():
            if reactome_id.find(reactome_tag) > -1:
                output.append(reactome_term.to_string())



    chebi_terms = {t: c.copy() for t, c in chebi_terms_all.items()}
    for chebi_term, pathway_terms in chebi_terms_organism.items():
        if chebi_term not in chebi_terms: continue
        chebi_terms[chebi_term].relations |= pathway_terms


    if tax_id in pathway_terms_organisms:
        pathway_terms = pathway_terms_organisms[tax_id]
        for pathway_id, pathway_term in pathway_terms.items():
            output.append(pathway_term.to_string())


    protein_terms = uniprot_terms_organisms[tax_id]
    for _, protein_term in protein_terms.items():
        output.append(protein_term.to_string())


    if tax_id in ensembl_terms_organisms:
        ensembl_terms = ensembl_terms_organisms[tax_id]
        for _, ensembl_term in ensembl_terms.items():
            output.append(ensembl_term.to_string())


    if tax_id in gene_term_organisms:
        written_genes = set()
        for _, gene_term in gene_term_organisms[tax_id].items():
            if gene_term in written_genes: continue
            written_genes.add(gene_term)
            output.append(gene_term.to_string())


    # create organism-specific rhea terms
    if tax_id in rhea_terms_organisms:
        for _, rhea_term in rhea_terms_organisms[tax_id].items():
            output.append(rhea_term.to_string())


    for _, chebi_term in chebi_terms.items():
        output.append(chebi_term.to_string())

    # diseases and phenotypes
    written_diseases = set()
    for _, disease_term in disease_and_phenotype_terms.items():
        if disease_term in written_diseases: continue
        written_diseases.add(disease_term)
        output.append(disease_term.to_string())


    for lipid_name, term in chebi_lipid_terms.items():
        output.append(term.to_string())


    for lipid_maps_id, term in lipid_maps_terms.items():
        output.append(term.to_string())


    for term_id, term in ontology_terms.items():
        output.append(term.to_string())


    with gzip.open(f"../ontology_{tax_id}.gz", "wb") as gz_output:
        gzip_out = "".join(output).encode("utf8")
        print("writing", len(gzip_out))
        gz_output.write(gzip_out)















