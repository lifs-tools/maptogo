from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidLevel import LipidLevel
import gzip
import sys

all_species = len(sys.argv) > 1 and sys.argv[1] == "all"
parser = LipidParser()

if all_species:
    print("Creating all ontologies")
    species = {
        'Homo sapiens': "9606",
        'Mus musculus': "10090",
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
        ("Data/Homo_sapiens.uniprot.tsv.gz", "9606"),
        ("Data/Mus_musculus.uniprot.tsv.gz", "10090"),
        ("Data/Saccharomyces_cerevisiae.uniprot.tsv.gz", "4932"),
        ("Data/Drosophila_melanogaster.uniprot.tsv.gz", "7227"),
        ("Data/Rattus_norvegicus.uniprot.tsv.gz", "10116"),
        ("Data/Bos_taurus.uniprot.tsv.gz", "9913"),
        ("Data/Caenorhabditis_elegans.uniprot.tsv.gz", "6239"),
    ]

else:
    print("Creating Human & Mouse ontology")
    species = {
        'Homo sapiens': "9606",
        'Mus musculus': "10090",
    }
    ensembl_files = [
        ("Data/Homo_sapiens.uniprot.tsv.gz", "9606"),
        ("Data/Mus_musculus.uniprot.tsv.gz", "10090"),
    ]

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
        

    def to_string(self):
        return f"{'|'.join(self.id)}\t{self.name}\t{'|'.join(self.relations)}\t{'|'.join(self.synonyms)}\t{'|'.join(self.namespace)}\t{'|'.join(self.categories)}\n"


    def copy(self):
        return Term(set(self.id), self.name, set(self.relations), set(self.synonyms), set(self.namespace), set(self.xref), _categories = set(self.categories))


    def merge(self, copy):
        self.id |= copy.id
        self.relations |= copy.relations
        self.synonyms |= copy.synonyms
        self.namespace |= copy.namespace
        self.xref |= copy.xref
        self.categories |= copy.categories



def get_terms(term_file, prefix, upper = False):
    ontology_terms = {}
    with open(term_file, "rt") as obo:
        term_id, name, relations, synonyms, namespace, is_obsolete, xref = "", "", set(), [], "", False, set()
        for i, line in enumerate(obo):
            line = line.strip()
            if len(line) == 0: continue
            if i > 0 and i % 10000 == 0: print(term_file, i)

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
    # import chebi to lipid dictionary
    chebi_to_lipid = {t.split("\t")[0]: {t.split("\t")[1]} for t in open("Data/chebi_to_lipid.csv", "rt").read().split("\n") if t.find("\t") > -1}

    for directory in ["ChEBI", "LipidMaps", "SwissLipids"]:
        for line in open(f"Data/{directory}_chebi_to_lipid.csv", "rt").read().split("\n"):
            if len(line) < 1: continue
            tokens = line.split("\t")
            if len(tokens) < 2: continue

            for t in tokens[0].split(" | "):
                if t not in chebi_to_lipid: chebi_to_lipid[t] = set()
                chebi_to_lipid[t] |= set(tokens[1:])



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
                            relations.add("LS:0000002")
                    except Exception as e:
                        pass


                elif line.startswith("is_a:"):
                    relation = line[5:].split(" ! ")[0].strip(" ")
                    relations.add(relation)
                    if relation == "LS:0000001": is_class = True
                    elif relation == "LS:0000002": is_species = True

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



    chebi_lipid_terms = {}
    ls_id = 1000
    UNDEFINED = "UNDEFINED"
    for chebi_id, lipid_names in chebi_to_lipid.items():
        for lipid_name in lipid_names:
            if lipid_name[:len(UNDEFINED)] == UNDEFINED or lipid_name == "": continue
            if lipid_name in lipid_class_terms:
                for lipid_class_term in lipid_class_terms[lipid_name]:
                    lipid_class_term.relations.add("CHEBI:" + chebi_id)

                continue

            lipid_translates = set()
            try:
                lipid = parser.parse(lipid_name)
                lipid_category = lipid.get_lipid_string(LipidLevel.CATEGORY)
                lipid_translates.add((lipid.get_lipid_string(), lipid_category))
                if lipid.lipid.info.level.value >= LipidLevel.COMPLETE_STRUCTURE.value:
                    lipid_translates.add((lipid.get_lipid_string(LipidLevel.COMPLETE_STRUCTURE), lipid_category))
                if lipid.lipid.info.level.value >= LipidLevel.FULL_STRUCTURE.value:
                    lipid_translates.add((lipid.get_lipid_string(LipidLevel.FULL_STRUCTURE), lipid_category))
                if lipid.lipid.info.level.value >= LipidLevel.STRUCTURE_DEFINED.value:
                    lipid_translates.add((lipid.get_lipid_string(LipidLevel.STRUCTURE_DEFINED), lipid_category))
                if lipid.lipid.info.level.value >= LipidLevel.SN_POSITION.value:
                    lipid_translates.add((lipid.get_lipid_string(LipidLevel.SN_POSITION), lipid_category))
                if lipid.lipid.info.level.value >= LipidLevel.MOLECULAR_SPECIES.value:
                    lipid.lipid.sort_fatty_acyl_chains()
                    lipid_translates.add((lipid.get_lipid_string(LipidLevel.MOLECULAR_SPECIES), lipid_category))
                if lipid.lipid.info.level.value >= LipidLevel.SPECIES.value:
                    lipid_translates.add((lipid.get_lipid_string(LipidLevel.SPECIES), lipid_category))

            except Exception as e:
                lipid_translates.add((lipid_name, ""))

            for lipid_translate, lipid_category in lipid_translates:
                if lipid_translate not in lipid_species_terms:
                    if lipid_translate not in chebi_lipid_terms:
                        term = Term(f"LS:{ls_id:07}", lipid_translate, {"LS:0000002"}, _categories = {lipid_category})
                        term.relations.add("CHEBI:" + chebi_id)
                        ls_id += 1
                        chebi_lipid_terms[lipid_translate] = term
                    else:
                        chebi_lipid_terms[lipid_translate].relations.add("CHEBI:" + chebi_id)

                else:
                    for lipid_species_term in lipid_species_terms[lipid_translate]:
                        lipid_species_term.relations.add("CHEBI:" + chebi_id)

    return ontology_terms, chebi_lipid_terms



# adding chebi terms with names and relations
chebi_terms_all = {}
with open("Data/chebi.csv", "rt") as infile:
    print("reading chebi")
    for line in infile:
        tokens = line.strip().split("\t")
        chebi_terms_all[tokens[0]] = Term("CHEBI:" + tokens[0], tokens[1], set(["CHEBI:" + t for t in tokens[2:]]) if len(tokens) > 2 else set())
        chebi_terms_all[tokens[0]].relations.add("LS:0000005")


chebi_obo = get_terms("Data/chebi_lite.obo", "CHEBI:")
max_depth = 2
for chebi_id, c_term in chebi_obo.items():
    if not chebi_id.startswith("CHEBI:"): continue
    chebi_id = chebi_id.replace("CHEBI:", "")
    chebi_term = chebi_terms_all[chebi_id]
    queue = [(rel, 1) for rel in c_term.relations]
    while queue:
        current_id, depth = queue.pop()
        current_term = chebi_obo[current_id]
        chebi_term.categories.add(current_term.name.replace("|", "/"))
        if depth + 1 > max_depth: continue
        for relation_id in current_term.relations:
            queue.append((relation_id, depth + 1))


with open("Data/rhea_to_chebi.csv", "rt") as infile:
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 2: continue
        rhea_id = "RHEA:" + tokens[0]
        for chebi_term in tokens[1:]:
            if chebi_term in chebi_terms_all:
                chebi_terms_all[chebi_term].relations.add(rhea_id)


chebi_terms_organisms = {}
with open("Data/chebi_to_pathway.csv", "rt") as infile:
    print("Readin chebi_to_pathway")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 3: continue

        if tokens[0] not in chebi_terms_all: continue
        if tokens[1] not in chebi_terms_organisms: chebi_terms_organisms[tokens[1]] = {}
        if tokens[0] not in chebi_terms_organisms[tokens[1]]: chebi_terms_organisms[tokens[1]][tokens[0]] = set()
        chebi_terms_organisms[tokens[1]][tokens[0]] |= set(tokens[2:])


# import uniprot accession to gene name
uniprot_data = {}
uniprot_terms_organisms = {}
uniprot_to_organism = {}
ensembl_terms_organisms = {organism: {} for _, organism in ensembl_files}
with gzip.open("Data/uniprot.csv.gz", "rt") as infile:
    print("Readin uniprot")
    for i, line in enumerate(infile):
        if i == 0: continue
        tokens = line.split("\t")

        if len(tokens) < 2: continue

        uniprot = tokens[0]
        organisms_id = tokens[1]
        protein_name = tokens[2] if len(tokens) > 2 and len(tokens[2]) > 0 else f"{uniprot} (Protein)"
        go_terms = {"LS:0000004"}
        categories = set(tokens[4].split(", ")) if len(tokens) > 4 and len(tokens[4]) > 0 else set()
        uniprot_term = "UNIPROT:" + uniprot

        if len(tokens) > 3 and len(tokens[3]) > 0 and tokens[3].find("[") > -1:
            for entry in tokens[3].split("["):
                go_term = entry.split("]")[0]
                if len(go_term) != 10 or go_term[:3] != "GO:": continue
                go_terms.add(go_term)

        term = Term(uniprot_term, protein_name, go_terms, _categories = categories)
        if organisms_id not in uniprot_terms_organisms: uniprot_terms_organisms[organisms_id] = {}
        uniprot_terms_organisms[organisms_id][uniprot] = term
        uniprot_data[uniprot] = term
        uniprot_to_organism[uniprot] = tokens[1]

        if len(tokens) < 8 or len(tokens[7]) == 0: continue
        ensembls = set()
        for cc in tokens[7].replace("\"", "").split(";"):
            ensembls |= set([tok.split(".")[0] for t in cc.split(" ") if ((tok := t.strip("., ")) and tok[:2] == "EN")])

        for ensembl in ensembls:
            if organisms_id not in ensembl_terms_organisms: ensembl_terms_organisms[organisms_id] = {}
            if ensembl not in ensembl_terms_organisms[organisms_id]:
                term = Term(ensembl, ensembl, {"LS:0000006"}, _categories = set(uniprot_data[tokens[0]].categories))
                ensembl_terms_organisms[organisms_id][ensembl] = term
            ensembl_terms_organisms[organisms_id][ensembl].relations.add(uniprot_term)


for ensembl_file, organism in ensembl_files:
    ensembl_organism = ensembl_terms_organisms[organism]
    print(f"Readin Ensembl {organism}")
    try:
        with gzip.open(ensembl_file) as infile:
            for i, line in enumerate(infile):
                if i == 0: continue
                tokens = line.decode("utf8").split("\t")
                if len(tokens) < 4: continue
                uniprot = tokens[3].split("-")[0]
                ENG, ENT, ENP = tokens[:3]

                relations = {"LS:0000006", ENG}
                uniprot_known = uniprot in uniprot_data
                if uniprot_known: relations.add(f"UNIPROT:{uniprot}")

                if ENG not in ensembl_organism:
                    ensembl_organism[ENG] = Term(ENG, ENG, {"LS:0000006"})
                    if uniprot_known: uniprot_data[uniprot].relations.add(ENG)

                if ENT not in ensembl_organism:
                    ensembl_organism[ENT] = Term(ENT, ENT, relations)

                if ENP != ENG and ENP not in ensembl_organism:
                    ensembl_organism[ENP] = Term(ENP, ENP, relations)


        with gzip.open(ensembl_file.replace(".uniprot.", ".all.")) as infile:
            for i, line in enumerate(infile):
                if i == 0: continue
                tokens = line.decode("utf8").split("\t")
                if len(tokens) < 4: continue

                ENG, ENT = tokens[2], tokens[3]
                if ENG not in ensembl_organism: continue

                if ENT not in ensembl_organism:
                    ensembl_organism[ENT] = Term(ENT, ENT, {"LS:0000006", ENG})

                if len(tokens) > 4 and len(tokens[4]) > 0:
                    ENP = tokens[4]
                    if ENP != ENG and ENP not in ensembl_organism:
                        ensembl_organism[ENP] = Term(ENP, ENP, {"LS:0000006", ENG})

    except Exception as e:
        print(e)
        exit(-1)


pathway_terms_organisms = {}
with open("Data/pathbank_to_uniprot.csv", "rt") as infile:
    print("Readin pathbank_to_uniprot")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) < 3: continue

        i = 2
        organism = None
        while i < len(tokens):
            if tokens[i] in uniprot_to_organism:
                organism = uniprot_to_organism[tokens[i]]
                break
            i += 1

        if organism == None: continue
        if organism not in pathway_terms_organisms: pathway_terms_organisms[organism] = {}
        pathway_terms_organisms[organism][tokens[0]] = Term(tokens[0], tokens[1], _namespace = "Metabolic and signalling pathway")
        for term_id in tokens[2:]:
            if term_id in uniprot_data:
                uniprot_data[term_id].relations.add(tokens[0])


rhea_terms_organisms = {tax_id: {} for tax_id in uniprot_terms_organisms}
with open("Data/rhea2uniprot_sprot.tsv", "rt") as infile:
    print("Readin rhea2 sprot")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) != 4: continue

        if tokens[3] not in uniprot_to_organism: continue
        organism = uniprot_to_organism[tokens[3]]

        if tokens[0] not in rhea_terms_organisms[organism]: rhea_terms_organisms[organism][tokens[0]] = set()
        rhea_terms_organisms[organism][tokens[0]].add(tokens[3])
        if tokens[2] not in rhea_terms_organisms[organism]: rhea_terms_organisms[organism][tokens[2]] = set()
        rhea_terms_organisms[organism][tokens[2]].add(tokens[3])


with gzip.open("Data/rhea2uniprot_trembl.tsv.gz", "rt") as infile:
    print("Readin rhea2 trembl")
    for line in infile:
        tokens = line.strip().split("\t")
        if len(tokens) != 4: continue

        if tokens[3] not in uniprot_to_organism: continue
        organism = uniprot_to_organism[tokens[3]]

        if tokens[0] not in rhea_terms_organisms[organism]: rhea_terms_organisms[organism][tokens[0]] = set()
        rhea_terms_organisms[organism][tokens[0]].add(tokens[3])
        if tokens[2] not in rhea_terms_organisms[organism]: rhea_terms_organisms[organism][tokens[2]] = set()
        rhea_terms_organisms[organism][tokens[2]].add(tokens[3])




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



gene_terms = {}
print("Readin genes_to_phenotype")
for line in open("Data/genes_to_phenotype.txt").read().split("\n"):
    if len(line) < 2: continue
    tokens = line.split("\t")
    if len(tokens) < 3 or tokens[1] == "-": continue
    ncbi_id, ncbi_name, hpo_id, orphanet_id = tokens[0], tokens[1], tokens[2], tokens[5].replace("ORPHA:", "Orphanet:")
    ncbi_id = "NCBI:" + ncbi_id
    if ncbi_id not in gene_terms: gene_terms[ncbi_id] = Term(ncbi_id, ncbi_name + " (Gene)", {hpo_id})
    else: gene_terms[ncbi_id].relations.add(hpo_id)
    if orphanet_id in orphanet_to_mondo: gene_terms[ncbi_id].relations.add(orphanet_to_mondo[orphanet_id])




print("Readin genes_to_disease")
for line in open("Data/genes_to_disease.txt").read().split("\n"):
    if len(line) < 2: continue
    tokens = line.split("\t")
    if len(tokens) < 3 or tokens[1] == "-": continue
    ncbi_id, ncbi_name, omim_id = tokens[0], tokens[1], tokens[3]
    if omim_id not in omim_to_mondo: continue
    ncbi_id = ncbi_id.replace("NCBIGene:", "NCBI:")
    if ncbi_id not in gene_terms: gene_terms[ncbi_id] = Term(ncbi_id, ncbi_name + " (Gene)", {omim_to_mondo[omim_id]})
    else: gene_terms[ncbi_id].relations.add(omim_to_mondo[omim_id])



with open("Data/HGNC.csv") as infile:
    ensembl_organism = ensembl_terms_organisms["9606"]
    print("Readin HGNC")
    for i, line in enumerate(infile):
        if i == 0: continue
        tokens = line.strip().split("\t")
        if len(tokens) < 2: continue
        hgnc_id = tokens[0]
        hgnc_name = tokens[1] + " (Gene)"
        hgnc_term = Term(hgnc_id, hgnc_name)

        uniprot = tokens[2] if len(tokens) > 2 and len(tokens[2]) > 0 else None
        ncbi_id = "NCBI:" + tokens[3] if len(tokens) > 3 and len(tokens[3]) > 0 else None
        ENG = tokens[4] if len(tokens) > 4 and len(tokens[4]) > 0 else None

        if uniprot != None and uniprot in uniprot_data:
            uniprot_data[uniprot].relations.add(hgnc_id)

        if ncbi_id != None and ncbi_id in gene_terms:
            gene_terms[ncbi_id].merge(hgnc_term)
            gene_terms[hgnc_id] = gene_terms[ncbi_id]

        else:
            gene_terms[hgnc_id] = hgnc_term
            if ncbi_id != None:
                hgnc_term.id.add(ncbi_id)
                gene_terms[ncbi_id] = hgnc_term

        if ENG != None:
            if ENG in ensembl_organism: ensembl_organism[ENG].relations.add(hgnc_id)

            else:
                ensembl_organism[ENG] = Term(ENG, ENG, {"LS:0000006", hgnc_id})


print("Readin hgnc_to_mondo")
for line in open("Data/hgnc_to_mondo.csv").read().split("\n"):
    if len(line) < 1: continue
    tokens = line.split("\t")
    if len(tokens) < 2: continue
    hgnc_id, mondo_ids = tokens[0], set(tokens[1:])
    if hgnc_id in gene_terms: gene_terms[hgnc_id].relations |= mondo_ids





go_terms = get_terms("Data/go-basic.obo", "GO:", True)
namespaces = {
    "biological_process": "Biological process",
    "cellular_component": "Cellular component",
    "molecular_function": "Molecular function"
}
for go_term_id, go_term in go_terms.items():
    go_term.namespace = {namespaces[n] if n in namespaces else n for n in go_term.namespace}

ontology_terms, chebi_lipid_terms = get_lipid_terms()



# do several organisms
for tax_name, tax_id in species.items():
    print(tax_id)

    output = []

    for go_id, term in go_terms.items():
        output.append(term.to_string())

    output.append("LS:0000004\tUniprot protein\t\t\n")
    output.append("LS:0000005\tChEBI metabolite\t\t\n")
    output.append("LS:0000006\tEnsembl transcript\t\t\n") # ENS

    # create organism-specific rhea terms
    rhea_terms = {}
    for rhea_term_id, uniprot_terms in rhea_terms_organisms[tax_id].items():
        rhea_terms[rhea_term_id] = Term(
            f"RHEA:{rhea_term_id}",
            f"Rhea reaction {rhea_term_id}",
            set(["UNIPROT:" + uniprot for uniprot in uniprot_terms])
        )


    chebi_terms = {t: c.copy() for t, c in chebi_terms_all.items()}
    for chebi_term, pathway_terms in chebi_terms_organisms[tax_id].items():
        if chebi_term not in chebi_terms: continue
        chebi_terms[chebi_term].relations |= pathway_terms

    protein_terms = uniprot_terms_organisms[tax_id]

    if tax_id in pathway_terms_organisms:
        pathway_terms = pathway_terms_organisms[tax_id]
        for pathway_id, pathway_term in pathway_terms.items():
            output.append(pathway_term.to_string())

    for _, protein_term in protein_terms.items():
        output.append(protein_term.to_string())


    if tax_id in ensembl_terms_organisms:
        ensembl_terms = ensembl_terms_organisms[tax_id]
        for _, ensembl_term in ensembl_terms.items():
            output.append(ensembl_term.to_string())


    written_genes = set()
    for _, gene_term in gene_terms.items():
        if gene_term in written_genes: continue
        written_genes.add(gene_term)
        output.append(gene_term.to_string())


    for _, rhea_term in rhea_terms.items():
        output.append(rhea_term.to_string())


    for _, chebi_term in chebi_terms.items():
        output.append(chebi_term.to_string())


    for lipid_name, term in chebi_lipid_terms.items():
        output.append(term.to_string())



    # diseases and phenotypes
    written_diseases = set()
    for _, disease_term in disease_and_phenotype_terms.items():
        if disease_term in written_diseases: continue
        written_diseases.add(disease_term)
        output.append(disease_term.to_string())



    for term_id, term in ontology_terms.items():
        output.append(term.to_string())


    with gzip.open(f"../ontology_{tax_id}.gz", "wb") as gz_output:
        gzip_out = "".join(output).encode("utf8")
        print("writing", len(gzip_out))
        gz_output.write(gzip_out)















