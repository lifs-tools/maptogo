######## OWN DATA ########
# Data/LION.obo
# Data/chebi_to_lipids.csv
# Data/ChEBI_chebi_to_lipids.csv
# Data/LipidMaps_chebi_to_lipids.csv
# Data/SwissLipids_chebi_to_lipids.csv


######## DOWNLOAD DATA ########
# Disease Ontology
#wget -O Data/doid-base.obo "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/refs/heads/main/src/ontology/doid-base.obo"

# Uniprot
#wget -O Data/uniprot.csv.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Corganism_id%2Cgene_primary%2Cgo%2Cprotein_families%2Cxref_mim_full%2Cxref_hgnc_full%2Cxref_ensembl_full&format=tsv&query=accession%3AA0A0C5B5G6+OR+accession%3AA0A1B0GTW7+OR+accession%3AA0JNW5+OR+accession%3AA0JP26"

# Gene Ontology
#wget -O Data/go-basic.obo https://current.geneontology.org/ontology/go-basic.obo

# ChEBI
#wget -O Data/chebi_lite.obo https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo

# RHEA
#wget -O Data/rhea-reactions.txt.gz https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz
#wget -O Data/rhea2uniprot_sprot.tsv https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot%5Fsprot.tsv

# PathBank
#wget -O Data/pathbank_all_metabolites.csv.zip https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip
#unzip -d Data Data/pathbank_all_metabolites.csv.zip
#wget -O Data/pathbank_all_proteins.csv.zip https://pathbank.org/downloads/pathbank_all_proteins.csv.zip
#unzip -d Data Data/pathbank_all_proteins.csv.zip

# Mondo
#wget -O Data/mondo.obo https://purl.obolibrary.org/obo/mondo.obo
#wget -O Data/causal_gene_to_disease_association.all.tsv.gz https://data.monarchinitiative.org/monarch-kg/latest/tsv/all_associations/causal_gene_to_disease_association.all.tsv.gz
#wget -O Data/correlated_gene_to_disease_association.all.tsv.gz https://data.monarchinitiative.org/monarch-kg/latest/tsv/all_associations/correlated_gene_to_disease_association.all.tsv.gz

# HP
#wget -O Data/hp.obo https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo
#wget -O Data/genes_to_phenotype.txt https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt
#wget -O Data/genes_to_disease.txt https://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt

# Ensembl
#wget -O Data/Homo_sapiens.uniprot.tsv.gz "https://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.114.uniprot.tsv.gz"

# HGNC
# wget -O Data/HGNC.csv "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=md_prot_id&col=gd_pub_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"





######## EXECUTE SCRIPTS ########
# Disease Ontology
# Scripts/uniprot_to_doid.py -> Data/uniprot_to_doid.csv

# Uniprot
# Scripts/uniprot_to_go.py -> Data/uniprot_to_go.csv

# ChEBI
# Scripts/chebi.py -> Data/chebi.csv

# RHEA
# Scripts/rhea_to_chebi.py -> Data/rhea_to_chebi.csv

# Pathbank
# Scripts/pathbank_to_uniprot.py -> Data/pathbank_to_uniprot.csv
# Scripts/chebi_to_pathway.py -> Data/chebi_to_pathway.csv

# MONDO
# Scritps/hgnc_to_mondo.py -> Data/hgnc_to_mondo.csv


