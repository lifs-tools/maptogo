######## OWN DATA ########
# Data/LION.obo
# Data/chebi_to_lipids.csv
# Data/ChEBI_chebi_to_lipids.csv
# Data/LipidMaps_chebi_to_lipids.csv
# Data/SwissLipids_chebi_to_lipids.csv


######## DOWNLOAD DATA ########
# Disease Ontology
# wget --tries=10 -O Data/doid-base.obo "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/refs/heads/main/src/ontology/doid-base.obo"
#
# # Uniprot
# wget --tries=10 -O Data/uniprot.csv.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Corganism_id%2Cgene_primary%2Cgo%2Cprotein_families%2Cxref_mim_full%2Cxref_hgnc_full%2Cxref_ensembl_full&format=tsv&query=%28reviewed%3Atrue%29"
#
# # Gene Ontology
# wget --tries=10 -O Data/go-basic.obo https://current.geneontology.org/ontology/go-basic.obo
#
# # ChEBI
# wget --tries=10 -O Data/chebi_lite.obo https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo
#
# # RHEA
# wget --tries=10 -O Data/rhea-reactions.txt.gz https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz
# wget --tries=10 -O Data/rhea2uniprot_sprot.tsv https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot%5Fsprot.tsv
#
# # PathBank
# wget --tries=10 -O Data/pathbank_all_metabolites.csv.zip https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip
# unzip -d Data Data/pathbank_all_metabolites.csv.zip
# wget --tries=10 -O Data/pathbank_all_proteins.csv.zip https://pathbank.org/downloads/pathbank_all_proteins.csv.zip
# unzip -d Data Data/pathbank_all_proteins.csv.zip
#
# # Mondo
# wget --tries=10 -O Data/mondo.obo https://purl.obolibrary.org/obo/mondo.obo
# wget --tries=10 -O Data/causal_gene_to_disease_association.all.tsv.gz https://data.monarchinitiative.org/monarch-kg/latest/tsv/all_associations/causal_gene_to_disease_association.all.tsv.gz
# wget --tries=10 -O Data/correlated_gene_to_disease_association.all.tsv.gz https://data.monarchinitiative.org/monarch-kg/latest/tsv/all_associations/correlated_gene_to_disease_association.all.tsv.gz
#
# # HP
# wget --tries=10 -O Data/hp.obo https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo
# wget --tries=10 -O Data/genes_to_phenotype.txt https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt
# wget --tries=10 -O Data/genes_to_disease.txt https://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt
#
# # Ensembl
# wget --tries=10 -O Data/Homo_sapiens.uniprot.tsv.gz "https://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.114.uniprot.tsv.gz"
#
# # HGNC
# wget --tries=10 -O Data/HGNC.csv "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=md_prot_id&col=gd_pub_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"





######## EXECUTE SCRIPTS ########
# Disease Ontology
# Scripts/uniprot_to_doid.py -> Data/uniprot_to_doid.csv
python3 Scripts/uniprot_to_doid.py

# Uniprot
# Scripts/uniprot_to_go.py -> Data/uniprot_to_go.csv
python3 Scripts/uniprot_to_go.py

# ChEBI
# Scripts/chebi.py -> Data/chebi.csv
python3 Scripts/chebi.py

# RHEA
# Scripts/rhea_to_chebi.py -> Data/rhea_to_chebi.csv
python3 Scripts/rhea_to_chebi.py

# Pathbank
# Scripts/pathbank_to_uniprot.py -> Data/pathbank_to_uniprot.csv
# Scripts/chebi_to_pathway.py -> Data/chebi_to_pathway.csv
python3 Scripts/pathbank_to_uniprot.py
python3 Scripts/chebi_to_pathway.py

# MONDO
# Scripts/hgnc_to_mondo.py -> Data/hgnc_to_mondo.csv
python3 Scripts/hgnc_to_mondo.py





######## CREATE ONTOLOGIES ########
python3 create_ontology.py
