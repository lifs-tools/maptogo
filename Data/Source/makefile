DOWNLOAD_FILES = Data/doid-base.obo Data/go-basic.obo Data/chebi_lite.obo Data/rhea-reactions.txt.gz Data/rhea2uniprot_sprot.tsv Data/pathbank_all_metabolites.csv Data/pathbank_all_proteins.csv Data/mondo.obo Data/causal_gene_to_disease_association.all.tsv.gz Data/correlated_gene_to_disease_association.all.tsv.gz Data/hp.obo Data/genes_to_phenotype.txt Data/genes_to_disease.txt Data/Homo_sapiens.uniprot.tsv.gz Data/HGNC.csv Data/uniprot.csv.gz

SCRIPT_FILES = Data/uniprot_to_doid.csv Data/uniprot_to_go.csv Data/chebi.csv Data/rhea_to_chebi.csv Data/pathbank_to_uniprot.csv Data/chebi_to_pathway.csv Data/hgnc_to_mondo.csv


####### OWN DATA ########
#Data/LION.obo
#Data/chebi_to_lipids.csv
#Data/ChEBI_chebi_to_lipids.csv
#Data/LipidMaps_chebi_to_lipids.csv
#Data/SwissLipids_chebi_to_lipids.csv


######## CREATE ONTOLOGIES ########
main: $(DOWNLOAD_FILES) $(SCRIPT_FILES)
	python3 create_ontology.py



####### DOWNLOAD DATA ########
#Disease Ontology
Data/doid-base.obo:
	wget --tries=10 -O Data/doid-base.obo "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/refs/heads/main/src/ontology/doid-base.obo"

Data/go-basic.obo:
	# Gene Ontology
	wget --tries=10 -O Data/go-basic.obo https://current.geneontology.org/ontology/go-basic.obo

Data/chebi_lite.obo:
	# ChEBI
	wget --tries=10 -O Data/chebi_lite.obo https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo

# RHEA
Data/rhea-reactions.txt.gz:
	wget --tries=10 -O Data/rhea-reactions.txt.gz https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz

Data/rhea2uniprot_sprot.tsv:
	wget --tries=10 -O Data/rhea2uniprot_sprot.tsv https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot%5Fsprot.tsv

# PathBank
Data/pathbank_all_metabolites.csv:
	wget --tries=10 -O Data/pathbank_all_metabolites.csv.zip https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip
	unzip -d Data Data/pathbank_all_metabolites.csv.zip

Data/pathbank_all_proteins.csv:
	wget --tries=10 -O Data/pathbank_all_proteins.csv.zip https://pathbank.org/downloads/pathbank_all_proteins.csv.zip
	unzip -d Data Data/pathbank_all_proteins.csv.zip

# Mondo
Data/mondo.obo:
	wget --tries=10 -O Data/mondo.obo https://purl.obolibrary.org/obo/mondo.obo

Data/causal_gene_to_disease_association.all.tsv.gz:
	wget --tries=10 -O Data/causal_gene_to_disease_association.all.tsv.gz https://data.monarchinitiative.org/monarch-kg/latest/tsv/all_associations/causal_gene_to_disease_association.all.tsv.gz

Data/correlated_gene_to_disease_association.all.tsv.gz:
	wget --tries=10 -O Data/correlated_gene_to_disease_association.all.tsv.gz https://data.monarchinitiative.org/monarch-kg/latest/tsv/all_associations/correlated_gene_to_disease_association.all.tsv.gz

# HP
Data/hp.obo:
	wget --tries=10 -O Data/hp.obo https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo

Data/genes_to_phenotype.txt:
	wget --tries=10 -O Data/genes_to_phenotype.txt https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt

Data/genes_to_disease.txt:
	wget --tries=10 -O Data/genes_to_disease.txt https://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt

# Ensembl
Data/Homo_sapiens.uniprot.tsv.gz:
	wget --tries=10 -O Data/Homo_sapiens.uniprot.tsv.gz "https://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.114.uniprot.tsv.gz"

# HGNC
Data/HGNC.csv:
	wget --tries=10 -O Data/HGNC.csv "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=md_prot_id&col=gd_pub_eg_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"

# Uniprot
Data/uniprot.csv.gz:
	wget --tries=10 -O Data/uniprot.csv.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Corganism_id%2Cgene_primary%2Cgo%2Cprotein_families%2Cxref_mim_full%2Cxref_hgnc_full%2Cxref_ensembl_full&format=tsv&query=%28reviewed%3Atrue%29"



######## EXECUTE SCRIPTS ########
# Disease Ontology
Data/uniprot_to_doid.csv: Data/uniprot.csv.gz Data/doid-base.obo
	python3 Scripts/uniprot_to_doid.py

# Uniprot
Data/uniprot_to_go.csv: Data/uniprot.csv.gz
	python3 Scripts/uniprot_to_go.py

# ChEBI
Data/chebi.csv: Data/chebi_lite.obo
	python3 Scripts/chebi.py

# RHEA
Data/rhea_to_chebi.csv: Data/rhea-reactions.txt.gz
	python3 Scripts/rhea_to_chebi.py

# Pathbank
Data/pathbank_to_uniprot.csv: Data/pathbank_all_proteins.csv
	python3 Scripts/pathbank_to_uniprot.py

Data/chebi_to_pathway.csv: Data/pathbank_all_metabolites.csv
	python3 Scripts/chebi_to_pathway.py

# MONDO
Data/hgnc_to_mondo.csv: Data/causal_gene_to_disease_association.all.tsv.gz Data/correlated_gene_to_disease_association.all.tsv.gz
	python3 Scripts/hgnc_to_mondo.py





clean:
	rm -f Data/causal_gene_to_disease_association.all.tsv.gz
	rm -f Data/chebi.csv
	rm -f Data/chebi_lite.obo
	rm -f Data/chebi_to_pathway.csv
	rm -f Data/correlated_gene_to_disease_association.all.tsv.gz
	rm -f Data/doid-base.obo
	rm -f Data/genes_to_disease.txt
	rm -f Data/genes_to_phenotype.txt
	rm -f Data/go-basic.obo
	rm -f Data/HGNC.csv
	rm -f Data/hgnc_to_mondo.csv
	rm -f Data/Homo_sapiens.uniprot.tsv.gz
	rm -f Data/hp.obo
	rm -f Data/mondo.obo
	rm -f Data/pathbank_all_metabolites.csv
	rm -f Data/pathbank_all_metabolites.csv.zip
	rm -f Data/pathbank_all_proteins.csv
	rm -f Data/pathbank_all_proteins.csv.zip
	rm -f Data/pathbank_to_uniprot.csv
	rm -f Data/rhea2uniprot_sprot.tsv
	rm -f Data/rhea-reactions.txt.gz
	rm -f Data/rhea_to_chebi.csv
	rm -f Data/uniprot.csv.gz
	rm -f Data/uniprot_to_doid.csv
	rm -f Data/uniprot_to_go.csv
