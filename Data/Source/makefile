DOWNLOAD_FILES = Data/doid-base.obo Data/go-basic.obo Data/chebi_lite.obo Data/chebi.obo.gz Data/rhea-reactions.txt.gz Data/rhea2uniprot_sprot.tsv Data/rhea2uniprot_trembl.tsv.gz Data/pathbank_all_metabolites.csv.zip Data/pathbank_all_proteins.csv.zip Data/mondo.obo Data/causal_gene_to_disease_association.all.tsv.gz Data/correlated_gene_to_disease_association.all.tsv.gz Data/hp.obo Data/genes_to_phenotype.txt Data/genes_to_disease.txt Data/Homo_sapiens.uniprot.tsv.gz Data/Homo_sapiens.all.tsv.gz Data/Mus_musculus.uniprot.tsv.gz Data/Mus_musculus.all.tsv.gz Data/HGNC.csv Data/uniprot.csv.gz Data/UniProt2Reactome_All_Levels.txt.gz Data/ChEBI2Reactome_All_Levels.txt.gz Data/Ensembl2Reactome_All_Levels.txt.gz Data/ReactomePathwaysRelation.txt

SCRIPT_FILES = Data/uniprot_to_doid.csv Data/chebi.csv Data/rhea_to_chebi.csv Data/pathbank_to_uniprot.csv Data/chebi_to_pathway.csv Data/hgnc_to_mondo.csv

ALL_SPECIES = Data/Saccharomyces_cerevisiae.uniprot.tsv.gz Data/Drosophila_melanogaster.uniprot.tsv.gz Data/Rattus_norvegicus.uniprot.tsv.gz Data/Bos_taurus.uniprot.tsv.gz Data/Caenorhabditis_elegans.uniprot.tsv.gz Data/Saccharomyces_cerevisiae.all.tsv.gz Data/Drosophila_melanogaster.all.tsv.gz Data/Rattus_norvegicus.all.tsv.gz Data/Bos_taurus.all.tsv.gz Data/Caenorhabditis_elegans.all.tsv.gz

######## OWN DATA ########
# Data/LION.obo
# Data/chebi_to_lipids.csv
# Data/ChEBI_chebi_to_lipids.csv
# Data/LipidMaps_chebi_to_lipids.csv
# Data/SwissLipids_chebi_to_lipids.csv


######## CREATE ONTOLOGIES ########
main: $(DOWNLOAD_FILES) $(SCRIPT_FILES)
	python3 create_ontology.py

all: $(DOWNLOAD_FILES) $(ALL_SPECIES) $(SCRIPT_FILES)
	python3 create_ontology.py all



####### DOWNLOAD DATA ########
#Disease Ontology
Data/doid-base.obo:
	wget --tries=10 -O Data/doid-base.obo "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/refs/heads/main/src/ontology/doid-base.obo"

Data/go-basic.obo:
	# Gene Ontology
	wget --tries=10 -O Data/go-basic.obo https://current.geneontology.org/ontology/go-basic.obo

Data/chebi_lite.obo:
	# ChEBI
	wget --tries=10 -O Data/chebi_lite.obo ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi_lite.obo

Data/chebi.obo.gz:
	# ChEBI
	wget --tries=10 -O Data/chebi.obo ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo
	gzip Data/chebi.obo

# RHEA
Data/rhea-reactions.txt.gz:
	wget --tries=10 -O Data/rhea-reactions.txt.gz ftp://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz

Data/rhea2uniprot_sprot.tsv:
	wget --tries=10 -O Data/rhea2uniprot_sprot.tsv ftp://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot%5Fsprot.tsv

Data/rhea2uniprot_trembl.tsv.gz:
	wget --tries=10 -O Data/rhea2uniprot_trembl.tsv.gz ftp://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_trembl.tsv.gz

# PathBank
Data/pathbank_all_metabolites.csv.zip:
	wget --tries=10 -O Data/pathbank_all_metabolites.csv.zip https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip

Data/pathbank_all_proteins.csv.zip:
	wget --tries=10 -O Data/pathbank_all_proteins.csv.zip https://pathbank.org/downloads/pathbank_all_proteins.csv.zip

# Reactome
Data/UniProt2Reactome_All_Levels.txt.gz:
	wget --tries=10 -O Data/UniProt2Reactome_All_Levels.txt https://reactome.org/download/current/UniProt2Reactome_All_Levels.txt
	gzip Data/UniProt2Reactome_All_Levels.txt

Data/ChEBI2Reactome_All_Levels.txt.gz:
	wget --tries=10 -O Data/ChEBI2Reactome_All_Levels.txt https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt
	gzip Data/ChEBI2Reactome_All_Levels.txt

Data/Ensembl2Reactome_All_Levels.txt.gz:
	wget --tries=10 -O Data/Ensembl2Reactome_All_Levels.txt https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt
	gzip Data/Ensembl2Reactome_All_Levels.txt


Data/ReactomePathwaysRelation.txt:
	wget --tries=10 -O Data/ReactomePathwaysRelation.txt https://reactome.org/download/current/ReactomePathwaysRelation.txt

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
	wget --tries=10 -O Data/Homo_sapiens.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.114.uniprot.tsv.gz"

Data/Homo_sapiens.all.tsv.gz:
	wget --tries=10 -O Data/Homo_sapiens.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"

Data/Mus_musculus.uniprot.tsv.gz:
	wget --tries=10 -O Data/Mus_musculus.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/mus_musculus/Mus_musculus.GRCm39.114.uniprot.tsv.gz"

Data/Mus_musculus.all.tsv.gz:
	wget --tries=10 -O Data/Mus_musculus.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz"

Data/Saccharomyces_cerevisiae.uniprot.tsv.gz:
	wget --tries=10 -O Data/Saccharomyces_cerevisiae.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.114.uniprot.tsv.gz"

Data/Saccharomyces_cerevisiae.all.tsv.gz:
	wget --tries=10 -O Data/Saccharomyces_cerevisiae.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz"

Data/Drosophila_melanogaster.uniprot.tsv.gz:
	wget --tries=10 -O Data/Drosophila_melanogaster.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.114.uniprot.tsv.gz"

Data/Drosophila_melanogaster.all.tsv.gz:
	wget --tries=10 -O Data/Drosophila_melanogaster.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.54.cds.all.fa.gz"

Data/Rattus_norvegicus.uniprot.tsv.gz:
	wget --tries=10 -O Data/Rattus_norvegicus.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/rattus_norvegicus/Rattus_norvegicus.GRCr8.114.uniprot.tsv.gz"

Data/Rattus_norvegicus.all.tsv.gz:
	wget --tries=10 -O Data/Rattus_norvegicus.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/rattus_norvegicus/cds/Rattus_norvegicus.GRCr8.cds.all.fa.gz"

Data/Bos_taurus.uniprot.tsv.gz:
	wget --tries=10 -O Data/Bos_taurus.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/bos_taurus/Bos_taurus.ARS-UCD2.0.114.uniprot.tsv.gz"

Data/Bos_taurus.all.tsv.gz:
	wget --tries=10 -O Data/Bos_taurus.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/bos_taurus/cds/Bos_taurus.ARS-UCD2.0.cds.all.fa.gz"

Data/Caenorhabditis_elegans.uniprot.tsv.gz:
	wget --tries=10 -O Data/Caenorhabditis_elegans.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.cds.all.fa.gz"

Data/Caenorhabditis_elegans.all.tsv.gz:
	wget --tries=10 -O Data/Caenorhabditis_elegans.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/caenorhabditis_elegans/cds/Caenorhabditis_elegans.WBcel235.cds.all.fa.gz"

# HGNC
Data/HGNC.csv:
	wget --tries=10 -O Data/HGNC.csv "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=md_prot_id&col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"

# Uniprot
Data/uniprot.csv.gz:
	# revieved = True
	wget --tries=10 -O Data/uniprot.csv.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Corganism_id%2Cgene_primary%2Cgo%2Cec%2Cprotein_families%2Cxref_mim_full%2Cxref_hgnc_full%2Cxref_ensembl_full%2Cxref_geneid%2Creviewed&format=tsv&query=%28%28taxonomy_id%3A9606%29+OR+%28taxonomy_id%3A10090%29+OR+%28taxonomy_id%3A4932%29+OR+%28taxonomy_id%3A562%29+OR+%28taxonomy_id%3A7227%29+OR+%28taxonomy_id%3A10116%29+OR+%28taxonomy_id%3A9913%29+OR+%28taxonomy_id%3A6239%29+OR+%28taxonomy_id%3A237561%29+OR+%28taxonomy_id%3A405534%29+OR+%28taxonomy_id%3A3702%29%29"



######## EXECUTE SCRIPTS ########
# Disease Ontology
Data/uniprot_to_doid.csv: Data/uniprot.csv.gz Data/doid-base.obo
	python3 Scripts/uniprot_to_doid.py

# ChEBI
Data/chebi.csv: Data/chebi_lite.obo
	python3 Scripts/chebi.py

# RHEA
Data/rhea_to_chebi.csv: Data/rhea-reactions.txt.gz
	python3 Scripts/rhea_to_chebi.py

# Pathbank
Data/pathbank_to_uniprot.csv: Data/pathbank_all_proteins.csv.zip
	python3 Scripts/pathbank_to_uniprot.py

Data/chebi_to_pathway.csv: Data/pathbank_all_metabolites.csv.zip
	python3 Scripts/chebi_to_pathway.py

# MONDO
Data/hgnc_to_mondo.csv: Data/causal_gene_to_disease_association.all.tsv.gz Data/correlated_gene_to_disease_association.all.tsv.gz
	python3 Scripts/hgnc_to_mondo.py




######## CLEAN WORKSPACE ########
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
	rm -f Data/pathbank_all_metabolites.csv.zip
	rm -f Data/pathbank_all_proteins.csv.zip
	rm -f Data/pathbank_to_uniprot.csv
	rm -f Data/rhea2uniprot_sprot.tsv
	rm -f Data/rhea-reactions.txt.gz
	rm -f Data/rhea_to_chebi.csv
	rm -f Data/uniprot.csv.gz
	rm -f Data/uniprot_to_doid.csv
