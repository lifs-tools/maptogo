DOWNLOAD_FILES = Data/doid-base.obo Data/go-basic.obo Data/chebi_lite.obo Data/chebi.obo.gz Data/rhea-reactions.txt.gz Data/rhea2uniprot_sprot.tsv Data/rhea2uniprot_trembl.tsv.gz Data/mondo.obo Data/causal_gene_to_disease_association.all.tsv.gz Data/correlated_gene_to_disease_association.all.tsv.gz Data/hp.obo Data/genes_to_phenotype.txt Data/genes_to_disease.txt Data/Homo_sapiens.uniprot.tsv.gz Data/Homo_sapiens.all.tsv.gz Data/Mus_musculus.uniprot.tsv.gz Data/Mus_musculus.all.tsv.gz Data/HGNC.csv Data/uniprot.csv.gz Data/UniProt2Reactome_All_Levels.txt.gz Data/ChEBI2Reactome_All_Levels.txt.gz Data/Ensembl2Reactome_All_Levels.txt.gz Data/ReactomePathwaysRelation.txt Data/ChEBI2ReactomeReactions.txt Data/UniProt2ReactomeReactions.txt # Data/pathbank_all_metabolites.csv.zip Data/pathbank_all_proteins.csv.zip

SCRIPT_FILES = Data/uniprot_to_doid.csv Data/chebi.csv Data/rhea_to_chebi.csv Data/hgnc_to_mondo.csv ../CHEBI_synonyms.csv.gz Data/goa_uniprot.csv.gz # Data/pathbank_to_uniprot.csv Data/chebi_to_pathway.csv

ALL_SPECIES = Data/Saccharomyces_cerevisiae.uniprot.tsv.gz Data/Saccharomyces_cerevisiae.all.tsv.gz Data/Drosophila_melanogaster.uniprot.tsv.gz Data/Drosophila_melanogaster.all.tsv.gz Data/Rattus_norvegicus.uniprot.tsv.gz Data/Rattus_norvegicus.all.tsv.gz Data/Bos_taurus.uniprot.tsv.gz Data/Bos_taurus.all.tsv.gz Data/Caenorhabditis_elegans.uniprot.tsv.gz Data/Caenorhabditis_elegans.all.tsv.gz Data/Danio_rerio.uniprot.tsv.gz Data/Danio_rerio.all.tsv.gz Data/Gallus_gallus.uniprot.tsv.gz Data/Gallus_gallus.all.tsv.gz Data/Sus_scrofa.uniprot.tsv.gz Data/Sus_scrofa.all.tsv.gz

SEMANTICS_FILES = Data/phenotype_descriptions.csv Data/reactome_descriptions.csv Data/enzyme_descriptions.csv Data/go_descriptions.csv

ENSEMBL_VERSION_ID = 115

######## OWN DATA ########
# Data/LION.obo
# Data/LM.csv
# Data/chebi_to_lipids.csv
# Data/ChEBI_chebi_to_lipids.csv
# Data/LipidMaps_chebi_to_lipids.csv
# Data/SwissLipids_chebi_to_lipids.csv


######## CREATE ONTOLOGIES ########
main: $(DOWNLOAD_FILES) $(SCRIPT_FILES)
	python3 create_ontology.py

all: $(DOWNLOAD_FILES) $(ALL_SPECIES) $(SCRIPT_FILES)
	python3 create_ontology.py all


semantics: $(SEMANTICS_FILES)
	echo "huhu"



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
# Data/pathbank_all_metabolites.csv.zip:
# 	wget --tries=10 -O Data/pathbank_all_metabolites.csv.zip https://pathbank.org/downloads/pathbank_all_metabolites.csv.zip
#
# Data/pathbank_all_proteins.csv.zip:
# 	wget --tries=10 -O Data/pathbank_all_proteins.csv.zip https://pathbank.org/downloads/pathbank_all_proteins.csv.zip



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

Data/ChEBI2ReactomeReactions.txt:
	wget --tries=10 -O Data/ChEBI2ReactomeReactions.txt https://download.reactome.org/96/ChEBI2ReactomeReactions.txt

Data/UniProt2ReactomeReactions.txt:
	wget --tries=10 -O Data/UniProt2ReactomeReactions.txt https://download.reactome.org/96/UniProt2ReactomeReactions.txt




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
	wget --tries=10 -O Data/genes_to_phenotype.txt https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-02-16/genes_to_phenotype.txt

Data/genes_to_disease.txt:
	wget --tries=10 -O Data/genes_to_disease.txt https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-02-16/genes_to_disease.txt



# Ensembl
Data/Homo_sapiens.uniprot.tsv.gz:
	wget --tries=10 -O Data/Homo_sapiens.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Homo_sapiens.all.tsv.gz:
	wget --tries=10 -O Data/Homo_sapiens.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"

Data/Mus_musculus.uniprot.tsv.gz:
	wget --tries=10 -O Data/Mus_musculus.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/mus_musculus/Mus_musculus.GRCm39.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Mus_musculus.all.tsv.gz:
	wget --tries=10 -O Data/Mus_musculus.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz"

Data/Saccharomyces_cerevisiae.uniprot.tsv.gz:
	wget --tries=10 -O Data/Saccharomyces_cerevisiae.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Saccharomyces_cerevisiae.all.tsv.gz:
	wget --tries=10 -O Data/Saccharomyces_cerevisiae.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz"

Data/Drosophila_melanogaster.uniprot.tsv.gz:
	wget --tries=10 -O Data/Drosophila_melanogaster.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Drosophila_melanogaster.all.tsv.gz:
	wget --tries=10 -O Data/Drosophila_melanogaster.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/drosophila_melanogaster/cds/Drosophila_melanogaster.BDGP6.54.cds.all.fa.gz"

Data/Rattus_norvegicus.uniprot.tsv.gz:
	wget --tries=10 -O Data/Rattus_norvegicus.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/rattus_norvegicus/Rattus_norvegicus.GRCr8.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Rattus_norvegicus.all.tsv.gz:
	wget --tries=10 -O Data/Rattus_norvegicus.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/rattus_norvegicus/cds/Rattus_norvegicus.GRCr8.cds.all.fa.gz"

Data/Bos_taurus.uniprot.tsv.gz:
	wget --tries=10 -O Data/Bos_taurus.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/bos_taurus/Bos_taurus.ARS-UCD2.0.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Bos_taurus.all.tsv.gz:
	wget --tries=10 -O Data/Bos_taurus.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/bos_taurus/cds/Bos_taurus.ARS-UCD2.0.cds.all.fa.gz"

Data/Caenorhabditis_elegans.uniprot.tsv.gz:
	wget --tries=10 -O Data/Caenorhabditis_elegans.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Caenorhabditis_elegans.all.tsv.gz:
	wget --tries=10 -O Data/Caenorhabditis_elegans.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/caenorhabditis_elegans/cds/Caenorhabditis_elegans.WBcel235.cds.all.fa.gz"

Data/Danio_rerio.uniprot.tsv.gz:
	wget --tries=10 -O Data/Danio_rerio.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/danio_rerio/Danio_rerio.GRCz11.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Danio_rerio.all.tsv.gz:
	wget --tries=10 -O Data/Danio_rerio.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/danio_rerio/cds/Danio_rerio.GRCz11.cds.all.fa.gz"

Data/Gallus_gallus.uniprot.tsv.gz:
	wget --tries=10 -O Data/Gallus_gallus.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Gallus_gallus.all.tsv.gz:
	wget --tries=10 -O Data/Gallus_gallus.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/gallus_gallus/cds/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cds.all.fa.gz"

Data/Sus_scrofa.uniprot.tsv.gz:
	wget --tries=10 -O Data/Sus_scrofa.uniprot.tsv.gz "ftp://ftp.ensembl.org/pub/current_tsv/sus_scrofa/Sus_scrofa.Sscrofa11.1.${ENSEMBL_VERSION_ID}.uniprot.tsv.gz"

Data/Sus_scrofa.all.tsv.gz:
	wget --tries=10 -O Data/Sus_scrofa.all.tsv.gz "ftp://ftp.ensembl.org/pub/current/fasta/sus_scrofa/cds/Sus_scrofa.Sscrofa11.1.cds.all.fa.gz"



# HGNC
Data/HGNC.csv:
	wget --tries=10 -O Data/HGNC.csv "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=md_prot_id&col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"



# Uniprot
Data/uniprot.csv.gz:
	# revieved = True
	wget --tries=10 -O Data/uniprot.csv.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Corganism_id%2Cgene_primary%2Cgo%2Cec%2Cprotein_families%2Cxref_mim_full%2Cxref_hgnc_full%2Cxref_ensembl_full%2Cxref_geneid%2Creviewed%2Ckeyword%2Ccc_function&format=tsv&query=%28%28taxonomy_id%3A9606%29+OR+%28taxonomy_id%3A10090%29+OR+%28taxonomy_id%3A4932%29+OR+%28taxonomy_id%3A562%29+OR+%28taxonomy_id%3A7227%29+OR+%28taxonomy_id%3A10116%29+OR+%28taxonomy_id%3A9913%29+OR+%28taxonomy_id%3A6239%29+OR+%28taxonomy_id%3A237561%29+OR+%28taxonomy_id%3A405534%29+OR+%28taxonomy_id%3A3702%29+OR+%28taxonomy_id%3A4550%29+OR+%28taxonomy_id%3A7955%29+OR+%28taxonomy_id%3A1423%29+OR+%28taxonomy_id%3A4896%29+OR+%28taxonomy_id%3A9031%29+OR+%28taxonomy_id%3A9823%29%29"



# GO
Data/goa_uniprot_gcrp.gpa.gz:
	wget --tries=10 -O Data/goa_uniprot_gcrp.gpa.gz "https://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_gcrp.gpa.gz"



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
# Data/pathbank_to_uniprot.csv: Data/pathbank_all_proteins.csv.zip
# 	python3 Scripts/pathbank_to_uniprot.py
#
# Data/chebi_to_pathway.csv: Data/pathbank_all_metabolites.csv.zip
# 	python3 Scripts/chebi_to_pathway.py

# MONDO
Data/hgnc_to_mondo.csv: Data/causal_gene_to_disease_association.all.tsv.gz Data/correlated_gene_to_disease_association.all.tsv.gz
	python3 Scripts/hgnc_to_mondo.py

# CHEBI
../CHEBI_synonyms.csv.gz: Data/chebi.obo.gz
	python3 Scripts/CHEBI_synonyms.py

# GO
Data/goa_uniprot.csv.gz: Data/goa_uniprot_gcrp.gpa.gz
	python3 Scripts/goa_uniprot.py


# Descriptions
Data/reactome_descriptions.csv:
	wget --tries=10 -O Data/gk_current.sql.gz https://download.reactome.org/96/databases/gk_current.sql.gz
	rm -f Data/gk_current.sql
	gunzip Data/gk_current.sql.gz
	wget --tries=10 -O Scripts/mysql2sqlite https://raw.githubusercontent.com/dumblob/mysql2sqlite/master/mysql2sqlite
	chmod +x Scripts/mysql2sqlite
	Scripts/mysql2sqlite Data/gk_current.sql | sqlite3 Data/gk_current.sqlite
	python3 Scripts/reactome_descriptions.py
	rm -f Data/gk_current.sql
	rm -f Data/gk_current.sqlite


Data/go_descriptions.csv:
	python3 Scripts/go_descriptions.py


Data/enzyme_descriptions.csv:
	python3 Scripts/enzyme_descriptions.py


Data/phenotype_descriptions.csv:
	wget --tries=10 -O Data/hp.obo https://github.com/obophenotype/human-phenotype-ontology/raw/refs/heads/master/hp.obo
	python3 Scripts/phenotype_descriptions.py


######## CLEAN WORKSPACE ########
clean:
	rm -f Data/causal_gene_to_disease_association.all.tsv.gz
	rm -f Data/chebi.csv
	rm -f Data/chebi_lite.obo
	rm -f Data/chebi.obo.gz
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
	rm -f Data/goa_uniprot.csv.gz
	rm -f Data/goa_uniprot_gcrp.gpa.gz
	rm -f ../CHEBI_synonyms.csv.gz
# 	rm -f Data/pathbank_all_metabolites.csv.zip
# 	rm -f Data/pathbank_all_proteins.csv.zip
# 	rm -f Data/pathbank_to_uniprot.csv
	#rm -f Data/chebi_to_pathway.csv
	rm -f Data/rhea2uniprot_sprot.tsv
	rm -f Data/rhea2uniprot_trembl.tsv.gz
	rm -f Data/rhea-reactions.txt.gz
	rm -f Data/rhea_to_chebi.csv
	rm -f Data/uniprot.csv.gz
	rm -f Data/uniprot_to_doid.csv
	rm -f Data/ReactomePathwaysRelation.txt
	rm -f Data/ChEBI2ReactomeReactions.txt
	rm -f Data/UniProt2ReactomeReactions.txt
	rm -f Data/UniProt2Reactome_All_Levels.txt.gz
	rm -f Data/ChEBI2Reactome_All_Levels.txt.gz
	rm -f Data/Ensembl2Reactome_All_Levels.txt.gz
	rm -f Data/Homo_sapiens.all.tsv.gz
	rm -f Data/Mus_musculus.uniprot.tsv.gz
	rm -f Data/Mus_musculus.all.tsv.gz
	rm -f Data/Saccharomyces_cerevisiae.uniprot.tsv.gz
	rm -f Data/Saccharomyces_cerevisiae.all.tsv.gz
	rm -f Data/Drosophila_melanogaster.uniprot.tsv.gz
	rm -f Data/Drosophila_melanogaster.all.tsv.gz
	rm -f Data/Rattus_norvegicus.uniprot.tsv.gz
	rm -f Data/Rattus_norvegicus.all.tsv.gz
	rm -f Data/Bos_taurus.uniprot.tsv.gz
	rm -f Data/Bos_taurus.all.tsv.gz
	rm -f Data/Caenorhabditis_elegans.uniprot.tsv.gz
	rm -f Data/Caenorhabditis_elegans.all.tsv.gz
	rm -f Data/Danio_rerio.uniprot.tsv.gz
	rm -f Data/Danio_rerio.all.tsv.gz
	rm -f Data/Gallus_gallus.uniprot.tsv.gz
	rm -f Data/Gallus_gallus.all.tsv.gz
	rm -f Data/Sus_scrofa.uniprot.tsv.gz
	rm -f Data/Sus_scrofa.all.tsv.gz
	rm -f Data/reactome_descriptions.csv
	rm -f Data/enzyme_descriptions.csv
	rm -f Data/go_descriptions.csv

