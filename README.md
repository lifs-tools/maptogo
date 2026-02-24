# MAPtoGO
Gene Ontology (GO) enrichment analysis is commonly used to interpret genomics, transcriptomics, or proteomics data by linking expression products to biological processes. However, comparable tools for lipidomics and metabolomics are limited, as small molecules typically lack direct GO annotations. **MAPtoGO** a GO-based multi-omics enrichment approach that connects lipids and metabolites to GO terms by mapping them to related proteins via curated enzymatic networks. This unified framework integrates lipidomics, metabolomics, proteomics, and transcriptomics, expanding GO coverage to small molecules and enhancing annotation accuracy.

## Running instance
A **MAPtoGO** instance with 11 predefined model organism knowledge graphs can be accessed from our portal [[https://lifs-tools.org/maptogo|https://lifs-tools.org/maptogo]].

## Building and running locally
To build and run the **MAPtoGO** Docker image locally, follow these steps:

1. Clone the repository:
   ```bash
   git clone
   cd maptogo
   ```
2. Build the Docker image:
   ```bash
   docker build -t maptogo:latest .
   ```
3. Run the Docker container:
   ```bash
   docker run -p 8040:8040 maptogo:latest
   ```
4. Access the application in your web browser at `http://localhost:8040`.
