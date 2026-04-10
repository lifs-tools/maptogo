# MAPtoGO
Gene Ontology (GO) enrichment analysis is commonly used to interpret genomics, transcriptomics, or proteomics data by linking expression products to biological processes. However, comparable tools for lipidomics and metabolomics are limited, as small molecules typically lack direct GO annotations. **MAPtoGO** is a GO-based multi-omics enrichment approach that connects lipids and metabolites to GO terms by mapping them to related proteins via curated enzymatic networks. This unified framework integrates lipidomics, metabolomics, proteomics, and transcriptomics, expanding GO coverage to small molecules and enhancing annotation accuracy.

## Running instance
A **MAPtoGO** instance with 11 predefined model organism knowledge graphs can be accessed from our portal [https://lifs-tools.org/maptogo](https://lifs-tools.org/maptogo).

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


## Use MAPtoGO with a REST API
The **MAPtoGO** instance at [https://maptogo.lifs-tools.org/](https://maptogo.lifs-tools.org/) has a built-in API. You can find the documentation [here](https://maptogo.lifs-tools.org/api/docs|here). We also provide an example shell script to send an API request:
```bash
cd examples
./api_request.sh
```
Please dissect the *api_request.sh* and adapt it for your necessities.


## Use MAPtoGO as an offline Python library
**MAPtoGO** can also be integrated directly into a Python script without the need for the API. We provide an example Python script:
```bash
cd examples
python3 use_MAPtoGO_library.py
```
Please dissect the *use_MAPtoGO_library.py* and adapt it for your necessities.


## Upload data to MAPtoGO from third-party tools
You can connect your biomolecule analysis tool with **MAPtoGO** by simply uploading the statistical results computed in your analysis tool. We provided two standalone code snippets that you can incorporate into your tool, one in Python:
```bash
cd examples
python3 open_maptogo_and_send_data.py
```

and one in C++:
```bash
cd examples
sudo apt install libcurl4-openssl-dev
g++ -std=c++17 -O3 -o open_maptogo_and_send_data open_maptogo_and_send_data.cpp -lcurl
./open_maptogo_and_send_data
```
