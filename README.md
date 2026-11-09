# MAPtoGO
Gene Ontology (GO) enrichment analysis is commonly used to interpret genomics, transcriptomics, or proteomics data by linking expression products to biological processes. However, comparable tools for lipidomics and metabolomics are limited, as small molecules typically lack direct GO annotations. **MAPtoGO** is a GO-based multi-omics enrichment approach that connects lipids and metabolites to GO terms by mapping them to related proteins via curated enzymatic networks. This unified framework integrates lipidomics, metabolomics, proteomics, and transcriptomics, expanding GO coverage to small molecules and enhancing annotation accuracy.

## Running instance
A **MAPtoGO** instance with 11 predefined model organism knowledge graphs can be accessed from our portal [https://lifs-tools.org/maptogo](https://lifs-tools.org/maptogo).

## Building and running locally
The **MAPtoGO** application and its knowledge graphs live in two separate repositories:

- [maptogo](https://github.com/lifs-tools/maptogo): the application code and Dockerfile
- [maptogo-data](https://github.com/lifs-tools/maptogo-data): the precomputed organism knowledge graphs (`ontology_<taxon id>.zst`), the ChEBI synonym table (`CHEBI_synonyms.csv.gz`), additional term links (`additional_links.csv`) and the example datasets (`examples.xlsx`)

The Docker image doesn't include the data. When the application starts, it loads every organism listed in `Config.py` from `/wd/Data` inside the container, so you have to mount your `maptogo-data` checkout to that path.

To build and run the **MAPtoGO** Docker image locally, follow these steps:

1. Clone both repositories next to each other:
   ```bash
   git clone https://github.com/lifs-tools/maptogo.git
   git clone https://github.com/lifs-tools/maptogo-data.git
   cd maptogo
   ```
   The `maptogo-data` repository is about 130 MB.
2. Build the Docker image from the `maptogo` directory:
   ```bash
   docker build -t maptogo:latest .
   ```
   The build compiles the C++ Fisher's exact test library and copies `Config.py.template` to `Config.py`. That file lists all 17 organisms that `maptogo-data` provides.
3. Run the Docker container with `maptogo-data` mounted read-only at `/wd/Data`:
   ```bash
   docker run --rm -it \
     -p 8040:8040 \
     -v "$(realpath ../maptogo-data)":/wd/Data:ro \
     --name maptogo \
     maptogo:latest
   ```
   The `-v` source path must be absolute, which is why the command uses `realpath`. If your `maptogo-data` checkout is somewhere else, replace `../maptogo-data` with its path. On Windows PowerShell, use `-v "${PWD}\..\maptogo-data:/wd/Data:ro"` instead.
4. Wait until the log shows that all organisms have loaded (one `Loading <organism> (<size> bytes)` line per organism) and gunicorn is listening on port 8040. Loading all knowledge graphs takes a while and needs several GB of memory, so give Docker Desktop enough RAM.
5. Open the application in your web browser at [http://localhost:8040](http://localhost:8040).

To stop the container, press `Ctrl+C` or run `docker stop maptogo` from another terminal.

### Troubleshooting
- **`FileNotFoundError: .../Data/ontology_<taxon id>.zst` at startup:** the data directory isn't mounted at `/wd/Data`, or the mount source path is wrong. Run `docker run --rm -v "$(realpath ../maptogo-data)":/wd/Data:ro maptogo:latest ls /wd/Data` and check that the `ontology_*.zst` files are listed.
- **The container exits without an error or gets killed while loading:** Docker probably ran out of memory. Increase the memory limit in Docker Desktop, or load fewer organisms (see below).
- **Loading fewer organisms:** to start faster or use less memory, edit `Config.py.template` before building the image and remove the organisms you don't need. Every organism that stays in the list needs a matching `ontology_<taxon id>.zst` file in `maptogo-data`.


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
