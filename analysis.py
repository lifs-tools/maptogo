from dash import Dash, dcc, html, Input, Output, State, callback, exceptions, no_update, MATCH, ALL, callback_context, clientside_callback
import dash_mantine_components as dmc
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import dash_ag_grid as dag
from dash_iconify import DashIconify
import json
import numpy as np
import pandas as pd
import io
import os
from statsmodels.stats.multitest import multipletests
from EnrichmentDataStructure import EnrichmentOntology, current_path
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import LipidParser
import logging
import pathlib
import time
import hashlib
import threading


hash_function = hashlib.new('sha256')
LINK_COLOR = "#2980B9"

logger = logging.getLogger(__name__)
logging.basicConfig(
    level = os.environ.get("LOGLEVEL", "INFO"),
    format = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
)

CC_LINK = html.A(
    "CC BY 4.0",
    href = "https://creativecommons.org/licenses/by/4.0/legalcode#s3a1",
    target = "_blank",
    style = {"color": LINK_COLOR},
)

species = {
    'Homo sapiens': '9606',
    'Mus musculus': '10090',
    # 'Saccharomyces cerevisiae': '4932',
    # 'Escherichia coli': '562',
    # 'Drosophila melanogaster': '7227',
    # 'Rattus norvegicus': '10116',
    # 'Bos taurus': '9913',
    # 'Caenorhabditis elegans': '6239',
    # 'Pseudomonas aeruginosa': '287',
    # 'Arabidopsis thaliana': '3702',
}
INIT_ORGANISM = "9606"

class SessionEntry:
    def __init__(self):
        self.time = time.time()
        self.data = None
        self.data_loaded = False

sessions, examples = {}, {}
xl = pd.ExcelFile(f"{current_path}/Data/examples.xlsx")
for worksheet_name in xl.sheet_names:
    if worksheet_name[0] == "_": continue
    df = xl.parse(worksheet_name)
    worksheet = {
        "bgl": [], "regl": [],
        "bgp": [], "regp": [],
        "bgm": [], "regm": [],
        "title": "", "desc": "", "comp": "", "doi": "", "org": INIT_ORGANISM
    }
    if "BackgroundLipids" in df: worksheet["bgl"] = list(df["BackgroundLipids"].dropna())
    if "RegulatedLipids" in df: worksheet["regl"] = list(df["RegulatedLipids"].dropna())
    if "BackgroundProteins" in df: worksheet["bgp"] = list(df["BackgroundProteins"].dropna())
    if "RegulatedProteins" in df: worksheet["regp"] = list(df["RegulatedProteins"].dropna())
    if "BackgroundMetabolites" in df: worksheet["bgm"] = list(df["BackgroundMetabolites"].dropna())
    if "RegulatedMetabolites" in df: worksheet["regm"] = list(df["RegulatedMetabolites"].dropna())
    if "Title" in df: worksheet["title"] = df.loc[0, "Title"]
    if "Description" in df: worksheet["desc"] = df.loc[0, "Description"]
    if "Comparison" in df: worksheet["comp"] = df.loc[0, "Comparison"]
    if "DOI" in df: worksheet["doi"] = df.loc[0, "DOI"]
    if "Organism" in df: worksheet["org"] = df.loc[0, "Organism"]
    if type(worksheet["title"]) in {float, np.float64} and np.isnan(worksheet["title"]): worksheet["title"] = ""
    if type(worksheet["desc"]) in {float, np.float64} and np.isnan(worksheet["desc"]): worksheet["desc"] = ""
    if type(worksheet["comp"]) in {float, np.float64} and np.isnan(worksheet["comp"]): worksheet["comp"] = ""
    if type(worksheet["doi"]) in {float, np.float64} and np.isnan(worksheet["doi"]): worksheet["doi"] = ""
    if type(worksheet["org"]) in {float, np.float64}:
        if np.isnan(worksheet["org"]): worksheet["org"] = INIT_ORGANISM
        else: worksheet["org"] = str(int(worksheet["org"]))
    examples[worksheet_name] = worksheet

plotly_config = {
    "scrollZoom": True,
    "modeBarButtonsToRemove": [
        "resetScale",
        "toImage",
        "zoom2d",
        "pan2d",
        "select2d",
        "lasso2d",
        "zoomIn2d",
        "zoomOut2d",
        "autoScale2d",
    ],
    "modeBarButtonsToAdd": [],
}


example_options = html.Div(
    dmc.ScrollArea(
        [dmc.Group(
            children = [
                dmc.Checkbox(
                    id = {"type": "checkbox_type", "index": example_name},
                    value = example_name,
                    checked = i == 0,
                    style = {
                        "padding": "10px",
                        "marginBottom": "5px",
                    },
                ),
                dmc.Paper(
                    children = [
                        dmc.Title(
                            example["title"],
                            order = 5,
                            style = {"marginBottom": "5px"},
                        ),
                        (
                            (
                                html.A(
                                    example["desc"],
                                    href = f"https://doi.org/{example['doi']}",
                                    target = "_blank",
                                    style = {"color": LINK_COLOR},
                                )
                            ) if "doi" in example and len(example["doi"]) > 0 else dmc.Text(example["desc"])
                        ),
                        dmc.Text(example["comp"]),
                    ],
                    shadow = "xs",
                    style = {
                        "padding": "10px",
                        "width": "90%",
                        "marginBottom": "15px"
                    },
                ),
            ],
        ) for i, (example_name, example) in enumerate(examples.items())],
        h = 600,
    ),
)

# Create the Dash app
app = Dash("app", update_title = None)
app.title = "GO lipids"

lipid_parser = LipidParser()
enrichment_ontologies = {}
for tax_name, tax_id in species.items():
    logger.info(f"loading {tax_name}")
    enrichment_ontologies[tax_id] = EnrichmentOntology(f"{current_path}/Data/ontology_{tax_id}.gz", lipid_parser = lipid_parser)


#enrichment_ontologies["10090"].set_background(protein_list = ["Q921I1"])
#enrichment_ontologies["10090"].enrichment_analysis(["LPA 18:0"], ["Biological process"])


def get_aggrid_modal(name, molecule):
    return dag.AgGrid(
        id = name,
        columnDefs = [
            {
                'field': "molecule",
                "headerName": molecule,
            },
            {
                'field': "regulated",
                "headerName": "Regulated",
                "maxWidth": 150,
            },
        ],
        rowData = [],
        columnSize = "responsiveSizeToFit",
        defaultColDef={
            "suppressMovable": True,
            "sortable": True,
            "filter": True,
        },
        dashGridOptions={
            "rowSelection": "multiple",
            "suppressMoveWhenRowDragging": True,
            "suppressRowClickSelection": True,
            "alwaysShowVerticalScroll": True,
            "noRowsOverlayComponent": "CustomNoRowsOverlay",
            "noRowsOverlayComponentParams": {
                "message": f"No {molecule.lower()} assigned to this term",
                "fontSize": 12,
            },
        },
    )


def session_timer_trigger(sessions):
    while True:
        time.sleep(10)
        current_time = time.time()
        sessions_to_delete = []
        for session_id, session_data in sessions.items():
            if current_time - session_data.time > 60 * 60: # one hour
                sessions_to_delete.append(session_id)
        for session_id in sessions_to_delete:
            del sessions[session_id]

thread = threading.Thread(target = session_timer_trigger, args = (sessions,), daemon = True)
thread.start()
logger.info("loaded")



def layout():
    current_time = time.time()
    hash_function.update(f"{current_time}t".encode())
    session_id = hash_function.hexdigest()
    sessions[session_id] = SessionEntry()

    return html.Div([
        dcc.Download(id = "download_data"),
        html.Div(session_id, id = "sessionid", style = {"display": "none"}),
        dcc.Loading(
            id = "loading_field",
            children = [html.Div([html.Div(id = "loading_output")])],
            type = "circle",
            fullscreen = True,
            style = {
                "opacity": "70%",
                "z-index": "100",
            },
        ),
        dmc.SimpleGrid([
            dmc.Group([
                dmc.Image(
                    src = f"/assets/golipids.png",
                    style = {"width": "48px"},
                ),
                html.Div(
                    dmc.Title("Gene Ontology (GO) lipidomics enrichment analysis"),
                ),
            ]),
            html.Div([
                html.A(
                    dmc.Text(
                        "Load example datasets",
                    ),
                    id = "load_examples_link",
                    style = {"color": LINK_COLOR, "cursor": "pointer"},
                ),
                dmc.Text("|", style = {"marginLeft": "20px", "marginRight": "20px"}),
                html.A(
                    dmc.Text(
                        "Description & disclaimer",
                    ),
                    id = "disclaimer_link",
                    style = {"color": LINK_COLOR, "cursor": "pointer"},
                )],
                style = {"height": "100%", "display": "flex", "alignItems": "flex-start", "justifyContent": "right"},
            ),
        ], cols = 2),
        html.Div(id = "background_lipids", style = {"display": "none"}),
        html.Div(id = "regulated_lipids", style = {"display": "none"}),
        html.Div(id = "background_proteins", style = {"display": "none"}),
        html.Div(id = "regulated_proteins", style = {"display": "none"}),
        html.Div(id = "background_metabolites", style = {"display": "none"}),
        html.Div(id = "regulated_metabolites", style = {"display": "none"}),
        dmc.Modal(
            title = "Load example datasets",
            id = "load_examples_modal",
            zIndex = 10000,
            children = [
                example_options,
                dmc.Space(h = 50),
                dmc.Group(
                    [
                        dmc.Button(
                            "Load example",
                            color = "blue.4",
                            id = "load_examples_modal_submit_button",
                        ),
                        dmc.Button(
                            "Cancel",
                            color = "red",
                            variant = "outline",
                            id = "load_examples_modal_close_button",
                        ),
                    ],
                    position="right",
                ),
            ],
            size = "60%",
        ),
        dmc.Modal(
            title = "Load",
            id = "term_molecules_modal",
            zIndex = 10000,
            children = [
                dmc.Tabs(
                    [
                        dmc.TabsList(
                            [
                                dmc.Tab(
                                    "Lipids",
                                    value = "lipid_tab_modal",
                                    id = "lipid_tab_modal_tab",
                                ),
                                dmc.Tab(
                                    "Proteins",
                                    value = "protein_tab_modal",
                                    id = "protein_tab_modal_tab",
                                ),
                                dmc.Tab(
                                    "Metabolites",
                                    value = "metabolite_tab_modal",
                                    id = "metabolite_tab_modal_tab",
                                ),
                            ],
                            grow = True,
                        ),
                        dmc.TabsPanel(
                            get_aggrid_modal("term_lipids_modal_grid", "Lipid"),
                            value = "lipid_tab_modal",
                        ),
                        dmc.TabsPanel(
                            get_aggrid_modal("term_proteins_modal_grid", "Protein"),
                            value = "protein_tab_modal",
                        ),
                        dmc.TabsPanel(
                            get_aggrid_modal("term_metabolites_modal_grid", "Metabolite"),
                            value = "metabolite_tab_modal",
                        ),
                    ],
                    id = "term_molecules_modal_tab",
                    color = "blue", # default is blue
                    orientation = "horizontal", # or "vertical"
                    variant = "outline", # or "outline" or "pills"
                    value = "lipid_tab_modal",
                ),
            ],
            size = "40%",
        ),
        dmc.Modal(
            #title = "Load",
            id = "barplot_terms_modal",
            zIndex = 10000,
            children = [
                dcc.Graph(
                    id = "barplot_terms",
                    config = plotly_config,
                ),
            ],
            size = "90%",
        ),
        dmc.Modal(
            [
                dmc.Text("", id = "info_modal_message"),
                dmc.Space(h = 50),
                dmc.Group(
                    dmc.Button(
                        "Close",
                        id = "info_modal_close",
                        color = "blue.4",
                        className = "ml-auto"
                    ),
                    position = "right",
                ),

            ],
            id = "info_modal",
            title = "Information",
            opened = False,
    ),
        dmc.Modal(
            #title = "Description & disclaimer",
            id = "disclaimer_modal",
            zIndex = 10000,
            children = dmc.ScrollArea([
                html.P([
                    dmc.Title("About GO lipids", order = 4),
                    dmc.Text(
                        "Go lipids is the first implementation that enables GO term enrichment analysis on lipids. Starting with a list of identified lipids (at least on the species level) as background and a list with only the (differentially) regulated lipids, each lipid is mapped to proteins involved in metabolic reactions. On this base, GO term analysis will then be performed with a statistical test. Depending on the selected domain (biological process, molecular function, cellular compartment, physical or chemical properties, or metabolic and signaling pathways), the analysis provides a sorted list of GO terms with a p-value for each term.",
                        style = {"textAlign": "justify"},
                    ),
                ]),
                dmc.Text("Responsible people for this database are:"),
                dmc.Text("Dominik Kopczynski: Implementation of both front- and backend"),
                dmc.Text("Cristina Coman: Method validation"),
                dmc.Text("Robert Ahrends: Project leader"),

                html.P([
                    dmc.Title("Disclaimer", order = 4),
                    dmc.Title("1. Liability limitation", order = 5),
                    dmc.Text(
                        "The usage of the downloadable data available on this web page takes place at the users own risk.",
                        style = {"textAlign": "justify"},
                    ),
                ]),

                html.P([
                    dmc.Title("2. External links", order = 5),
                    dmc.Text(
                        "This web page contains links forwarding to external web pages. For all external web pages we disclaim liability. During the linking no statutory violation was evident. As soon as a violation emerges, we delete these links.",
                        style = {"textAlign": "justify"},
                    ),
                ]),

                html.P([
                    dmc.Title("3. Copyright", order = 5),
                    dmc.Text(
                        "The content of this web page is subject to the Austrian copyright law and ancillary copyright. Every inadmissible utilization defined by the Austrian copyright law and ancillary copyright is forbidden. This includes the disposal of the downloadable data. Download and utilization of the available data for researching purposes is explicitly allowed.",
                        style = {"textAlign": "justify"},
                    ),
                ]),

                html.P([
                    dmc.Title("4. Privacy", order = 5),
                    dmc.Text(
                        "For statistical purposes, we record the number of visits to this web page as well as the number of downloads. To ensure an unbiased record, we store a hash value of your IP address for several minutes. Since this method is not reversible, restoring your IP address is not possible, hence this value does not belong to personal-related data. We do not sell or hand out our collected statistical data to third parties.",
                        style = {"textAlign": "justify"},
                    ),
                ]),

                html.P([
                    dmc.Text("GO lipids is using data from third-party databases. All databases are listed with their according licenses or permission:"),
                    dmc.Text([
                        "- ",
                        html.A(
                            "ChEBI",
                            href = "https://www.ebi.ac.uk/chebi/aboutChebiForward.do",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                        ": ",
                        CC_LINK
                    ]),
                    dmc.Text([
                        "- ",
                        html.A(
                            "Gene Ontology (GO)",
                            href = "https://www.geneontology.org/docs/go-citation-policy/",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                        ": ",
                        CC_LINK
                    ]),
                    dmc.Text([
                        "- ",
                        html.A(
                            "LION",
                            href = "https://martijnmolenaar.github.io/lipidontology.com/faq.html",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                        " via ",
                        html.A(
                            "BioPortal",
                            href = "https://www.bioontology.org/terms/",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                        ": freely available for public use",
                    ]),
                    dmc.Text([
                        "- ",
                        html.A(
                            "LIPID MAPS",
                            href = "https://www.lipidmaps.org/databases/lmsd/download",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                        ": ",
                        CC_LINK
                    ]),
                    dmc.Text([
                        "- ",
                        html.A(
                            "Pathbank",
                            href = "https://pathbank.org/about",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                        ": ",
                        html.A(
                            "Open database license",
                            href = "https://opendatacommons.org/licenses/odbl/1-0/",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                    ]),
                    dmc.Text([
                        "- ",
                        html.A(
                            "Rhea",
                            href = "https://www.rhea-db.org/help/license-disclaimer",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                        ": ",
                        CC_LINK
                    ]),
                    dmc.Text([
                        "- ",
                        html.A(
                            "SwissLipids",
                            href = "https://www.swisslipids.org/#/downloads",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                        ": ",
                        CC_LINK
                    ]),
                    dmc.Text([
                        "- ",
                        html.A(
                            "UniProt",
                            href = "https://www.uniprot.org/help/license",
                            target = "_blank",
                            style = {"color": LINK_COLOR},
                        ),
                        ": ",
                        CC_LINK
                    ]),
                ]),
            ], h = 450, offsetScrollbars = True, scrollHideDelay = 0),
            size = "50%",
        ),
        dmc.SimpleGrid(
            cols = 2,
            style = {"marginTop": "5px"},
            children = [
                html.Div([
                    dmc.Tabs(
                        [
                            dmc.TabsList(
                                [
                                    dmc.Tab("Lipids", value="lipid_tab"),
                                    dmc.Tab("Proteins", value="protein_tab"),
                                    dmc.Tab("Metabolites", value="metabolites_tab"),
                                ],
                                grow = True,
                            ),
                            dmc.TabsPanel([
                                dmc.SimpleGrid([
                                    dmc.Title(
                                        "All lipid names in experiment (background)",
                                        order = 5,
                                        style = {"marginTop": "10px"},
                                    ),
                                    dmc.Title(
                                        "All regulated lipid names in experiment",
                                        order = 5,
                                        style = {"marginTop": "10px"},
                                    ),
                                ], cols = 2),
                                dmc.SimpleGrid([
                                    dmc.Textarea(
                                        id = "textarea_all_lipids",
                                        style = {"height": "100%", "display": "inline"},
                                        minRows = 15,
                                    ),
                                    dmc.Textarea(
                                        id = "textarea_regulated_lipids",
                                        minRows = 15,
                                    ),
                                ], cols = 2)],
                                value="lipid_tab",
                            ),
                            dmc.TabsPanel([
                                dmc.SimpleGrid([
                                    dmc.Group([
                                        dmc.Title(
                                            "All protein accessions in experiment (background)",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.ActionIcon(
                                            DashIconify(icon="simple-icons:helix", width = 16),
                                            id = "predefined_bg_proteome",
                                            title = "Select predefined background proteome",
                                            style = {"display": "flex", "alignItems": "flex-end"},
                                        ),
                                    ]),
                                    dmc.Title(
                                        "All regulated protein accession in experiment",
                                        order = 5,
                                        style = {"marginTop": "10px"},
                                    ),
                                ], cols = 2),
                                dmc.SimpleGrid([
                                    dmc.Textarea(
                                        id = "textarea_all_proteins",
                                        style = {"height": "100%", "display": "inline"},
                                        minRows = 15,
                                    ),
                                    dmc.Textarea(
                                        id = "textarea_regulated_proteins",
                                        minRows = 15,
                                    ),
                                ], cols = 2)],
                                value="protein_tab",
                            ),
                            dmc.TabsPanel([
                                dmc.SimpleGrid([
                                    dmc.Title(
                                        "All metabolite ChEBI Ids in experiment (background)",
                                        order = 5,
                                        style = {"marginTop": "10px"},
                                    ),
                                    dmc.Title(
                                        "All regulated metabolite ChEBI Ids in experiment",
                                        order = 5,
                                        style = {"marginTop": "10px"},
                                    ),
                                ], cols = 2),
                                dmc.SimpleGrid([
                                    dmc.Textarea(
                                        id = "textarea_all_metabolites",
                                        style = {"height": "100%", "display": "inline"},
                                        minRows = 15,
                                    ),
                                    dmc.Textarea(
                                        id = "textarea_regulated_metabolites",
                                        minRows = 15,
                                    ),
                                ], cols = 2)],
                                value="metabolites_tab",
                            ),
                        ],
                        color = "blue", # default is blue
                        orientation = "horizontal", # or "vertical"
                        variant = "outline", # or "outline" or "pills"
                        value = "lipid_tab"
                    ),
                    dmc.Title(
                        "Analysis parameters",
                        order = 5,
                        style = {"marginTop": "20px"},
                    ),
                    html.Div(
                        [
                            dmc.SimpleGrid(
                                cols = 2,
                                children = [
                                    dmc.Select(
                                        id = "select_organism",
                                        data = [
                                            {"value": species[key], "label": key} for key in sorted(species.keys())
                                        ],
                                        value = INIT_ORGANISM,
                                        label = "Select organism:",
                                    ),
                                    dmc.MultiSelect(
                                        id = "select_domains",
                                        data = sorted(list(enrichment_ontologies[INIT_ORGANISM].domains)),
                                        value = ["Biological process"],
                                        label = "Select domain(s):",
                                    ),
                                    dmc.Select(
                                        id = "select_test_method",
                                        label = "Select method for p-value correction:",
                                        data=[
                                            {"value": "no", "label": "No correction"},
                                            {"value": "bonferroni", "label": "Bonferroni"},
                                            {"value": "fdr_bh", "label": "Benjamini Hochberg"},
                                        ],
                                        value = "fdr_bh",
                                    ),
                                    dmc.Select(
                                        id = "select_term_regulation",
                                        data = [
                                            {"value": "two-sided", "label": "Term regulated"},
                                            {"value": "less", "label": "Term down-regulated"},
                                            {"value": "greater", "label": "Term up-regulated"},
                                        ],
                                        value = "two-sided",
                                        label = "Term regulation:",
                                    ),
                                    html.Div(
                                        dmc.Switch(
                                            id = "switch_ignore_unparsable_lipids",
                                            checked = False,
                                            label = "Ignore unrecognizable molecules",
                                            style = {"paddingBottom": "8px"},
                                        ),
                                        style = {"height": "100%", "display": "flex", "alignItems": "flex-end", "paddingLeft": "10px"},
                                    ),
                                    html.Div(
                                        dmc.Switch(
                                            id = "switch_ignore_unknown_regulated_lipids",
                                            checked = False,
                                            label = "Ignore regulated molecules that aren't in background",
                                            style = {"paddingBottom": "8px"},
                                        ),
                                        style = {"height": "100%", "display": "flex", "alignItems": "flex-end", "paddingLeft": "10px"},
                                    ),
                                ],
                            ),

                            html.Div([
                                dmc.Group([
                                    dmc.Checkbox(
                                        label = "Use lipids",
                                        id = "checkbox_use_lipids",
                                        checked = False,
                                    ),
                                    dmc.Checkbox(
                                        label = "Use proteins",
                                        id = "checkbox_use_proteins",
                                        checked = False,
                                    ),
                                    dmc.Checkbox(
                                        label = "Use metabolites",
                                        id = "checkbox_use_metabolites",
                                        checked = False,
                                    ),
                                ]),
                                dmc.Button(
                                    "Run enrichment analysis",
                                    id = "button_run_enrichment",
                                    style = {"margin-left": "auto"},
                                ),
                            ], style={"display": "flex", "width": "100%", "marginTop": "20px"}),
                        ],
                        style = {"border": "1px solid #dddddd", "radius": "10px", "padding": "10px"},
                    ),
                ]),
                html.Div([
                    dmc.SimpleGrid([
                        dmc.Title(
                            "Results",
                            order = 5,
                            style = {"marginTop": "10px"},
                        ),
                        html.Div(
                            html.Span(
                                dmc.Group([
                                    dmc.ActionIcon(
                                        DashIconify(icon="tdesign:chart-column-filled", width = 20),
                                        id = "chart_results",
                                        title = "Show chart",
                                        disabled = True,
                                    ),
                                    dmc.ActionIcon(
                                        DashIconify(icon="material-symbols:download-rounded", width = 20),
                                        id = "icon_download_results",
                                        title = "Download table",
                                        disabled = True,
                                    ),
                                ], spacing = 0),
                            ),
                            style = {"height": "100%", "display": "flex", "alignItems": "flex-end", "justifyContent": "right"},
                        ),
                    ], cols = 2),
                    dag.AgGrid(
                        id = "graph_enrichment_results",
                        columnDefs = [
                            {
                                'field': "domain",
                                "headerName": "Domain",
                                "maxWidth": 220,
                                "checkboxSelection": True,
                            },
                            {
                                'field': "termid",
                                "headerName": "Term ID",
                                "cellRenderer": "TermIDRenderer",
                                "maxWidth": 150,
                            },
                            {
                                'field': "term",
                                "headerName": "Term",
                                "cellRenderer": "TermRenderer",
                            },
                            {
                                'field': "pvalue",
                                "headerName": "pValue",
                                "width": 150,
                            },
                        ],
                        rowData = [],
                        columnSize = "responsiveSizeToFit",
                        defaultColDef={
                            "suppressMovable": True,
                            "resizable": True,
                            "sortable": True,
                            "filter": True,
                            "headerCheckboxSelection": {
                                "function": "params.column == params.columnApi.getAllDisplayedColumns()[0]"
                            },
                        },
                        style = {
                            "height": "100%",
                        },
                        getRowId = "params.data.termid",
                        dashGridOptions={
                            "rowSelection": "multiple",
                            "suppressMoveWhenRowDragging": True,
                            "suppressRowClickSelection": True,
                            "alwaysShowVerticalScroll": True,
                            "noRowsOverlayComponent": "CustomNoRowsOverlay",
                            "noRowsOverlayComponentParams": {
                                "message": "Run the enrichment analysis to obtain results",
                                "fontSize": 12,
                            },
                        },
                    ),
                ]),
            ],
        ),
        html.Div(style = {"height": "500px"}),
    ])

app.layout = layout



@callback(
    Output("loading_output", "children", allow_duplicate = True),
    Output("graph_enrichment_results", "rowData", allow_duplicate = True),
    #Output("alert_enrichment", "hide", allow_duplicate = True),
    #Output("alert_enrichment", "children", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("background_lipids", "children", allow_duplicate = True),
    Output("regulated_lipids", "children", allow_duplicate = True),
    Output("background_proteins", "children", allow_duplicate = True),
    Output("regulated_proteins", "children", allow_duplicate = True),
    Output("background_metabolites", "children", allow_duplicate = True),
    Output("regulated_metabolites", "children", allow_duplicate = True),
    Input("button_run_enrichment", "n_clicks"),
    State("textarea_all_lipids", "value"),
    State("textarea_regulated_lipids", "value"),
    State("textarea_all_proteins", "value"),
    State("textarea_regulated_proteins", "value"),
    State("textarea_all_metabolites", "value"),
    State("textarea_regulated_metabolites", "value"),
    State("select_organism", "value"),
    State("select_domains", "value"),
    State("switch_ignore_unparsable_lipids", "checked"),
    State("select_test_method", "value"),
    State("switch_ignore_unknown_regulated_lipids", "checked"),
    State("select_term_regulation", "value"),
    State("sessionid", "children"),
    State("checkbox_use_lipids", "checked"),
    State("checkbox_use_proteins", "checked"),
    State("checkbox_use_metabolites", "checked"),
    prevent_initial_call = True,
)
def run_enrichment(
    n_clicks,
    all_lipids_list,
    regulated_lipids_list,
    all_proteins_list,
    regulated_proteins_list,
    all_metabolites_list,
    regulated_metabolites_list,
    organism,
    domains,
    ignore_unrecognizable_molecules,
    correction_method,
    ignore_unknown,
    term_regulation,
    session_id,
    with_lipids,
    with_proteins,
    with_metabolites,
):
    if session_id not in sessions:
        raise exceptions.PreventUpdate

    do_activate_alert = True

    if not with_lipids and not with_proteins and not with_metabolites:
        return "", [], do_activate_alert, "No omics data is inserted.", "", "", "", "", "", ""

    ontology = enrichment_ontologies[organism]

    target_set = set()
    lipidome, regulated_lipids = {}, set()
    proteome, regulated_proteins = set(), set()
    metabolome, regulated_metabolites = set(), set()
    if with_lipids:
        if len(all_lipids_list) == 0:
            return "", [], do_activate_alert, "Please paste lipid names into the first text area.", "", "", "", "", "", ""

        if len(regulated_lipids_list) == 0:
            return "", [], do_activate_alert, "Please paste lipid names into the second text area.", "", "", "", "", "", ""

        for lipid_name in all_lipids_list.split("\n"):
            if len(lipid_name) == 0: continue
            try:
                lipid = lipid_parser.parse(lipid_name)
                if lipid.lipid.info.level.value > LipidLevel.MOLECULAR_SPECIES.value:
                    lipid.lipid.info.level = LipidLevel.MOLECULAR_SPECIES
                lipid.sort_fatty_acyl_chains()
                lipidome[lipid_name] = lipid

            except Exception as e:
                if not ignore_unrecognizable_molecules:
                    return "", [], do_activate_alert, f"Lipid name '{lipid_name}' unrecognizable in background list! Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", ""

        for lipid_name in regulated_lipids_list.split("\n"):
            if len(lipid_name) == 0: continue
            if lipid_name not in lipidome:
                if ignore_unknown: continue
                return "", [], do_activate_alert, f"The regulated lipid '{lipid_name}' does not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", "", "", "", "", "", ""


            try:
                regulated_lipids.add(lipidome[lipid_name].get_lipid_string())
            except Exception as e:
                print(e)
                if not ignore_unrecognizable_molecules:
                    return "", [], do_activate_alert, f"Lipid name '{lipid_name}' unrecognizable in regulated list! Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", ""

        if len(lipidome) == 0:
            return "", [], do_activate_alert, "No background lipid left after lipid recognition.", "", "", "", "", "", ""

        if len(regulated_lipids) == 0:
            return "", [], do_activate_alert, "No regulated lipid left after lipid recognition.", "", "", "", "", "", ""

        if len(regulated_lipids) > len(lipidome):
            return "", [], do_activate_alert, "Length of regulated lipid list must be smaller than background list.", "", "", "", "", "", ""

        left_lipids = regulated_lipids - set([v.get_lipid_string() for k, v in lipidome.items()])
        if len(left_lipids) > 0:
            if ignore_unknown:
                for lipid_name in left_lipids:
                    del lipidome[lipid_name]
            else:
                return "", [], do_activate_alert, "The regulated lipid" + (' ' if len(left_lipids) == 1 else 's ') + "'" + "', '".join(left_lipids) + ("' does" if len(left_lipids) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", "", "", "", "", "", ""

        target_set |= regulated_lipids

    if with_proteins:
        if len(all_proteins_list) == 0:
            return "", [], do_activate_alert, "Please paste protein accession into the first text area.", "", "", "", "", "", ""

        if len(regulated_proteins_list) == 0:
            return "", [], do_activate_alert, "Please paste protein accessions into the second text area.", "", "", "", "", "", ""

        proteome = set(protein for protein in all_proteins_list.split("\n") if len(protein) > 0)
        regulated_proteins = set(protein for protein in regulated_proteins_list.split("\n") if len(protein) > 0)

        left_proteins = proteome - ontology.clean_protein_ids
        if len(left_proteins) > 0:
            if ignore_unrecognizable_molecules:
                proteome -= left_proteins
            else:
                return "", [], do_activate_alert, "The protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' is" if len(left_proteins) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", ""

        left_proteins = regulated_proteins - ontology.clean_protein_ids
        if len(left_proteins) > 0:
            if ignore_unrecognizable_molecules:
                regulated_proteins -= left_proteins
            else:
                return "", [], do_activate_alert, "The protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' is" if len(left_proteins) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", ""

        left_proteins = regulated_proteins - proteome
        if len(left_proteins) > 0:
            if ignore_unknown:
                proteome -= left_proteins
            else:
                return "", [], do_activate_alert, "The regulated protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' does" if len(left_proteins) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", "", "", "", "", "", ""

        if len(proteome) == 0:
            return "", [], do_activate_alert, "No background protein left after protein recognition.", "", "", "", "", "", ""

        if len(regulated_proteins) == 0:
            return "", [], do_activate_alert, "No regulated protein left after protein recognition.", "", "", "", "", "", ""

        if len(regulated_proteins) > len(proteome):
            return "", [], do_activate_alert, "Length of regulated protein list must be smaller than background list.", "", "", "", "", "", ""

        proteome = set([f"UNIPROT:{protein}" for protein in proteome])
        regulated_proteins = set([f"UNIPROT:{protein}" for protein in regulated_proteins])
        target_set |= regulated_proteins

    if with_metabolites:
        if len(all_metabolites_list) == 0:
            return "", [], do_activate_alert, "Please paste metabolite ChEBI Ids into the first text area.", "", "", "", "", "", ""

        if len(regulated_metabolites_list) == 0:
            return "", [], do_activate_alert, "Please paste metabolite ChEBI Ids into the second text area.", "", "", "", "", "", ""

        metabolome = set(metabolite for metabolite in all_metabolites_list.split("\n") if len(metabolite) > 0)
        regulated_metabolites = set(metabolite for metabolite in regulated_metabolites_list.split("\n") if len(metabolite) > 0)

        left_metabolites = metabolome - ontology.clean_metabolite_ids - ontology.metabolites.keys()
        if len(left_metabolites) > 0:
            if ignore_unrecognizable_molecules:
                metabolome -= left_metabolites
            else:
                return "", [], do_activate_alert, "The metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' is" if len(left_metabolites) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", ""

        left_metabolites = regulated_metabolites - ontology.clean_metabolite_ids - ontology.metabolites.keys()
        if len(left_metabolites) > 0:
            if ignore_unrecognizable_molecules:
                regulated_metabolites -= left_metabolites
            else:
                return "", [], do_activate_alert, "The metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' is" if len(left_metabolites) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", ""

        left_metabolites = regulated_metabolites - metabolome
        if len(left_metabolites) > 0:
            if ignore_unknown:
                metabolome -= left_metabolites
            else:
                return "", [], do_activate_alert, "The regulated metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' does" if len(left_metabolites) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", "", "", "", "", "", ""

        if len(metabolome) == 0:
            return "", [], do_activate_alert, "No background metabolite left after metabolite recognition.", "", "", "", "", "", ""

        if len(regulated_metabolites) == 0:
            return "", [], do_activate_alert, "No regulated metabolite left after metabolite recognition.", "", "", "", "", "", ""

        if len(regulated_metabolites) > len(metabolome):
            return "", [], do_activate_alert, "Length of regulated metabolite list must be smaller than background list.", "", "", "", "", "", ""

        metabolome = set([f"CHEBI:{metabolite}" if type(metabolite) == int or metabolite[:6] != "CHEBI:" else metabolite for metabolite in metabolome])
        regulated_metabolites = set([f"CHEBI:{metabolite}" if type(metabolite) == int or metabolite[:6] != "CHEBI:" else metabolite for metabolite in regulated_metabolites])
        target_set |= regulated_metabolites

    if len(domains) == 0:
        return "", [], do_activate_alert, "No domain(s) selected.", "", "", "", "", "", ""

    if sessions[session_id].data_loaded == False:
        ontology.set_background(lipid_list = lipidome, protein_list = proteome, metabolite_list = metabolome)
        sessions[session_id].data_loaded = True

    results = ontology.enrichment_analysis(target_set, domains, term_regulation)

    if correction_method != "no" and len(results) > 0:
        pvalues = [r.pvalue for r in results]
        pvalues = multipletests(pvalues, method = correction_method)[1]
        for pvalue, r in zip(pvalues, results): r.pvalue = pvalue

    results.sort(key = lambda row: row.pvalue)
    data = [
        {
            "domain": result.term.domain,
            "term": result.term.name,
            "termid": result.term.term_id,
            "pvalue": result.pvalue
        } for result in results
    ]

    sessions[session_id].time = time.time()
    sessions[session_id].data = {result.term.term_id: result for result in results}

    return (
        "",
        data,
        not do_activate_alert,
        "",
        "|".join(lipidome),
        "|".join(regulated_lipids),
        "|".join(proteome) if with_proteins else "",
        "|".join(regulated_proteins) if with_proteins else "",
        "|".join(metabolome) if with_metabolites else "",
        "|".join(regulated_metabolites) if with_metabolites else "",
    )



@callback(
    Output("textarea_all_lipids", "id", allow_duplicate = True),
    Input("textarea_all_lipids", "value"),
    Input("textarea_regulated_lipids", "value"),
    Input("textarea_all_proteins", "value"),
    Input("textarea_regulated_proteins", "value"),
    Input("textarea_all_metabolites", "value"),
    Input("textarea_regulated_metabolites", "value"),
    Input("checkbox_use_lipids", "checked"),
    Input("checkbox_use_proteins", "checked"),
    Input("checkbox_use_metabolites", "checked"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def update_background(
    all_lipids_list,
    regulated_lipids_list,
    all_proteins_list,
    regulated_proteins_list,
    all_metabolites_list,
    regulated_metabolites_list,
    with_lipids,
    with_proteins,
    with_metabolites,
    session_id,
):
    if session_id not in sessions:
        raise exceptions.PreventUpdate

    sessions[session_id].data_loaded = False

    raise exceptions.PreventUpdate



@callback(
    Output("loading_output", "children", allow_duplicate = True),
    Output("download_data", "data", allow_duplicate = True),
    Input("icon_download_results", "n_clicks"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    State("background_lipids", "children"),
    State("regulated_lipids", "children"),
    State("background_proteins", "children"),
    State("regulated_proteins", "children"),
    State("background_metabolites", "children"),
    State("regulated_metabolites", "children"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def download_table(
    _,
    graph_enrichment_results,
    selected_rows,
    background_lipids,
    regulated_lipids,
    background_proteins,
    regulated_proteins,
    background_metabolites,
    regulated_metabolites,
    session_id,
):
    if session_id not in sessions:
        raise exceptions.PreventUpdate

    domains = []
    term_ids = []
    terms = []
    pvalues = []

    selected_term_ids = {row["termid"] for row in selected_rows}
    if len(graph_enrichment_results) == 0 or len(selected_term_ids) == 0:
        raise exceptions.PreventUpdate

    for row in graph_enrichment_results:
        term_id = row["termid"]
        if term_id not in selected_term_ids: continue
        domains.append(row["domain"])
        term_ids.append(term_id)
        terms.append(row["term"])
        pvalues.append(row["pvalue"])

    df = pd.DataFrame({"Domain": domains, "Term ID": term_ids, "Term": terms, "pValue": pvalues})
    data = sessions[session_id].data

    with_lipids = len(background_lipids) > 0 or len(regulated_lipids) > 0
    with_proteins = len(background_proteins) > 0 or len(regulated_proteins) > 0
    with_metabolites = len(background_metabolites) > 0 or len(regulated_metabolites) > 0

    if with_lipids:
        background_lipids = sorted(background_lipids.split("|"))
        regulated_lipids = sorted(regulated_lipids.split("|"))
        regulated_lipid_set = set(regulated_lipids)
        associated_lipids = {}

    if with_proteins:
        background_proteins = sorted(background_proteins.split("|"))
        regulated_proteins = sorted(regulated_proteins.split("|"))
        regulated_protein_set = set(regulated_proteins)
        associated_proteins = {}

    if with_metabolites:
        background_metabolites = sorted(background_metabolites.split("|"))
        regulated_metabolites = sorted(regulated_metabolites.split("|"))
        regulated_metabolite_set = set(regulated_metabolites)
        associated_metabolites = {}

    for term_id, term in zip(term_ids, terms):
        if term_id not in data: continue
        molecules = data[term_id].term.term_paths.keys()

        if with_lipids:
            lipids = sorted(list(molecules & set(background_lipids)))
            associated_lipids[term] = lipids
            associated_lipids[f"{term}, regulated"] = ["X" if lipid in regulated_lipid_set else "" for lipid in lipids]

        if with_proteins:
            proteins = sorted(list(molecules & set(background_proteins)))
            associated_proteins[term] = [protein.replace("UNIPROT:", "") for protein in proteins]
            associated_proteins[f"{term}, regulated"] = ["X" if protein in regulated_protein_set else "" for protein in proteins]

        if with_metabolites:
            metabolites = sorted(list(molecules & set(background_metabolites)))
            associated_metabolites[term] = [metabolite for metabolite in metabolites]
            associated_metabolites[f"{term}, regulated"] = ["X" if metabolite in regulated_metabolite_set else "" for metabolite in metabolites]

    output = io.BytesIO()
    writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
    df.to_excel(writer, sheet_name = "GO multiomics results", index = False)
    if with_lipids:
        pd.DataFrame.from_dict(associated_lipids, orient = "index").transpose().to_excel(writer, sheet_name = "Associated lipids", index = False)
        pd.DataFrame({"LipidName": background_lipids}).to_excel(writer, sheet_name = "Background lipids", index = False)
        pd.DataFrame({"LipidName": regulated_lipids}).to_excel(writer, sheet_name = "Regulated lipids", index = False)

    if with_proteins:
        background_proteins = [protein.replace("UNIPROT:", "") for protein in background_proteins]
        regulated_proteins = [protein.replace("UNIPROT:", "") for protein in regulated_proteins]
        pd.DataFrame.from_dict(associated_proteins, orient = "index").transpose().to_excel(writer, sheet_name = "Associated proteins", index = False)
        pd.DataFrame({"Accession": background_proteins}).to_excel(writer, sheet_name = "Background proteins", index = False)
        pd.DataFrame({"Accession": regulated_proteins}).to_excel(writer, sheet_name = "Regulated proteins", index = False)

    if with_metabolites:
        background_metabolites = [metabolite for metabolite in background_metabolites]
        regulated_metabolites = [metabolite for metabolite in regulated_metabolites]
        pd.DataFrame.from_dict(associated_metabolites, orient = "index").transpose().to_excel(writer, sheet_name = "Associated metabolites", index = False)
        pd.DataFrame({"ChEBI": background_metabolites}).to_excel(writer, sheet_name = "Background metabolites", index = False)
        pd.DataFrame({"ChEBI": regulated_metabolites}).to_excel(writer, sheet_name = "Regulated metabolites", index = False)
    writer._save()

    return "", dcc.send_bytes(output.getvalue(), "GO_multiomics_results.xlsx")



@callback(
    Output("load_examples_modal", "opened", allow_duplicate = True),
    Input("load_examples_link", "n_clicks"),
    prevent_initial_call = True,
)
def open_load_examples_modal(n_clicks):
    return True



@callback(
    Output("load_examples_modal", "opened", allow_duplicate = True),
    Output("textarea_all_lipids", "value", allow_duplicate = True),
    Output("textarea_regulated_lipids", "value", allow_duplicate = True),
    Output("textarea_all_proteins", "value", allow_duplicate = True),
    Output("textarea_regulated_proteins", "value", allow_duplicate = True),
    Output("textarea_all_metabolites", "value", allow_duplicate = True),
    Output("textarea_regulated_metabolites", "value", allow_duplicate = True),
    Output("select_organism", "value", allow_duplicate = True),
    Output("checkbox_use_lipids", "checked", allow_duplicate = True),
    Output("checkbox_use_proteins", "checked", allow_duplicate = True),
    Output("checkbox_use_metabolites", "checked", allow_duplicate = True),
    Input("load_examples_modal_submit_button", "n_clicks"),
    State({"type": "checkbox_type", "index": ALL}, "checked"),
    State({"type": "checkbox_type", "index": ALL}, "id"),
    prevent_initial_call = True,
)
def submit_load_examples_modal(n_clicks, checked, checkbox_ids):
    index = ""
    for c, checkbox_id in zip(checked, checkbox_ids):
        if c:
            index = checkbox_id["index"]
            break

    return (
        False,
        "\n".join(examples[index]["bgl"]),
        "\n".join(examples[index]["regl"]),
        "\n".join(examples[index]["bgp"]),
        "\n".join(examples[index]["regp"]),
        "\n".join(examples[index]["bgm"]),
        "\n".join(examples[index]["regm"]),
        examples[index]["org"],
        len(examples[index]["bgl"]) > 0 or len(examples[index]["regl"]) > 0,
        len(examples[index]["bgp"]) > 0 or len(examples[index]["regp"]) > 0,
        len(examples[index]["bgm"]) > 0 or len(examples[index]["regm"]) > 0,
    )



@callback(
    Output("load_examples_modal", "opened", allow_duplicate = True),
    Input("load_examples_modal_close_button", "n_clicks"),
    prevent_initial_call = True,
)
def close_load_examples_modal(n_clicks):
    return False



@callback(
    Output({"type": "checkbox_type", "index": ALL}, "checked", allow_duplicate = True),
    Input({"type": "checkbox_type", "index": ALL}, "checked"),
    State({"type": "checkbox_type", "index": ALL}, "id"),
    prevent_initial_call = True,
)
def checkbox_checks(_, checkbox_ids):
    index = json.loads(callback_context.triggered[0]["prop_id"].split(".")[0])["index"]
    return [checkbox_id["index"] == index for checkbox_id in checkbox_ids]



@callback(
    Output("disclaimer_modal", "opened", allow_duplicate = True),
    Input("disclaimer_link", "n_clicks"),
    prevent_initial_call = True,
)
def disclaimer_clicked(n_clicks):
    return True




@callback(
    Output("chart_results", "disabled", allow_duplicate = True),
    Output("icon_download_results", "disabled", allow_duplicate = True),
    Input("graph_enrichment_results", "selectedRows"),
    prevent_initial_call = True,
)
def update_action_icons(selected_rows):
    return len(selected_rows) == 0, len(selected_rows) == 0



@callback(
    Output("term_molecules_modal", "opened", allow_duplicate = True),
    Output("term_molecules_modal", "title", allow_duplicate = True),
    Output("term_lipids_modal_grid", "rowData", allow_duplicate = True),
    Output("term_proteins_modal_grid", "rowData", allow_duplicate = True),
    Output("term_metabolites_modal_grid", "rowData", allow_duplicate = True),
    Output("lipid_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("protein_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("metabolite_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("term_molecules_modal_tab", "value", allow_duplicate = True),
    Input("graph_enrichment_results", "cellRendererData"),
    State("sessionid", "children"),
    State("background_lipids", "children"),
    State("regulated_lipids", "children"),
    State("background_proteins", "children"),
    State("regulated_proteins", "children"),
    State("background_metabolites", "children"),
    State("regulated_metabolites", "children"),
    prevent_initial_call = True,
)
def open_term_window(
    row_data,
    session_id,
    background_lipids,
    regulated_lipids,
    background_proteins,
    regulated_proteins,
    background_metabolites,
    regulated_metabolites,
):
    if session_id not in sessions or "rowId" not in row_data:
        raise exceptions.PreventUpdate

    term_id = row_data["rowId"]
    if term_id not in sessions[session_id].data:
        raise exceptions.PreventUpdate

    with_lipids = len(background_lipids) > 0 or len(regulated_lipids) > 0
    with_proteins = len(background_proteins) > 0 or len(regulated_proteins) > 0
    with_metabolites = len(background_metabolites) > 0 or len(regulated_metabolites) > 0

    lipid_table = []
    protein_table = []
    metabolite_table = []

    result = sessions[session_id].data[term_id]
    molecules = sorted(list(result.term.term_paths.keys()))

    if with_lipids:
        background_lipids = set(background_lipids.split("|"))
        regulated_lipids = set(regulated_lipids.split("|"))
        lipid_table = [{"molecule": molecule, "regulated": ("X" if molecule in regulated_lipids else "")} for molecule in molecules if molecule in background_lipids]

    if with_proteins:
        background_proteins = set(background_proteins.split("|"))
        regulated_proteins = set(regulated_proteins.split("|"))
        protein_table = [{"molecule": molecule.replace("UNIPROT:", ""), "regulated": ("X" if molecule in regulated_proteins else "")} for molecule in molecules if molecule in background_proteins]

    if with_metabolites:
        background_metabolites = set(background_metabolites.split("|"))
        regulated_metabolites = set(regulated_metabolites.split("|"))
        metabolite_table = [{"molecule": molecule, "regulated": ("X" if molecule in regulated_metabolites else "")} for molecule in molecules if molecule in background_metabolites]

    # for term_id, term_path in result.term.term_paths.items():
    #     print(term_id, term_path)

    tab_value = "lipid_tab_modal"
    if not with_lipids:
        tab_value = "protein_tab_modal"
        if not with_proteins:
            tab_value = "metabolite_tab_modal"

    return (
        True,
        f"Molecules for '{result.term.name}'",
        lipid_table,
        protein_table,
        metabolite_table,
        not with_lipids,
        not with_proteins,
        not with_metabolites,
        tab_value,
    )



@callback(
    Output("barplot_terms_modal", "opened", allow_duplicate = True),
    Output("barplot_terms", "figure", allow_duplicate = True),
    Input("chart_results", "n_clicks"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    prevent_initial_call = True,
)
def open_barplot(n_clicks, row_data, selected_rows):
    fig = go.Figure()
    selected_term_ids = set([row["termid"] for row in selected_rows])
    pvalues = -np.log10([row["pvalue"] for row in row_data if row["termid"] in selected_term_ids])
    categories = [row["term"] for row in row_data if row["termid"] in selected_term_ids]
    # Add bar trace with red color
    fig.add_trace(go.Bar(
        x = pvalues,
        y = categories,
        orientation = "h",  # Horizontal bars
        marker = dict(color = "#1f77b4")  # Set bars to red
    ))

    # Layout customization
    fig.update_layout(
        xaxis = dict(title = "-log<sub>10</sub>(p-value)", fixedrange = True),
        yaxis = dict(title = "Term", categoryarray = categories[::-1]),
        template = "plotly_white",
        dragmode = "pan",
        margin = dict(l = 20, r = 20, t = 20, b = 20),
    )
    return True, fig



clientside_callback(
    """
    function (graph_config) {
        graph_config.modeBarButtonsToAdd.push({
            name: 'download_png',
            title: 'Download as png',
            icon: {
                width: 24,
                height: 24,
                path: 'M22.71 6.29a1 1 0 0 0-1.42 0L20 7.59V2a1 1 0 0 0-2 0v5.59l-1.29-1.3a1 1 0 0 0-1.42 1.42l3 3a1 1 0 0 0 .33.21a.94.94 0 0 0 .76 0a1 1 0 0 0 .33-.21l3-3a1 1 0 0 0 0-1.42M19 13a1 1 0 0 0-1 1v.38l-1.48-1.48a2.79 2.79 0 0 0-3.93 0l-.7.7l-2.48-2.48a2.85 2.85 0 0 0-3.93 0L4 12.6V7a1 1 0 0 1 1-1h8a1 1 0 0 0 0-2H5a3 3 0 0 0-3 3v12a3 3 0 0 0 3 3h12a3 3 0 0 0 3-3v-5a1 1 0 0 0-1-1M5 20a1 1 0 0 1-1-1v-3.57l2.9-2.9a.79.79 0 0 1 1.09 0l3.17 3.17l4.3 4.3Zm13-1a.9.9 0 0 1-.18.53L13.31 15l.7-.7a.77.77 0 0 1 1.1 0L18 17.21Z',
            },
            click: function(gd) {
                Plotly.downloadImage(gd, {format: "png", filename: "barplot"});
            }
        });

        graph_config.modeBarButtonsToAdd.push({
            name: 'download_svg',
            title: 'Download as svg',
            icon: {
                width: 24,
                height: 24,
                path: 'M22.71 6.29a1 1 0 0 0-1.42 0L20 7.59V2a1 1 0 0 0-2 0v5.59l-1.29-1.3a1 1 0 0 0-1.42 1.42l3 3a1 1 0 0 0 .33.21a.94.94 0 0 0 .76 0a1 1 0 0 0 .33-.21l3-3a1 1 0 0 0 0-1.42M19 13a1 1 0 0 0-1 1v.38l-1.48-1.48a2.79 2.79 0 0 0-3.93 0l-.7.7l-2.48-2.48a2.85 2.85 0 0 0-3.93 0L4 12.6V7a1 1 0 0 1 1-1h8a1 1 0 0 0 0-2H5a3 3 0 0 0-3 3v12a3 3 0 0 0 3 3h12a3 3 0 0 0 3-3v-5a1 1 0 0 0-1-1M5 20a1 1 0 0 1-1-1v-3.57l2.9-2.9a.79.79 0 0 1 1.09 0l3.17 3.17l4.3 4.3Zm13-1a.9.9 0 0 1-.18.53L13.31 15l.7-.7a.77.77 0 0 1 1.1 0L18 17.21Z',
            },
            click: function(gd) {
                Plotly.downloadImage(gd, {format: "svg", filename: "barplot"});
            }
        });
        return window.dash_clientside.no_update;
    }
    """,
    Output("barplot_terms", "config"),
    Input("barplot_terms", "config"),
)



@callback(
    Output("info_modal", "opened", allow_duplicate = True),
    Input("info_modal_close", "n_clicks"),
    prevent_initial_call = True
)
def close_info_modal(close_clicks):
    return False



@callback(
    Output("textarea_all_proteins", "value", allow_duplicate = True),
    Input("predefined_bg_proteome", "n_clicks"),
    State("select_organism", "value"),
    prevent_initial_call = True
)
def select_predefined_modal(n_clicks, selected_organism):
    if selected_organism not in enrichment_ontologies:
        raise exceptions.PreventUpdate

    ontology = enrichment_ontologies[selected_organism]
    return "\n".join([protein.replace("UNIPROT:", "") for protein in ontology.proteins])


