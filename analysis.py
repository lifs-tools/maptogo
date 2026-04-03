"""
MIT License

Copyright (c) 2026 Dominik Kopczynski  -  dominik.kopczynski {at} univie.ac.at
                   Cristina Coman  -  cristina.coman {at} univie.ac.at
                   Nils Hoffmann  -  n.hoffmann {at} fz-juelich.de
                   Robert Ahrends  -  robert.ahrends {at} univie.ac.at

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import logging
import os
import uuid

logger = logging.getLogger(__name__)
logging.basicConfig(
    level = os.environ.get("LOGLEVEL", "INFO"),
    format = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
)
logger.info("Started enrichment server")
logging.getLogger("flask_swagger_ui").setLevel(logging.ERROR)
logging.getLogger("swagger_ui").setLevel(logging.ERROR)

DEFAULT_VERSION_NUMBER = "1.0.0"

import dash
from dash import Dash, dcc, html, Input, Output, State, callback, exceptions, no_update, MATCH, ALL, callback_context, clientside_callback, dash_table, ctx
from flask import request, session as uuid_session
import dash_mantine_components as dmc
import plotly.graph_objs as go
import dash_ag_grid as dag
from dash_iconify import DashIconify
import json
import numpy as np
import pandas as pd
import io
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import squareform
from EnrichmentDataStructure import *
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidAdduct import LipidAdduct
import pathlib
import time
import freetype
from math import ceil
from flask_restx import Api, Resource, fields
import traceback
import requests
import threading
from datetime import datetime
from collections import defaultdict, deque
import configparser
import gc


INIT_ORGANISM = "NCBITaxon:10090"
APPLICATION_SHORT_TITLE = "MAPtoGO"
APPLICATION_TITLE = f"{APPLICATION_SHORT_TITLE} - Multiomics Analysis Platform towards Gene Ontology"
TEXTFIELD_HEIGHT = "350px"

if os.path.exists(f"{current_path}/maptogo.ini"):
    ini_config = configparser.ConfigParser()
    ini_config.read(f"{current_path}/maptogo.ini")
else:
    ini_config = {
        "git": {
            "commit": "-",
            "branch": "-",
            "timestamp": "-",
            "version": DEFAULT_VERSION_NUMBER,
        }
    }

class SessionEntry:
    def __init__(self):
        self.time = time.time()
        self.data = None
        self.results = []
        self.data_loaded = False
        self.search_terms = {}
        self.all_parent_nodes = {}
        self.num_background = 0
        self.use_bounded_fatty_acyls = False
        self.ontology = None
        self.domains = None
        self.background_lipids = None
        self.regulated_lipids = None
        self.background_proteins = None
        self.regulated_proteins = None
        self.background_metabolites = None
        self.regulated_metabolites = None
        self.background_transcripts = None
        self.regulated_transcripts = None
        self.background_list = []
        self.min_pvalue = "0.0001"
        self.max_pvalue = "0.05"
        self.ui = {}
        self.sankey_data = None
        self.separate_updown_switch = False



def get_git_info():
    return f"Built on {ini_config['git']['timestamp']} from commit {ini_config['git']['commit']} on branch {ini_config['git']['branch']}"


def get_latest_tag():
    return ini_config['git']['version']



def graph_enrichment_results_def(separate_updown_switch, multiple_domains, qvalue = True):
    return [
        {
            'field': "termid",
            "headerName": "Term ID",
            "cellRenderer": "TermIDRenderer",
            "maxWidth": 160,
            "checkboxSelection": True,
        },
        {
            'field': "domain",
            "headerName": "Domain",
            "maxWidth": 160,
            "hide": not multiple_domains,
        },
        {
            'field': "term",
            "headerName": "Term",
            "cellRenderer": "TermRenderer",
        },
        {
            'field': "direction",
            "headerName": "Direction",
            "hide": not separate_updown_switch,
            "width": 60,
            "headerTooltip": "Up: term enhanced; Down: term suppressed; Both: term modulated / adjusting",
        },
        {
            'field': "count",
            "headerName": "Count",
            "width": 110,
            "headerTooltip": "Number of regulated associated molecules (Expected number of regulated associated molecules) / Number of all associated molecules",
        },
        {
            'field': "pvalue",
            "headerName": "q-value" if qvalue else "p-value",
            "width": 80,
            "valueFormatter": {"function": "params.value != null ? Number(params.value).toPrecision(5) : ''"},
            "filter": "agNumberColumnFilter",
            "headerTooltip": "The p-value is the statistical significance, the q-value is the adjusted p-value after multiple testing correction",
        },
        {
            'field': "log_odds_ratio",
            "headerName": "Log odds",
            "width": 80,
            "valueFormatter": {} if separate_updown_switch else {"function": "params.value != null ? Number(params.value).toPrecision(5) : ''"},
            "filter": "" if separate_updown_switch else "agNumberColumnFilter",
            "headerTooltip": "The 'log odds ratio' determines how strong the association is (effect size). Positive: enriched; Zero: no enrichment; Negative: depleted.",
        },
    ]



correction_list = [
    {"value": "no", "label": "No correction"},
    {"value": "bonferroni", "label": "Bonferroni"},
    {"value": "fdr_bh", "label": "Benjamini/Hochberg (default)"},
    {"value": "fdr_by", "label": "Benjamini/Yekutieli"},
    {"value": "holm-sidak", "label": "Step down method using Sidak adjustment"},
    {"value": "holm", "label": "Step-down method using Bonferroni adjustments"},
    {"value": "simes-hochberg", "label": "Step-up method (independent)"},
    {"value": "fdr_tsbh", "label": "Two stage fdr correction"},
]

try:
    from Config import *
except Exception as e:
    logger.warning("No config file found, using defaults")
    organisms = {
        'Mus musculus': 'NCBITaxon:10090',
        # 'Homo sapiens': 'NCBITaxon:9606',
        # 'Bacillus cereus': "NCBITaxon:405534",
        # 'Saccharomyces cerevisiae': 'NCBITaxon:4932',
        # 'Escherichia coli': 'NCBITaxon:562',
        # 'Drosophila melanogaster': 'NCBITaxon:7227',
        # 'Rattus norvegicus': 'NCBITaxon:10116',
        # 'Bos taurus': 'NCBITaxon:9913',
        # 'Caenorhabditis elegans': 'NCBITaxon:6239',
        # 'Pseudomonas aeruginosa': 'NCBITaxon:287',
        # 'Arabidopsis thaliana': 'NCBITaxon:3702',
    }


def shorten_label(label, max_len = 10):
    return label if len(label) <= max_len else label[:max_len].rsplit(" ", 1)[0] + "..."



class SunburstTerm:
    def __init__(self, _term, _entry_point, _pvalue = None):
        self.term = _term
        self.entry_point = _entry_point
        self.parent = ""
        self.pvalue = _pvalue


LIBRARY_TO_DOMAINS = {
    "GO": {"Biological process", "Cellular component", "Molecular function"},
    "SNP": {"Metabolic and signalling pathway"},
    "R-": {"Metabolic and signalling pathway"},
    "LION": {"Physical or chemical properties"},
    "DOID": {"Disease"},
    "MONDO": {"Disease"},
    "HP": {"Phenotype"},
}


face = freetype.Face(f"{current_path}/assets/DejaVuSans.ttf")
char_sizes = []
for font_size in range(1, 23):
    char_sizes_font = []
    char_sizes.append(char_sizes_font)
    face.set_char_size(14 * 64)
    for char in range(256):
        face.load_char(chr(char))
        char_sizes_font.append((face.glyph.advance.x >> 6) * 0.1 * font_size)


LINK_COLOR = "#2980B9"
SESSION_DURATION_TIME = 60 * 60 * 2 # two hours
domain_colors = {
    "Biological process": ["#FFB3BA", "#F6C0C5"],
    "Cellular component": ["#BAFFC9", "#C7F6D4"],
    "Disease": ["#BAE1FF", "#C7E2F6"],
    "Metabolic and signalling pathway": ["#FFD1BA", "#F6D8C7"],
    "Molecular function": ["#BAFFF1", "#C7F6EE"],
    "Phenotype": ["#D1BAFF", "#D8C7F6"],
    "Physical or chemical properties": ["#FFFFBA", "#F6F6C7"],
    "Enzymatic activity (Swiss-Prot)": ["#F6D4C7", "#f0c9f3"],
    "Enzymatic activity (Swiss-Prot + TrEMBL)": ["#BAFFD1", "#C7F6D8"],
}
REGULATED_COLOR = "#F49E4C"
ASSOCIATED_COLOR = "#87BCDE"
INNER_CIRCLE = 300
CIRCLE_WIDTH = 300
BAR_SORTING_PVALUE = "pvalue"
BAR_SORTING_SIMILARITY = "sim"

LIPIDS_COLOR = "#A3DC9A"
PROTEINS_COLOR = "#DEE791"
METABOLITES_COLOR = "#FFF9BD"
TRANSCRIPTS_COLOR = "#FFD6BA"


regulated_molecule_handling = [
    {"value": MOLECULE_HANDLING_ERROR, "label": "Halt analysis and report"},
    {"value": MOLECULE_HANDLING_REMOVE, "label": "Remove for analysis"},
]

molecule_handling = [
    {"value": MOLECULE_HANDLING_ERROR, "label": "Halt analysis and report"},
    {"value": MOLECULE_HANDLING_REMOVE, "label": "Remove for analysis"},
    {"value": MOLECULE_HANDLING_IGNORE, "label": "Keep and continue"},
]

term_representation = [
    {"value": "two-sided", "label": "Term over/under-represented"},
    {"value": "less", "label": "Term under-represented"},
    {"value": "greater", "label": "Term over-represented"},
]


def get_path(nodes, node):
    if node not in nodes: return []
    path = []
    current_node = node
    while current_node != None:
        path.append(current_node)
        current_node = nodes[current_node]
    return path[::-1]




class ItemCounter:
    def __init__(self, first_item = None):
        self.counter = {}
        if first_item != None: self.counter[first_item] = 1

    def add(self, item):
        if item not in self.counter: self.counter[item] = 1
        else: self.counter[item] += 1

    def __getitem__(self, item):
        if item in self.counter: return self.counter[item]
        return None

    def merge(self, foreign_counter):
        for item, count in foreign_counter.counter.items():
            if item in self.counter: self.counter[item] += count
            else: self.counter[item] = count


def get_term_link(term_id):
    if term_id.startswith("GO:"):
        return "https://www.ebi.ac.uk/QuickGO/term/" + term_id

    elif term_id.startswith("SMP"):
        return "https://pathbank.org/view/" + term_id

    elif term_id.startswith("LION:") or term_id.startswith("CAT:"):
        return "https://bioportal.bioontology.org/ontologies/LION?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F" + term_id.replace(":", "_")

    elif term_id.startswith("RHEA:"):
        return f"https://www.rhea-db.org/rhea/{term_id[5:]}"

    elif term_id.startswith("UNIPROT:"):
        return f"https://www.uniprot.org/uniprotkb/{term_id[8:]}/entry"

    elif term_id.startswith("CHEBI:"):
        return f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={term_id}"

    elif term_id.startswith("DOID:"):
        return f"https://disease-ontology.org/?id={term_id}"

    elif term_id.startswith("MONDO:"):
        return f"https://monarchinitiative.org/{term_id}"

    elif term_id.startswith("HGNC:"):
        return f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{term_id}"

    elif term_id.startswith("HP:"):
        return f"https://hpo.jax.org/browse/term/{term_id}"

    elif term_id.startswith("NCBI:"):
        return f"https://www.ncbi.nlm.nih.gov/gene/{term_id.split(':')[1]}"

    elif term_id.startswith("ENS") or term_id.startswith("WBGene") or term_id.startswith("FBgn"):
        return f"https://www.ensembl.org/id/{term_id}"

    elif term_id.startswith("R-"):
        return f"https://reactome.org/PathwayBrowser/#/{term_id}"

    elif term_id.startswith("LM"):
        return f"https://www.lipidmaps.org/databases/lmsd/{term_id}"

    return ""

def analytics(action):
    def send_action(recorded_action):
        try:
            url = f"https://lifs-tools.org/matomo/matomo.php?idsite=17&rec=1&e_c=MAPtoGO-{ini_config["git"]["version"]}&e_a={recorded_action}"
            logger.debug(f"Sending analytics data to Matomo without token: {url}")
            response = requests.get(url, timeout = 2)

        except Exception as e:
            logger.warning(e)
            pass

    logger.info(f"Recording action: {action}")
    thread1 = threading.Thread(target = send_action, args = (action,))
    thread1.start()



def annotate_arc(
    figure,
    annotations,
    text_to_add,
    mid_angle,
    angle_width,
    distance,
    font_size = 14,
    max_radius = 400,
    color = "black",
    max_distance = -1,
    shorten = False,
):
    text_space = int(-20 * distance / max_distance + 19) if max_distance > -1 else 0
    text_width = sum(char_sizes[font_size][ord(char)] for char in text_to_add) + text_space * len(text_to_add)
    if 360 / (2 * max_radius * np.pi) * text_width > angle_width:
        new_text_width = angle_width / 360 * (2 * max_radius * np.pi)
        if shorten:
            text_to_add = shorten_label(text_to_add, int(len(text_to_add) / text_width * new_text_width * 0.9))
            text_width = sum(char_sizes[font_size][ord(char)] for char in text_to_add) + text_space * len(text_to_add)
        else:
            text_width = new_text_width

    # Calculate positions of each character along the arc
    start_pos = (2 * max_radius * np.pi) / 360 * mid_angle + text_width / 2
    # Compute the positions (x, y) for each character
    x, y, rotation_angles, text = [], [], [], []
    for i, char in enumerate(text_to_add):
        char_size = char_sizes[font_size][ord(char)]
        current_pos = start_pos - char_size / 2
        start_pos -= char_size + text_space

        if char == " ": continue
        rotation_angle = np.mod(360 / (2 * max_radius * np.pi) * current_pos, 360)
        # Convert polar coordinates (r, theta) to Cartesian (x, y)
        x = distance * np.cos(np.radians(rotation_angle))
        y = distance * np.sin(np.radians(rotation_angle))

        annotations.append(
            dict(
                x = x, y = y,
                text = f"<b>{char}</b>",
                showarrow = False,
                textangle = -rotation_angle + 90,
                font = dict(color = color, size = font_size),
                xanchor = "center",
                yanchor = "middle"
            )
        )



def add_arc(figure, start_angle, end_angle, arc_inner_radius, arc_outer_radius, color, hoverinfo = None, degree_factor = 1):
    theta_deg = np.linspace(start_angle, end_angle, int(abs(end_angle - start_angle) * degree_factor))
    theta_rad = np.deg2rad(theta_deg)

    x_outer = arc_outer_radius * np.cos(theta_rad)
    y_outer = arc_outer_radius * np.sin(theta_rad)
    x_inner = arc_inner_radius * np.cos(theta_rad[::-1])
    y_inner = arc_inner_radius * np.sin(theta_rad[::-1])
    x_arc = np.concatenate([x_outer, x_inner])
    y_arc = np.concatenate([y_outer, y_inner])

    scatter_dict = dict(
        x = x_arc,
        y = y_arc,
        fill = 'toself',
        fillcolor = color,
        line = dict(color = 'rgba(0, 0, 0, 0)'),
        mode = 'lines',
        showlegend = False,
    )
    if hoverinfo is not None:
        scatter_dict["text"] = hoverinfo
        scatter_dict["hoverinfo"] = "text"
    else:
        scatter_dict["hoverinfo"] = "skip"
    figure.add_trace(go.Scatter(**scatter_dict))






CC4_LINK = html.A(
    "CC BY 4.0",
    href = "https://creativecommons.org/licenses/by/4.0/legalcode#s3a1",
    target = "_blank",
    style = {"color": LINK_COLOR},
)

CC0_LINK = html.A(
    "CC0",
    href = "https://creativecommons.org/public-domain/cc0/",
    target = "_blank",
    style = {"color": LINK_COLOR},
)



sessions, examples = {}, {}
logger.info("Loading examples table")
xl = pd.ExcelFile(f"{current_path}/Data/examples.xlsx")
for worksheet_name in xl.sheet_names:
    if worksheet_name[0] in {"_", "#"}: continue
    df = xl.parse(worksheet_name)
    worksheet = {
        "bgl": [], "regl": [], "upregl": [], "downregl": [],
        "bgp": [], "regp": [], "upregp": [], "downregp": [],
        "bgm": [], "regm": [], "upregm": [], "downregm": [],
        "bgt": [], "regt": [], "upregt": [], "downregt": [],
        "title": "", "desc": "", "comp": "", "doi": "", "org": INIT_ORGANISM
    }
    if "BackgroundLipids" in df: worksheet["bgl"] = list(df["BackgroundLipids"].dropna())
    if "RegulatedLipids" in df: worksheet["regl"] = list(df["RegulatedLipids"].dropna())
    if "UpregulatedLipids" in df: worksheet["upregl"] = list(df["UpregulatedLipids"].dropna())
    if "DownregulatedLipids" in df: worksheet["downregl"] = list(df["DownregulatedLipids"].dropna())
    if "BackgroundProteins" in df: worksheet["bgp"] = list(df["BackgroundProteins"].dropna())
    if "RegulatedProteins" in df: worksheet["regp"] = list(df["RegulatedProteins"].dropna())
    if "UpregulatedProteins" in df: worksheet["upregp"] = list(df["UpregulatedProteins"].dropna())
    if "DownregulatedProteins" in df: worksheet["downregp"] = list(df["DownregulatedProteins"].dropna())
    if "BackgroundMetabolites" in df: worksheet["bgm"] = list(df["BackgroundMetabolites"].dropna())
    if "RegulatedMetabolites" in df: worksheet["regm"] = list(df["RegulatedMetabolites"].dropna())
    if "UpregulatedMetabolites" in df: worksheet["upregm"] = list(df["UpregulatedMetabolites"].dropna())
    if "DownregulatedMetabolites" in df: worksheet["downregm"] = list(df["DownregulatedMetabolites"].dropna())
    if "BackgroundTranscripts" in df: worksheet["bgt"] = list(df["BackgroundTranscripts"].dropna())
    if "RegulatedTranscripts" in df: worksheet["regt"] = list(df["RegulatedTranscripts"].dropna())
    if "UpregulatedTranscripts" in df: worksheet["upregt"] = list(df["UpregulatedTranscripts"].dropna())
    if "DownregulatedTranscripts" in df: worksheet["downregt"] = list(df["DownregulatedTranscripts"].dropna())
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
    "scrollZoom": False,
    "displayModeBar": True,
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
app.server.secret_key = "ce7a4618ff121b96faea2c896421ba5e"
app.title = APPLICATION_SHORT_TITLE


@app.server.before_request
def set_udata():
    if "uid" not in uuid_session:
        uuid_session["uid"] = str(uuid.uuid4())

api = Api(
    app.server,
    version = "1.0",
    title = f"{APPLICATION_TITLE} - REST API",
    description = "Application programming interface for the GO multiomics analysis platform",
    doc = app.get_relative_path("/api/docs"),
    prefix = app.get_relative_path("")
)
api.static_url_path = app.get_relative_path("/swaggerui")

logger.info(f"Deploying swagger-ui at {app.get_relative_path('/api/docs')}, setting prefix to: {app.get_relative_path('')}, static url path to: {app.get_relative_path('/swaggerui')}")


gc.disable()
enrichment_ontologies = {}
for tax_name, tax_id in organisms.items():
    logger.info(f"Loading {tax_name}")
    tax_id_number = tax_id.replace("NCBITaxon:", "")
    enrichment_ontologies[tax_id] = EnrichmentOntology(f"{current_path}/Data/ontology_{tax_id_number}.gz", tax_name)
gc.enable()
gc.collect()


def get_aggrid_modal(name, molecule):
    return dag.AgGrid(
        id = name,
        columnDefs = [
            {
                'field': "molecule_id",
                "cellDataType": False,
                "hide": True,
            },
            {
                'field': "molecule",
                "headerName": molecule,
                "cellRenderer": "MoleculeRenderer",
            },
            {
                'field': "regulated",
                "headerName": "Regulated",
                "maxWidth": 150,
            },
        ],
        rowData = [],
        columnSize = "responsiveSizeToFit",
        defaultColDef = {
            "suppressMovable": True,
            "sortable": True,
            "filter": True,
        },
        dashGridOptions = {
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
        style = {"height": "500px"},
    )


def session_timer_trigger(sessions):
    while True:
        time.sleep(10)
        current_time = time.time()
        sessions_to_delete = []
        for session_id, session_data in sessions.items():
            if current_time - session_data.time > SESSION_DURATION_TIME:
                sessions_to_delete.append(session_id)
        for session_id in sessions_to_delete:
            del sessions[session_id]
            logger.info(f"Deleting session: {session_id}")

thread = threading.Thread(target = session_timer_trigger, args = (sessions,), daemon = True)
thread.start()
logger.info("loaded")



def layout():
    ontology = enrichment_ontologies[INIT_ORGANISM]
    predefined_proteins = [
        len(ontology.reviewed_proteins),
        len(set(ontology.proteins.keys()) - ontology.reviewed_proteins),
        len(ontology.proteins),
    ]
    link_margin = "14px"

    return html.Div([
        dcc.Download(id = "download_data"),
        dmc.TextInput(id = "html_download_trigger", style = {"display": "none"}),
        html.Div("", id = "sessionid", style = {"display": "none"}),
        html.Div("", id = "analytics", style = {"display": "none"}),
        dcc.Location(id="url"),
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
        #dmc.SimpleGrid([
        html.Div([
            dmc.Group([
                dmc.Image(
                    src = dash.get_app().get_asset_url(f"golipids.png"),
                    style = {"width": "48px"},
                ),
                html.Div(
                    dmc.Title(APPLICATION_TITLE),
                ),
            ]),
            html.Div([
                html.A(
                    dmc.Text(
                        "LIFS Tools",
                    ),
                    id = "lifs_tools_link",
                    href = "https://lifs-tools.org/",
                    target = "_blank",
                    style = {"textDecoration": "none", "color": LINK_COLOR, "cursor": "pointer"},
                ),
                dmc.Text("|", style = {"marginLeft": link_margin, "marginRight": link_margin}),
                html.A(
                    dmc.Text(
                        "Load example datasets",
                    ),
                    id = "load_examples_link",
                    style = {"color": LINK_COLOR, "cursor": "pointer"},
                ),
                dmc.Text("|", style = {"marginLeft": link_margin, "marginRight": link_margin}),
                html.A(
                    dmc.Text(
                        "Description & disclaimer",
                    ),
                    id = "disclaimer_link",
                    style = {"color": LINK_COLOR, "cursor": "pointer"},
                ),
                dmc.Text("|", style = {"marginLeft": link_margin, "marginRight": link_margin}),
                html.A(
                    dmc.Text(
                        "REST API",
                    ),
                    id = "api_link",
                    href = dash.get_app().get_relative_path("/api/docs"),
                    target = "_blank",
                    style = {"textDecoration": "none", "color": LINK_COLOR, "cursor": "pointer"},
                ),
                dmc.Text(" ", style = {"marginLeft": "5px"}),
                ], style = {"height": "100%", "display": "flex", "alignItems": "flex-start", "justifyContent": "right"},
            ),
        ],  style={
                "display": "flex",
                "justifyContent": "space-between",
                "width": "100%",
            }
        ),
        #], cols = 2),
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
                html.Div(
                    "",
                    id = "term_molecules_modal_id",
                    style = {"display": "none"},
                ),
                dmc.Grid(
                    children = [
                        dmc.Col(
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
                                            dmc.Tab(
                                                "Transcripts",
                                                value = "transcript_tab_modal",
                                                id = "transcript_tab_modal_tab",
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
                                    dmc.TabsPanel(
                                        get_aggrid_modal("term_transcripts_modal_grid", "Transcript"),
                                        value = "transcript_tab_modal",
                                    ),
                                ],
                                id = "term_molecules_modal_tab",
                                color = "blue", # default is blue
                                orientation = "horizontal", # or "vertical"
                                variant = "outline", # or "outline" or "pills"
                                value = "lipid_tab_modal",
                            ),
                            span = 6,
                        ),
                        dmc.Col(
                            html.Div(
                                children = [
                                    dmc.Center(
                                        dmc.Text("Term path"),
                                    ),
                                    dmc.ScrollArea(
                                        html.Div(
                                            id = "term_path_area",
                                        ),
                                        h = 414,
                                    )
                                ]
                            ),
                            span = 6,
                        ),
                    ],
                    gutter = "xl",
                ),
                html.Div(
                    id = "contingency_table",
                    style = {"marginTop": "10px"},
                ),
            ],
            size = "60%",
        ),
        dmc.Modal(
            id = "barplot_terms_modal",
            zIndex = 10000,
            children = [
                html.Div(
                    dcc.Graph(
                        id = "barplot_terms",
                        config = {**plotly_config, **{"responsive": True}},
                        style = {
                            "height": "100%",
                            "width": "100%",
                        },
                    ),
                    id = "barplot_terms_wrapper",
                    style = {
                        "resize": "both",
                        "overflow": "hidden",
                        "width": "100%",
                        "height": "70vh",
                    }
                ),
                html.Div(
                    dmc.SimpleGrid(
                        [
                            dmc.NumberInput(
                                id = "barplot_numberinput_font_size",
                                label = "Font size",
                                value = 14,
                                min = 1,
                                max = 20,
                            ),
                            dmc.NumberInput(
                                id = "barplot_numberinput_connect_ths",
                                label = "Connectivity Threshold",
                                value = 80,
                                min = 0,
                                max = 100,
                            ),
                            dmc.Select(
                                id = "barplot_select_name",
                                label = "Select bar label:",
                                data = [{"value": "id", "label": "Term ID"}, {"value": "name", "label": "Term name"}],
                                value = "id",
                            ),
                            dmc.Select(
                                id = "barplot_select_sorting",
                                label = "Select sorting:",
                                data = [{"value": BAR_SORTING_PVALUE, "label": "p-value"}, {"value": BAR_SORTING_SIMILARITY, "label": "Molecule similarity in term"}],
                                value = BAR_SORTING_PVALUE,
                            ),
                        ],
                        cols = 2,
                    ),
                    id = "barplot_controls",
                    style = {"marginTop": "10px"},
                ),
                html.Div(
                    dmc.SimpleGrid(
                        [
                            dmc.Text("Maximum p-value", size = "xs", fw = 500),
                            dmc.Text("Minimum p-value", size = "xs", fw = 500),
                            dcc.Input(
                                id = "barplot_numberinput_max_pvalue",
                                value = "0.05",
                                debounce = True,
                                style = {
                                    "borderRadius": "4px",
                                    "paddingLeft": "12px",
                                    "paddingRight": "12px",
                                    "fontSize": "14px",
                                    "border": "1px solid rgb(206, 212, 218)",
                                    "height": "32px",
                                    "fontFamily": "-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, Helvetica, Arial, sans-serif, Apple Color Emoji, Segoe UI Emoji",
                                },
                            ),
                            dcc.Input(
                                id = "barplot_numberinput_min_pvalue",
                                value = "0.0001",
                                debounce = True,
                                style = {
                                    "borderRadius": "4px",
                                    "paddingLeft": "12px",
                                    "paddingRight": "12px",
                                    "fontSize": "14px",
                                    "border": "1px solid rgb(206, 212, 218)",
                                    "height": "32px",
                                    "fontFamily": "-apple-system, BlinkMacSystemFont, Segoe UI, Roboto, Helvetica, Arial, sans-serif, Apple Color Emoji, Segoe UI Emoji",
                                },
                            ),
                        ],
                        cols = 2,
                        verticalSpacing = "2px",
                        spacing = "xs",
                    ),
                    id = "sunburst_controls",
                    style = {"marginTop": "10px"},
                ),
            ],
            size = "90%",
        ),
        dmc.Modal(
            id = "sankey_modal",
            zIndex = 10000,
            children = [
                html.Div(
                    style={
                        "display": "flex",
                        "flexDirection": "row",
                        #"height": "100vh",   # full height
                    },
                    children=[
                        html.Div(
                            dcc.Graph(
                                id = "sankey_graph",
                                config = {**plotly_config, **{"responsive": True}},
                                style = {
                                    "height": "100%",
                                    "width": "100%",
                                },
                            ),
                            id = "sankey_graph_wrapper",
                            style = {
                                "resize": "both",
                                "overflow": "hidden",
                                "width": "100%",
                                "height": "70vh",
                                "flex": "1",
                            }
                        ),
                        html.Div(
                            style={
                                "flex": "1",
                                "padding": "10px",
                                "minWidth": "250px",
                            },
                            children=[
                                html.Div(
                                    style={
                                        "display": "flex",
                                        "justifyContent": "space-between",
                                        "alignItems": "center",
                                        "marginBottom": "0px",
                                    },
                                    children=[
                                        dmc.Title(
                                            "Node entries",
                                            order = 5,
                                            style = {"margin": "0px"},
                                        ),
                                        dmc.ActionIcon(
                                            DashIconify(icon = "material-symbols:download-rounded", width = 20),
                                            id = "icon_download_sankey_entries",
                                            title = "Download table",
                                            disabled = True,
                                        ),
                                    ],
                                ),

                                dag.AgGrid(
                                    id = "sankey_entries",
                                    columnDefs= [
                                        {
                                            "headerName": "Entry",
                                            "field": "termid",
                                            "cellRenderer": "TermIDRenderer",
                                        },
                                        {
                                            "headerName": "Count",
                                            "field": "count",
                                            "filter": "agNumberColumnFilter",
                                            "maxWidth": 100,
                                        },
                                    ],
                                    rowData = [],
                                    columnSize = "responsiveSizeToFit",
                                    defaultColDef={"sortable": True, "filter": True, "resizable": True},
                                    style={"width": "100%", "height": "100%"},
                                )
                            ],
                        ),
                    ],
                ),
                dmc.RadioGroup(
                    children = [
                        dmc.Group([
                            dmc.Radio("All molecules", value = "sankey_all"),
                            dmc.Radio("All non-regulated molecules", value = "sankey_non"),
                            dmc.Radio("All regulated molecules", value = "sankey_regulated"),
                            dmc.Radio(
                                "All up-regulated molecules",
                                value = "sankey_upregulated",
                                id = "sankey_upregulated",
                                style = {"display": "none"},
                            ),
                            dmc.Radio(
                                "All down-regulated molecules",
                                value = "sankey_downregulated",
                                id = "sankey_downregulated",
                                style = {"display": "none"},
                            ),
                        ], my = 10)
                    ],
                    id = "radiogroup_sankey",
                    value = "sankey_all",
                    size = "sm",
                    mb = 10,
                    style = {"paddingBottom": "8px", "paddingRight": "2px", "display": "block"},
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
        children = dmc.ScrollArea(html.Div([
            html.P([
                dmc.Title(f"About {APPLICATION_SHORT_TITLE}", order = 4),
                dmc.Text(
                    f"{APPLICATION_SHORT_TITLE} is the first implementation that enables an integrated GO term enrichment analysis on lipids, proteins, metabolites, and transcripts simultaneously. Starting with a list of identified lipids (at least on the species level) as background and a list with only the (differentially) regulated lipids, each lipid is mapped to proteins involved in metabolic reactions. On this base, GO term analysis will then be performed with a statistical test. Depending on the selected domain (biological process, molecular function, cellular compartment, physical or chemical properties, or metabolic and signaling pathways), the analysis provides a sorted list of GO terms with a p-value for each term.",
                    style = {"textAlign": "justify"},
                ),
            ]),

            html.P([
                dmc.Text("Responsible people for this database are:"),
                dmc.Text("Dominik Kopczynski: Implementation of both front- and backend"),
                dmc.Text("Cristina Coman: Method validation"),
                dmc.Text("Nils Hoffmann: Application coordination and testing"),
                dmc.Text("Robert Ahrends: Project leader"),
            ]),

            html.P([
                dmc.Title("Version", order = 4),
                dmc.Text(f"Current version: {get_latest_tag()} ({get_git_info()})"),
            ]),

            html.P([
                dmc.Title("Disclaimer", order = 4),
                dmc.Title("Legal form and governing bodies", order = 5),
                dmc.Text(
                    "The University of Vienna is a legal entity under public law in accordance with section 4 of the 2002 Universities Act. In accordance with section 20 of the 2002 Universities Act, the senior governing bodies of the University of Vienna are the University Board, the Rectorate, the Rector and the Senate.",
                    style = {"textAlign": "justify"},
                ),
            ]),

            html.P([
                dmc.Title("Basic purpose of the website", order = 5),
                dmc.Text(
                    "Providing a web service for multiomics enrichment analyses.",
                    style = {"textAlign": "justify"},
                ),
            ]),

            html.P([
                dmc.Title("Copyright", order = 5),
                dmc.Text(
                    "© The contents of the website are published in the World Wide Web for online access. Copyright for texts, images, graphic design and source code is owned by the Department of Analytical Chemistry, University of Vienna and subject to legal protection. The production, use and non-commercial distribution of copies in electronic or printed form is allowed provided that the contents are not altered and the source is mentioned (Source: Department of Analytical Chemistry, University of Vienna).",
                    style = {"textAlign": "justify"},
                ),
            ]),

            html.P([
                dmc.Title("Liability", order = 5),
                dmc.Text(
                    "The text provided on this homepage has been reviewed with great care. However, the Department of Analytical Chemistry or the University of Vienna cannot guarantee the accuracy, completeness or validity of the information provided. Therefore the Department of Analytical Chemistry and the University of Vienna accept no liability for the contents provided. Links to other websites have been carefully selected. However, the Department of Analytical Chemistry and the University of Vienna are not responsible for contents on any other websites.",
                    style = {"textAlign": "justify"},
                ),
            ]),

            html.P([
                dmc.Title("Responsibility for contents and editing", order = 5),
                dmc.Text("Department of Analytical Chemistry"),
                dmc.Text("University of Vienna"),
                dmc.Text("Sensengasse 8 / TOP 12"),
                dmc.Text("A-1090 Vienna"),
                dmc.Text("Austria"),
                html.Br(),
                dmc.Text("Tel: +43-1-4277-52307"),
                dmc.Text([
                    "Email: ",
                    html.A(
                        "dominik.kopczynski@univie.ac.at",
                        href = "mailto:dominik.kopczynski@univie.ac.at",
                        style = {"color": LINK_COLOR},
                    ),
                ]),
            ]),

            html.P([
                dmc.Title("Third-party databeses", order = 5),
                dmc.Text("GO multiomics is using data from third-party databases. All databases are listed with their according licenses or permission:"),
                dmc.Text([
                    "- ",
                    html.A(
                        "ChEBI",
                        href = "https://www.ebi.ac.uk/chebi/aboutChebiForward.do",
                        target = "_blank",
                        style = {"color": LINK_COLOR},
                    ),
                    ": ",
                    CC4_LINK
                ]),
                dmc.Text([
                    "- ",
                    html.A(
                        "Disease Ontology",
                        href = "https://github.com/DiseaseOntology/HumanDiseaseOntology/tree/main",
                        target = "_blank",
                        style = {"color": LINK_COLOR},
                    ),
                    ": ",
                    CC0_LINK
                ]),
                dmc.Text([
                    "- ",
                    html.A(
                        "Ensembl",
                        href = "https://www.ensembl.org/info/about/legal/disclaimer.html",
                        target = "_blank",
                        style = {"color": LINK_COLOR},
                    ),
                    ": no restriction (details on clicking link)"
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
                    CC4_LINK
                ]),
                dmc.Text([
                    "- ",
                    html.A(
                        "HUGO Gene Nomenclature Committee (HGNC)",
                        href = "https://www.genenames.org/about/license/",
                        target = "_blank",
                        style = {"color": LINK_COLOR},
                    ),
                    ": ",
                    CC0_LINK
                ]),
                dmc.Text([
                    "- ",
                    html.A(
                        "Human Phenotype Ontology",
                        href = "https://hpo.jax.org/",
                        target = "_blank",
                        style = {"color": LINK_COLOR},
                    ),
                    ": ",
                    html.A(
                        "License",
                        href = "https://hpo.jax.org/license",
                        target = "_blank",
                        style = {"color": LINK_COLOR},
                    ),
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
                    CC4_LINK
                ]),
                dmc.Text([
                    "- ",
                    html.A(
                        "Mondo Disease Ontology",
                        href = "https://mondo.monarchinitiative.org/pages/download/",
                        target = "_blank",
                        style = {"color": LINK_COLOR},
                    ),
                    ": ",
                    CC4_LINK
                ]),
                # dmc.Text([
                #     "- ",
                #     html.A(
                #         "Pathbank",
                #         href = "https://pathbank.org/about",
                #         target = "_blank",
                #         style = {"color": LINK_COLOR},
                #     ),
                #     ": ",
                #     html.A(
                #         "Open database license",
                #         href = "https://opendatacommons.org/licenses/odbl/1-0/",
                #         target = "_blank",
                #         style = {"color": LINK_COLOR},
                #     ),
                # ]),
                dmc.Text([
                    "- ",
                    html.A(
                        "Reactome",
                        href = "https://reactome.org/license",
                        target = "_blank",
                        style = {"color": LINK_COLOR},
                    ),
                    ": ",
                    CC0_LINK
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
                    CC4_LINK
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
                    CC4_LINK
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
                    CC4_LINK
                ]),
            ]),
        ], style = {"width": "98%"}), h = 450, offsetScrollbars = True, scrollHideDelay = 0),
        size = "50%",
    ),
    dmc.SimpleGrid(
        cols = 2,
        style = {"marginTop": "5px", "gridTemplateColumns": "44% 55%"},
        children = [
            html.Div([
                dmc.Switch(
                    "Separate between up / down regulated molecules",
                    id = "separate_updown_switch",
                    style = {"marginBottom": "5px"},
                ),
                dmc.Tabs(
                    [
                        dmc.TabsList(
                            [
                                dmc.Tab("Lipids", value="lipid_tab"),
                                dmc.Tab("Proteins", value="protein_tab"),
                                dmc.Tab("Metabolites", value="metabolites_tab"),
                                dmc.Tab("Transcripts", value="transcripts_tab"),
                            ],
                            grow = True,
                        ),
                        dmc.TabsPanel([
                            dmc.SimpleGrid([
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All lipid names in experiment (background)",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_all_lipids",
                                            value = "",
                                            style = {"width": "100%", "height": "100%"},
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_all_lipids",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {"height": TEXTFIELD_HEIGHT, "display": "flex", "flexDirection": "column"},
                                    className = "textarea-expand-container",
                                ),
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All regulated lipid names in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_regulated_lipids",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_regulated_lipids",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {
                                        "height": TEXTFIELD_HEIGHT,
                                        "display": "flex",
                                        "flexDirection": "column",
                                    },
                                    className = "textarea-expand-container",
                                    id = "field_regulated_lipids",
                                ),
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All up-regulated lipid names in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_upregulated_lipids",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_upregulated_lipids",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),

                                        dmc.Title(
                                            "All down-regulated lipid names in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_downregulated_lipids",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_downregulated_lipids",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {
                                        "height": TEXTFIELD_HEIGHT,
                                        "display": "none",
                                        "flexDirection": "column",
                                    },
                                    className = "textarea-expand-container",
                                    id = "field_up_down_regulated_lipids",
                                ),
                            ], cols = 2)],
                            value="lipid_tab",
                        ),

                        dmc.TabsPanel([
                            dmc.SimpleGrid([
                                html.Div(
                                    [
                                        dmc.Group([
                                            dmc.Title(
                                                "All protein accessions in experiment (background)",
                                                order = 5,
                                                style = {"marginTop": "10px"},
                                            ),
                                            dmc.Menu(
                                                [
                                                    dmc.MenuTarget(
                                                        dmc.ActionIcon(
                                                            DashIconify(icon = "simple-icons:helix", width = 16),
                                                            id = "predefined_bg_proteome",
                                                            title = "Select predefined background proteome",
                                                            style = {"display": "flex", "alignItems": "flex-end"},
                                                        ),
                                                    ),
                                                    dmc.MenuDropdown(
                                                        [
                                                            dmc.MenuItem(
                                                                f"Reviewed proteins only ({predefined_proteins[0]})",
                                                                id = "button_load_reviewed_proteins",
                                                                n_clicks = 0,
                                                            ),
                                                            dmc.MenuItem(
                                                                f"Unreviewed proteins only ({predefined_proteins[1]})",
                                                                id = "button_load_unreviewed_proteins",
                                                                n_clicks = 0,
                                                            ),
                                                            dmc.MenuItem(
                                                                f"All registered proteins ({predefined_proteins[2]})",
                                                                id = "button_load_all_proteins",
                                                                n_clicks = 0,
                                                            ),
                                                        ]
                                                    ),
                                                ]
                                            ),
                                        ]),
                                        dmc.Textarea(
                                            id = "textarea_all_proteins",
                                            style = {"width": "100%", "height": "100%"},
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_all_proteins",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {"height": TEXTFIELD_HEIGHT, "display": "flex", "flexDirection": "column"},
                                    className = "textarea-expand-container",
                                ),
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All regulated protein accession in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_regulated_proteins",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_regulated_proteins",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {
                                        "height": TEXTFIELD_HEIGHT,
                                        "display": "flex",
                                        "flexDirection": "column",
                                    },
                                    className = "textarea-expand-container",
                                    id = "field_regulated_proteins",
                                ),
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All up-regulated protein accession in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_upregulated_proteins",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_upregulated_proteins",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                        dmc.Title(
                                            "All down-regulated protein accession in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_downregulated_proteins",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_downregulated_proteins",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {
                                        "height": TEXTFIELD_HEIGHT,
                                        "display": "none",
                                        "flexDirection": "column",
                                    },
                                    className = "textarea-expand-container",
                                    id = "field_up_down_regulated_proteins",
                                ),
                            ], cols = 2)],
                            value="protein_tab",
                        ),

                        dmc.TabsPanel([
                            dmc.SimpleGrid([
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All metabolite ChEBI Ids in experiment (background)",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_all_metabolites",
                                            style = {"width": "100%", "height": "100%"},
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_all_metabolites",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),

                                    ],
                                    style = {"height": TEXTFIELD_HEIGHT, "display": "flex", "flexDirection": "column"},
                                    className = "textarea-expand-container",
                                ),
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All regulated metabolite ChEBI Ids in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_regulated_metabolites",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_regulated_metabolites",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {
                                        "height": TEXTFIELD_HEIGHT,
                                        "display": "flex",
                                        "flexDirection": "column",
                                    },
                                    className = "textarea-expand-container",
                                    id = "field_regulated_metabolites",
                                ),
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All up-regulated metabolite ChEBI Ids in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_upregulated_metabolites",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_upregulated_metabolites",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                        dmc.Title(
                                            "All down-regulated metabolite ChEBI Ids in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_downregulated_metabolites",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_downregulated_metabolites",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {
                                        "height": TEXTFIELD_HEIGHT,
                                        "display": "none",
                                        "flexDirection": "column",
                                    },
                                    className = "textarea-expand-container",
                                    id = "field_up_down_regulated_metabolites",
                                ),

                            ], cols = 2)],
                            value = "metabolites_tab",
                        ),

                        dmc.TabsPanel([
                            dmc.SimpleGrid([
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All ensembl Ids in experiment (background)",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_all_transcripts",
                                            style = {"width": "100%", "height": "100%"},
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_all_transcripts",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {"height": TEXTFIELD_HEIGHT, "display": "flex", "flexDirection": "column"},
                                    className = "textarea-expand-container",
                                ),
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All regulated ensembl Ids in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_regulated_transcripts",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_regulated_transcripts",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {
                                        "height": TEXTFIELD_HEIGHT,
                                        "display": "flex",
                                        "flexDirection": "column",
                                    },
                                    className = "textarea-expand-container",
                                    id = "field_regulated_transcripts",
                                ),
                                html.Div(
                                    [
                                        dmc.Title(
                                            "All up-regulated ensembl Ids in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_upregulated_transcripts",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_upregulated_transcripts",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                        dmc.Title(
                                            "All down-regulated ensembl Ids in experiment",
                                            order = 5,
                                            style = {"marginTop": "10px"},
                                        ),
                                        dmc.Textarea(
                                            id = "textarea_downregulated_transcripts",
                                            value = "",
                                        ),
                                        dmc.Group(
                                            dmc.Text(
                                                "Entries: 0",
                                                id = "num_downregulated_transcripts",
                                                style = {"color": "#808080"},
                                                size = "12px",
                                            ),
                                            position = "right",
                                        ),
                                    ],
                                    style = {
                                        "height": TEXTFIELD_HEIGHT,
                                        "display": "none",
                                        "flexDirection": "column",
                                    },
                                    className = "textarea-expand-container",
                                    id = "field_up_down_regulated_transcripts",
                                ),
                            ], cols = 2)],
                            value = "transcripts_tab",
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
                    style = {"marginTop": "0px"},
                ),
                html.Div(
                    [
                        dmc.SimpleGrid(
                            cols = 2,
                            children = [
                                dmc.Select(
                                    id = "select_organism",
                                    data = [
                                        {"value": organisms[key], "label": key} for key in sorted(organisms.keys())
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
                                    data = correction_list,
                                    value = "fdr_bh",
                                ),
                                html.Div([
                                    dmc.Group([
                                        dmc.Text(
                                            "Term representation:",
                                            fw = 500,
                                            size = 14,
                                            color = "rgb(33, 37, 41)",
                                        ),
                                        dmc.HoverCard(
                                            withArrow = True,
                                            width = 400,
                                            zIndex = 1000,
                                            shadow = "md",
                                            children = [
                                                dmc.HoverCardTarget(
                                                    DashIconify(
                                                        icon="material-symbols:live-help",
                                                        style={
                                                            "position": "relative",
                                                            "height": "100%",
                                                            "top": "-8px",
                                                            "left": "-10px",
                                                        },
                                                    )
                                                ),
                                                dmc.HoverCardDropdown(
                                                    dmc.Text("Over-representation: more biomolecules are regulated within a term, suggesting term is enriched. Under-representation: less biomolecules are regulated within a term, suggesting term is stable within experimental study."),
                                                )
                                            ],
                                        ),
                                    ]),
                                    dmc.Select(
                                        id = "select_term_representation",
                                        data = term_representation,
                                        value = "greater",
                                    ),
                                ]),
                                html.Div(
                                    dmc.Select(
                                        id = "select_molecule_handling",
                                        data = molecule_handling,
                                        value = MOLECULE_HANDLING_REMOVE,
                                        label = "Handling of unrecognizable molecules:",
                                    ),
                                ),
                                html.Div(
                                    dmc.Select(
                                        id = "select_regulated_molecule_handling",
                                        data = regulated_molecule_handling,
                                        value = MOLECULE_HANDLING_REMOVE,
                                        label = "Handling of non-background regulated molecules:",
                                    ),
                                ),
                                html.Div([
                                        dmc.Switch(
                                            id = "use_bounded_fatty_acyls",
                                            checked = False,
                                            label = "Use bounded fatty acyls for analysis, too",
                                            style = {"paddingBottom": "8px", "paddingRight": "2px"},
                                        ),
                                        dmc.HoverCard(
                                            withArrow = True,
                                            width = 400,
                                            zIndex = 1000,
                                            shadow = "md",
                                            children = [
                                                dmc.HoverCardTarget(
                                                    DashIconify(
                                                        icon="material-symbols:live-help",
                                                        style={
                                                            "position": "relative",
                                                            "height": "100%",
                                                            "top": "-8px",
                                                        },
                                                    )
                                                ),
                                                dmc.HoverCardDropdown(
                                                    dmc.Text("When activated, bounded fatty chains in lipids such as the 20:4 (AA) in 'PI 18:0/20:4' will be considered for the analysis, too."),
                                                )
                                            ],
                                        ),
                                    ],
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
                                dmc.Checkbox(
                                    label = "Use transcripts",
                                    id = "checkbox_use_transcripts",
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
                html.Div([
                    dmc.Title(
                        "Results",
                        order = 5,
                        style = {"marginTop": "10px"},
                    ),
                    html.Div(
                        html.Span(
                            dmc.Group([
                                dmc.MultiSelect(
                                    id = "multiselect_filter_molecules",
                                    placeholder = "Filter terms for molecules",
                                    style = {"marginRight": "15px", "width": "400px", "maxWidth": 400},
                                    searchable = True,
                                    className = "truncate-multiselect",
                                ),
                                dmc.ActionIcon(
                                    DashIconify(icon="ion:bar-chart-sharp", width = 20),
                                    id = "histogram_results",
                                    title = "Show p-value histogram",
                                    disabled = True,
                                    style = {"marginRight": "5px"},
                                ),
                                dmc.ActionIcon(
                                    DashIconify(icon = "tdesign:chart-column-filled", width = 20),
                                    id = "chart_results",
                                    title = "Show p-value chart",
                                    disabled = True,
                                    style = {"marginRight": "5px"},
                                ),
                                dmc.ActionIcon(
                                    DashIconify(icon = "carbon:chart-sunburst", width = 20),
                                    id = "sunburst_results",
                                    title = "Show term hierarchy",
                                    disabled = True,
                                    style = {"marginRight": "5px"},
                                ),
                                dmc.ActionIcon(
                                    DashIconify(icon = "carbon:sankey-diagram-alt", width = 20),
                                    id = "sankey_results",
                                    title = "Show term flow",
                                    disabled = True,
                                    style = {"marginRight": "5px"},
                                ),
                                dmc.ActionIcon(
                                    DashIconify(icon = "material-symbols:download-rounded", width = 20),
                                    id = "icon_download_results",
                                    title = "Download table",
                                    disabled = True,
                                ),
                            ], spacing = 0),
                        ),
                        style = {"height": "100%", "display": "flex", "alignItems": "flex-end", "justifyContent": "right"},
                    ),
                ], style = {
                    "display": "flex",
                    "justifyContent": "space-between",
                    #"width": "100%",  # or set fixed width like "400px"
                    #"padding": "0 10px",
                }),
                html.Div([
                    dag.AgGrid(
                        id = "graph_enrichment_results",
                        columnDefs = graph_enrichment_results_def(False, False),
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
                        className = "ag-theme-alpine terms-theme",
                        style = {
                            "height": "95%",
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
                    dmc.Text(
                        "Entries: 0",
                        id = "num_enrichment_terms",
                        style = {"color": "#808080"},
                        size = "12px",
                    ),
                ], className = "terms-overlay",
                ),
            ]),
        ],
    ),
    html.Div(style = {"height": "500px"}),
], style={
    "height": "100vh",
    #"background": "linear-gradient(to bottom, white 0%, white 0%, #e4eff7 8%, #e4eff7 100%)",
    "margin": "0",
    "padding": "0"
})

app.layout = layout



@callback(
    Output("select_domains", "value", allow_duplicate = True),
    Output("select_domains", "data", allow_duplicate = True),
    Output("button_load_reviewed_proteins", "children", allow_duplicate = True),
    Output("button_load_unreviewed_proteins", "children", allow_duplicate = True),
    Output("button_load_all_proteins", "children", allow_duplicate = True),
    Input("select_organism", "value"),
    State("select_domains", "value"),
    prevent_initial_call = True,
)
def organism_changed(organism, domain_values):
    if organism not in enrichment_ontologies:
        raise exceptions.PreventUpdate

    ontology = enrichment_ontologies[organism]
    domains = ontology.domains
    data = sorted(list(domains))
    values = [v for v in domain_values if v in domains]

    predefined_proteins = [
        len(ontology.reviewed_proteins),
        len(set(ontology.proteins.keys()) - ontology.reviewed_proteins),
        len(ontology.proteins),
    ]

    return (
        values,
        data,
        f"Reviewed proteins only ({predefined_proteins[0]})",
        f"Unreviewed proteins only ({predefined_proteins[1]})",
        f"All registered proteins ({predefined_proteins[2]})",
    )



@app.callback(
    Output("sessionid", "children", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("separate_updown_switch", "checked", allow_duplicate = True),
    Output("num_all_lipids", "children", allow_duplicate = True),
    Output("num_regulated_lipids", "children", allow_duplicate = True),
    Output("num_upregulated_lipids", "children", allow_duplicate = True),
    Output("num_downregulated_lipids", "children", allow_duplicate = True),
    Output("num_all_proteins", "children", allow_duplicate = True),
    Output("num_regulated_proteins", "children", allow_duplicate = True),
    Output("num_upregulated_proteins", "children", allow_duplicate = True),
    Output("num_downregulated_proteins", "children", allow_duplicate = True),
    Output("num_all_metabolites", "children", allow_duplicate = True),
    Output("num_regulated_metabolites", "children", allow_duplicate = True),
    Output("num_upregulated_metabolites", "children", allow_duplicate = True),
    Output("num_downregulated_metabolites", "children", allow_duplicate = True),
    Output("num_all_transcripts", "children", allow_duplicate = True),
    Output("num_regulated_transcripts", "children", allow_duplicate = True),
    Output("num_upregulated_transcripts", "children", allow_duplicate = True),
    Output("num_downregulated_transcripts", "children", allow_duplicate = True),
    Output("textarea_all_lipids", "value", allow_duplicate = True),
    Output("textarea_regulated_lipids", "value", allow_duplicate = True),
    Output("textarea_upregulated_lipids", "value", allow_duplicate = True),
    Output("textarea_downregulated_lipids", "value", allow_duplicate = True),
    Output("textarea_all_proteins", "value", allow_duplicate = True),
    Output("textarea_regulated_proteins", "value", allow_duplicate = True),
    Output("textarea_upregulated_proteins", "value", allow_duplicate = True),
    Output("textarea_downregulated_proteins", "value", allow_duplicate = True),
    Output("textarea_all_metabolites", "value", allow_duplicate = True),
    Output("textarea_regulated_metabolites", "value", allow_duplicate = True),
    Output("textarea_upregulated_metabolites", "value", allow_duplicate = True),
    Output("textarea_downregulated_metabolites", "value", allow_duplicate = True),
    Output("textarea_all_transcripts", "value", allow_duplicate = True),
    Output("textarea_regulated_transcripts", "value", allow_duplicate = True),
    Output("textarea_upregulated_transcripts", "value", allow_duplicate = True),
    Output("textarea_downregulated_transcripts", "value", allow_duplicate = True),
    Output("use_bounded_fatty_acyls", "checked", allow_duplicate = True),
    Output("checkbox_use_lipids", "checked", allow_duplicate = True),
    Output("checkbox_use_proteins", "checked", allow_duplicate = True),
    Output("checkbox_use_metabolites", "checked", allow_duplicate = True),
    Output("checkbox_use_transcripts", "checked", allow_duplicate = True),
    Output("select_organism", "value", allow_duplicate = True),
    Output("select_molecule_handling", "value", allow_duplicate = True),
    Output("select_regulated_molecule_handling", "value", allow_duplicate = True),
    Output("multiselect_filter_molecules", "value", allow_duplicate = True),
    Output("select_term_representation", "value", allow_duplicate = True),
    Output("select_test_method", "value", allow_duplicate = True),
    Output("select_domains", "value", allow_duplicate = True),
    Input("url", "pathname"),
    prevent_initial_call = True,
)
def load_uid(_):
    session_id = uuid_session.get("uid")
    if session_id not in sessions:
        sessions[session_id] = SessionEntry()
        logger.info(f"New session: {session_id}")
        return (session_id, *([no_update] * 47))

    logger.info(f"Reopen session: {session_id}")
    session = sessions[session_id]
    ui = session.ui
    use_bounded_fatty_acyls = session.use_bounded_fatty_acyls
    all_lipids_list = ui["all_lipids_list"] if "all_lipids_list" in ui else ""
    regulated_lipids_list = ui["regulated_lipids_list"] if "regulated_lipids_list" in ui else ""
    upregulated_lipids_list = ui["upregulated_lipids_list"] if "upregulated_lipids_list" in ui else ""
    downregulated_lipids_list = ui["downregulated_lipids_list"] if "downregulated_lipids_list" in ui else ""
    all_proteins_list = ui["all_proteins_list"] if "all_proteins_list" in ui else ""
    regulated_proteins_list = ui["regulated_proteins_list"] if "regulated_proteins_list" in ui else ""
    upregulated_proteins_list = ui["upregulated_proteins_list"] if "upregulated_proteins_list" in ui else ""
    downregulated_proteins_list = ui["downregulated_proteins_list"] if "downregulated_proteins_list" in ui else ""
    all_metabolites_list = ui["all_metabolites_list"] if "all_metabolites_list" in ui else ""
    regulated_metabolites_list = ui["regulated_metabolites_list"] if "regulated_metabolites_list" in ui else ""
    upregulated_metabolites_list = ui["upregulated_metabolites_list"] if "upregulated_metabolites_list" in ui else ""
    downregulated_metabolites_list = ui["downregulated_metabolites_list"] if "downregulated_metabolites_list" in ui else ""
    all_transcripts_list = ui["all_transcripts_list"] if "all_transcripts_list" in ui else ""
    regulated_transcripts_list = ui["regulated_transcripts_list"] if "regulated_transcripts_list" in ui else ""
    upregulated_transcripts_list = ui["upregulated_transcripts_list"] if "upregulated_transcripts_list" in ui else ""
    downregulated_transcripts_list = ui["downregulated_transcripts_list"] if "downregulated_transcripts_list" in ui else ""
    separate_updown_switch = ui["separate_updown_switch"] if "separate_updown_switch" in ui else False
    with_lipids = ui["with_lipids"] if "with_lipids" in ui else False
    with_proteins = ui["with_proteins"] if "with_proteins" in ui else False
    with_metabolites = ui["with_metabolites"] if "with_metabolites" in ui else False
    with_transcripts = ui["with_transcripts"] if "with_transcripts" in ui else False
    select_organism = ui["select_organism"] if "select_organism" in ui else INIT_ORGANISM
    select_molecule_handling = ui["select_molecule_handling"] if "select_molecule_handling" in ui else MOLECULE_HANDLING_REMOVE
    select_regulated_molecule_handling = ui["select_regulated_molecule_handling"] if "select_regulated_molecule_handling" in ui else MOLECULE_HANDLING_REMOVE
    multiselect_filter_molecules = ui["multiselect_filter_molecules"] if "multiselect_filter_molecules" in ui else []
    select_term_representation = ui["select_term_representation"] if "select_term_representation" in ui else "greater"
    select_test_method = ui["select_test_method"] if "select_test_method" in ui else "fdr_bh"
    select_domains = ui["select_domains"] if "select_domains" in ui else ["Biological process"]

    num_all_lipids = sum(len(line) > 0 for line in all_lipids_list.split("\n"))
    num_regulated_lipids = sum(len(line) > 0 for line in regulated_lipids_list.split("\n"))
    num_upregulated_lipids = sum(len(line) > 0 for line in upregulated_lipids_list.split("\n"))
    num_downregulated_lipids = sum(len(line) > 0 for line in downregulated_lipids_list.split("\n"))
    num_all_proteins = sum(len(line) > 0 for line in all_proteins_list.split("\n"))
    num_regulated_proteins = sum(len(line) > 0 for line in regulated_proteins_list.split("\n"))
    num_upregulated_proteins = sum(len(line) > 0 for line in upregulated_proteins_list.split("\n"))
    num_downregulated_proteins = sum(len(line) > 0 for line in downregulated_proteins_list.split("\n"))
    num_all_metabolites = sum(len(line) > 0 for line in all_metabolites_list.split("\n"))
    num_regulated_metabolites = sum(len(line) > 0 for line in regulated_metabolites_list.split("\n"))
    num_upregulated_metabolites = sum(len(line) > 0 for line in upregulated_metabolites_list.split("\n"))
    num_downregulated_metabolites = sum(len(line) > 0 for line in downregulated_metabolites_list.split("\n"))
    num_all_transcripts = sum(len(line) > 0 for line in all_transcripts_list.split("\n"))
    num_regulated_transcripts = sum(len(line) > 0 for line in regulated_transcripts_list.split("\n"))
    num_upregulated_transcripts = sum(len(line) > 0 for line in upregulated_transcripts_list.split("\n"))
    num_downregulated_transcripts = sum(len(line) > 0 for line in downregulated_transcripts_list.split("\n"))

    return (
        session_id,
        no_update,
        no_update,
        separate_updown_switch,
        f"Entries: {num_all_lipids}",
        f"Entries: {num_regulated_lipids}",
        f"Entries: {num_upregulated_lipids}",
        f"Entries: {num_downregulated_lipids}",
        f"Entries: {num_all_proteins}",
        f"Entries: {num_regulated_proteins}",
        f"Entries: {num_upregulated_proteins}",
        f"Entries: {num_downregulated_proteins}",
        f"Entries: {num_all_metabolites}",
        f"Entries: {num_regulated_metabolites}",
        f"Entries: {num_upregulated_metabolites}",
        f"Entries: {num_downregulated_metabolites}",
        f"Entries: {num_all_transcripts}",
        f"Entries: {num_regulated_transcripts}",
        f"Entries: {num_upregulated_transcripts}",
        f"Entries: {num_downregulated_transcripts}",
        all_lipids_list,
        regulated_lipids_list,
        upregulated_lipids_list,
        downregulated_lipids_list,
        all_proteins_list,
        regulated_proteins_list,
        upregulated_proteins_list,
        downregulated_proteins_list,
        all_metabolites_list,
        regulated_metabolites_list,
        upregulated_metabolites_list,
        downregulated_metabolites_list,
        all_transcripts_list,
        regulated_transcripts_list,
        upregulated_transcripts_list,
        downregulated_transcripts_list,
        use_bounded_fatty_acyls,
        with_lipids,
        with_proteins,
        with_metabolites,
        with_transcripts,
        select_organism,
        select_molecule_handling,
        select_regulated_molecule_handling,
        multiselect_filter_molecules,
        select_term_representation,
        select_test_method,
        select_domains,
    )



@callback(
    Output("loading_output", "children", allow_duplicate = True),
    Output("graph_enrichment_results", "rowData", allow_duplicate = True),
    Output("graph_enrichment_results", "selectedRows", allow_duplicate = True),
    Output("graph_enrichment_results", "filterModel", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("histogram_results", "disabled", allow_duplicate = True),
    Output("multiselect_filter_molecules", "data", allow_duplicate = True),
    Output("multiselect_filter_molecules", "value", allow_duplicate = True),
    Output("num_enrichment_terms", "children", allow_duplicate = True),
    Output("graph_enrichment_results", "columnDefs", allow_duplicate = True),
    Output("analytics", "children", allow_duplicate = True),
    Input("button_run_enrichment", "n_clicks"),
    State("textarea_all_lipids", "value"),
    State("textarea_regulated_lipids", "value"),
    State("textarea_upregulated_lipids", "value"),
    State("textarea_downregulated_lipids", "value"),
    State("textarea_all_proteins", "value"),
    State("textarea_regulated_proteins", "value"),
    State("textarea_upregulated_proteins", "value"),
    State("textarea_downregulated_proteins", "value"),
    State("textarea_all_metabolites", "value"),
    State("textarea_regulated_metabolites", "value"),
    State("textarea_upregulated_metabolites", "value"),
    State("textarea_downregulated_metabolites", "value"),
    State("textarea_all_transcripts", "value"),
    State("textarea_regulated_transcripts", "value"),
    State("textarea_upregulated_transcripts", "value"),
    State("textarea_downregulated_transcripts", "value"),
    State("select_organism", "value"),
    State("select_domains", "value"),
    State("select_molecule_handling", "value"),
    State("select_test_method", "value"),
    State("select_regulated_molecule_handling", "value"),
    State("select_term_representation", "value"),
    State("sessionid", "children"),
    State("checkbox_use_lipids", "checked"),
    State("checkbox_use_proteins", "checked"),
    State("checkbox_use_metabolites", "checked"),
    State("checkbox_use_transcripts", "checked"),
    State("separate_updown_switch", "checked"),
    prevent_initial_call = True,
)
def run_enrichment(
    n_clicks,
    all_lipids_list,
    regulated_lipids_list,
    upregulated_lipids_list,
    downregulated_lipids_list,
    all_proteins_list,
    regulated_proteins_list,
    upregulated_proteins_list,
    downregulated_proteins_list,
    all_metabolites_list,
    regulated_metabolites_list,
    upregulated_metabolites_list,
    downregulated_metabolites_list,
    all_transcripts_list,
    regulated_transcripts_list,
    upregulated_transcripts_list,
    downregulated_transcripts_list,
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
    with_transcripts,
    separate_updown_switch,
):
    histogram_disabled = True
    num_enrichment_terms = "Entries: 0"

    if session_id not in sessions:
        return "", [], [], {}, True, "Your session has expired. Please refresh the website.",  histogram_disabled, [], [], no_update, num_enrichment_terms

    logger.info(f"Enrichment session: {session_id}")
    session = sessions[session_id]
    session.time = time.time()

    if not with_lipids and not with_proteins and not with_metabolites and not with_transcripts:
        return "", [], [], {}, True, "No omics data is selected.", histogram_disabled, [], [], num_enrichment_terms, no_update, no_update

    if len(domains) == 0:
        return "", [], [], {}, True, "No domain(s) selected.", histogram_disabled, [], [], num_enrichment_terms, no_update, no_update

    if organism is None or organism not in enrichment_ontologies:
        return "", [], [], {}, True, "No organism selected.", histogram_disabled, [], [], num_enrichment_terms, no_update, no_update
    
    ontology = enrichment_ontologies[organism]
    target_set = set()
    lipidome, regulated_lipids = {}, set()
    proteome, regulated_proteins = set(), set()
    metabolome, regulated_metabolites = set(), set()
    transcriptome, regulated_transcripts = set(), set()
    background_list = []
    data_not_loaded = not session.data_loaded

    if data_not_loaded:
        omics_included = [with_lipids, with_proteins, with_metabolites, with_transcripts]

        omics_lists = [
            all_lipids_list,
            regulated_lipids_list,
            upregulated_lipids_list,
            downregulated_lipids_list,
            all_proteins_list,
            regulated_proteins_list,
            upregulated_proteins_list,
            downregulated_proteins_list,
            all_metabolites_list,
            regulated_metabolites_list,
            upregulated_metabolites_list,
            downregulated_metabolites_list,
            all_transcripts_list,
            regulated_transcripts_list,
            upregulated_transcripts_list,
            downregulated_transcripts_list,
        ]

        try:
            (
                target_set,
                lipidome,
                regulated_lipids,
                upregulated_lipids,
                downregulated_lipids,
                proteome,
                regulated_proteins,
                upregulated_proteins,
                downregulated_proteins,
                metabolome,
                regulated_metabolites,
                upregulated_metabolites,
                downregulated_metabolites,
                transcriptome,
                regulated_transcripts,
                upregulated_transcripts,
                downregulated_transcripts,
                background_list,
            ) = check_user_input(
                separate_updown_switch,
                omics_included,
                omics_lists,
                ontology,
                ignore_unrecognizable_molecules,
                ignore_unknown,
            )
        except Exception as error_message:
            return "", [], [], {}, True, str(error_message), histogram_disabled, [], [], num_enrichment_terms, no_update, no_update

        (
            session.search_terms,
            session.all_parent_nodes,
        ) = ontology.set_background(
            lipid_dict = lipidome,
            protein_set = proteome,
            metabolite_set = metabolome,
            transcript_set = transcriptome,
            use_bounded_fatty_acyls = session.use_bounded_fatty_acyls,
        )
        session.num_background = len(lipidome) + len(proteome) + len(metabolome) + len(transcriptome)
        session.ontology = ontology
        session.data_loaded = True

        session.background_lipids = lipidome if with_lipids else None
        session.regulated_lipids = regulated_lipids if with_lipids else None
        session.upregulated_lipids = upregulated_lipids if with_lipids else None
        session.downregulated_lipids = downregulated_lipids if with_lipids else None
        session.background_proteins = proteome if with_proteins else None
        session.regulated_proteins = regulated_proteins if with_proteins else None
        session.upregulated_proteins = upregulated_proteins if with_proteins else None
        session.downregulated_proteins = downregulated_proteins if with_proteins else None
        session.background_metabolites = metabolome if with_metabolites else None
        session.regulated_metabolites = regulated_metabolites if with_metabolites else None
        session.upregulated_metabolites = upregulated_metabolites if with_metabolites else None
        session.downregulated_metabolites = downregulated_metabolites if with_metabolites else None
        session.background_transcripts = transcriptome if with_transcripts else None
        session.regulated_transcripts = regulated_transcripts if with_transcripts else None
        session.upregulated_transcripts = upregulated_transcripts if with_transcripts else None
        session.downregulated_transcripts = downregulated_transcripts if with_transcripts else None
        background_list.sort(key = lambda row: row["value"])
        session.background_list = background_list

    else:
        target_set = set()
        if with_lipids: target_set |= session.regulated_lipids
        if with_proteins: target_set |= session.regulated_proteins
        if with_metabolites: target_set |= session.regulated_metabolites
        if with_transcripts: target_set |= session.regulated_transcripts
        background_list = no_update

    session.domains = set(domains)
    session.separate_updown_switch = separate_updown_switch
    results = ontology.enrichment_analysis(separate_updown_switch, session.search_terms, session.num_background, target_set, domains, term_regulation, correction_method)
    session.result = results

    data = []
    session.data = {}
    session.results = []
    result_map = {}
    if separate_updown_switch:
        for result in results:
            expected_up = round((result.fisher_data[0][0] + result.fisher_data[0][1]) * (result.fisher_data[0][0] + result.fisher_data[0][2]) / sum(result.fisher_data[0]))
            expected_down = round((result.fisher_data[1][0] + result.fisher_data[1][1]) * (result.fisher_data[1][0] + result.fisher_data[1][2]) / sum(result.fisher_data[1]))
            direction = ""
            if result.pvalue_corrected[0] < 0.05:
                if 0 < result.pvalue_corrected[1] < 0.05 and 0 < result.pvalue_corrected[2] < 0.05:
                    direction = "⇅"
                elif 0 < result.pvalue_corrected[1] < 0.05:
                    direction = "↑"
                elif 0 < result.pvalue_corrected[2] < 0.05:
                    direction = "↓"

            row = {
                "domain": " | ".join(result.term.domain),
                "term": result.term.name,
                "termid": result.term.term_id_str,
                "count": f"{result.fisher_data[0][0]}, {result.fisher_data[1][0]} ({expected_up}, {expected_down}) / {len(result.source_terms)}",
                "direction": direction,
                "pvalue": result.pvalue_corrected[0],
                "log_odds_ratio": f"{result.lor[0]:.3f} / {result.lor[1]:.3f}",
            }
            data.append(row)
            session.results.append((result, row))

    else:
        for result in results:
            expected = round((result.fisher_data[0] + result.fisher_data[1]) * (result.fisher_data[0] + result.fisher_data[2]) / sum(result.fisher_data))
            row = {
                "domain": " | ".join(result.term.domain),
                "term": result.term.name,
                "termid": result.term.term_id_str,
                "direction": "",
                "count": f"{result.fisher_data[0]} ({expected}) / {len(result.source_terms)}",
                "pvalue": result.pvalue_corrected,
                "log_odds_ratio": result.lor,
            }
            data.append(row)
            session.results.append((result, row))

    histogram_disabled = False
    for result in results:
        for term_id in result.term.term_id:
            session.data[term_id] = result

    results_def = graph_enrichment_results_def(separate_updown_switch, len(domains) > 1, correction_method != "no")
    num_enrichment_terms = f"Entries: {len(results)}"

    return (
        "",
        data,
        [],
        {},
        no_update,
        no_update,
        histogram_disabled,
        background_list,
        [],
        num_enrichment_terms,
        results_def,
        "enrichment_analysis",
    )



@callback(
    Output("graph_enrichment_results", "rowData", allow_duplicate = True),
    Output("graph_enrichment_results", "selectedRows", allow_duplicate = True),
    Output("graph_enrichment_results", "filterModel", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("num_enrichment_terms", "children", allow_duplicate = True),
    Input("multiselect_filter_molecules", "value"),
    State("multiselect_filter_molecules", "data"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def filter_result_table(multiselect_values, multiselect_data, session_id):
    if multiselect_data == None or len(multiselect_data) == 0 or multiselect_values == None or session_id == None:
        raise exceptions.PreventUpdate

    num_enrichment_terms = "Entries: 0"
    if session_id not in sessions:
        return (
            no_update,
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            num_enrichment_terms,
        )

    session = sessions[session_id]
    session.time = time.time()

    if session.results is None or len(session.results) == 0:
        raise exceptions.PreventUpdate

    if len(multiselect_values) > 0:
        multiselect_values = set(multiselect_values)
        data = [row for result, row in session.results if not set(result.source_terms).isdisjoint(multiselect_values)]
    else:
        data = [row for _, row in session.results]
    return data, [], {}, no_update, no_update, f"Entries: {len(data)}"



@callback(
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("num_all_lipids", "children", allow_duplicate = True),
    Output("num_regulated_lipids", "children", allow_duplicate = True),
    Output("num_upregulated_lipids", "children", allow_duplicate = True),
    Output("num_downregulated_lipids", "children", allow_duplicate = True),
    Output("num_all_proteins", "children", allow_duplicate = True),
    Output("num_regulated_proteins", "children", allow_duplicate = True),
    Output("num_upregulated_proteins", "children", allow_duplicate = True),
    Output("num_downregulated_proteins", "children", allow_duplicate = True),
    Output("num_all_metabolites", "children", allow_duplicate = True),
    Output("num_regulated_metabolites", "children", allow_duplicate = True),
    Output("num_upregulated_metabolites", "children", allow_duplicate = True),
    Output("num_downregulated_metabolites", "children", allow_duplicate = True),
    Output("num_all_transcripts", "children", allow_duplicate = True),
    Output("num_regulated_transcripts", "children", allow_duplicate = True),
    Output("num_upregulated_transcripts", "children", allow_duplicate = True),
    Output("num_downregulated_transcripts", "children", allow_duplicate = True),
    Input("textarea_all_lipids", "value"),
    Input("textarea_regulated_lipids", "value"),
    Input("textarea_upregulated_lipids", "value"),
    Input("textarea_downregulated_lipids", "value"),
    Input("textarea_all_proteins", "value"),
    Input("textarea_regulated_proteins", "value"),
    Input("textarea_upregulated_proteins", "value"),
    Input("textarea_downregulated_proteins", "value"),
    Input("textarea_all_metabolites", "value"),
    Input("textarea_regulated_metabolites", "value"),
    Input("textarea_upregulated_metabolites", "value"),
    Input("textarea_downregulated_metabolites", "value"),
    Input("textarea_all_transcripts", "value"),
    Input("textarea_regulated_transcripts", "value"),
    Input("textarea_upregulated_transcripts", "value"),
    Input("textarea_downregulated_transcripts", "value"),
    Input("use_bounded_fatty_acyls", "checked"),
    Input("checkbox_use_lipids", "checked"),
    Input("checkbox_use_proteins", "checked"),
    Input("checkbox_use_metabolites", "checked"),
    Input("checkbox_use_transcripts", "checked"),
    Input("select_organism", "value"),
    Input("select_molecule_handling", "value"),
    Input("select_regulated_molecule_handling", "value"),
    Input("multiselect_filter_molecules", "value"),
    Input("select_term_representation", "value"),
    Input("select_test_method", "value"),
    Input("select_domains", "value"),
    Input("separate_updown_switch", "checked"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def update_background(
    all_lipids_list,
    regulated_lipids_list,
    upregulated_lipids_list,
    downregulated_lipids_list,
    all_proteins_list,
    regulated_proteins_list,
    upregulated_proteins_list,
    downregulated_proteins_list,
    all_metabolites_list,
    regulated_metabolites_list,
    upregulated_metabolites_list,
    downregulated_metabolites_list,
    all_transcripts_list,
    regulated_transcripts_list,
    upregulated_transcripts_list,
    downregulated_transcripts_list,
    use_bounded_fatty_acyls,
    with_lipids,
    with_proteins,
    with_metabolites,
    with_transcripts,
    select_organism,
    select_molecule_handling,
    select_regulated_molecule_handling,
    multiselect_filter_molecules,
    select_term_representation,
    select_test_method,
    select_domains,
    separate_updown_switch,
    session_id,
):
    num_all_lipids = sum(len(line) > 0 for line in all_lipids_list.split("\n"))
    num_regulated_lipids = sum(len(line) > 0 for line in regulated_lipids_list.split("\n"))
    num_upregulated_lipids = sum(len(line) > 0 for line in upregulated_lipids_list.split("\n"))
    num_downregulated_lipids = sum(len(line) > 0 for line in downregulated_lipids_list.split("\n"))
    num_all_proteins = sum(len(line) > 0 for line in all_proteins_list.split("\n"))
    num_regulated_proteins = sum(len(line) > 0 for line in regulated_proteins_list.split("\n"))
    num_upregulated_proteins = sum(len(line) > 0 for line in upregulated_proteins_list.split("\n"))
    num_downregulated_proteins = sum(len(line) > 0 for line in downregulated_proteins_list.split("\n"))
    num_all_metabolites = sum(len(line) > 0 for line in all_metabolites_list.split("\n"))
    num_regulated_metabolites = sum(len(line) > 0 for line in regulated_metabolites_list.split("\n"))
    num_upregulated_metabolites = sum(len(line) > 0 for line in upregulated_metabolites_list.split("\n"))
    num_downregulated_metabolites = sum(len(line) > 0 for line in downregulated_metabolites_list.split("\n"))
    num_all_transcripts = sum(len(line) > 0 for line in all_transcripts_list.split("\n"))
    num_regulated_transcripts = sum(len(line) > 0 for line in regulated_transcripts_list.split("\n"))
    num_upregulated_transcripts = sum(len(line) > 0 for line in upregulated_transcripts_list.split("\n"))
    num_downregulated_transcripts = sum(len(line) > 0 for line in downregulated_transcripts_list.split("\n"))

    if session_id not in sessions:
        return (
            True,
            "Your session has expired. Please refresh the website.",
            f"Entries: {num_all_lipids}",
            f"Entries: {num_regulated_lipids}",
            f"Entries: {num_upregulated_lipids}",
            f"Entries: {num_downregulated_lipids}",
            f"Entries: {num_all_proteins}",
            f"Entries: {num_regulated_proteins}",
            f"Entries: {num_upregulated_proteins}",
            f"Entries: {num_downregulated_proteins}",
            f"Entries: {num_all_metabolites}",
            f"Entries: {num_regulated_metabolites}",
            f"Entries: {num_upregulated_metabolites}",
            f"Entries: {num_downregulated_metabolites}",
            f"Entries: {num_all_transcripts}",
            f"Entries: {num_regulated_transcripts}",
            f"Entries: {num_upregulated_transcripts}",
            f"Entries: {num_downregulated_transcripts}",
        )
    session = sessions[session_id]
    session.time = time.time()
    session.data_loaded = False
    session.use_bounded_fatty_acyls = use_bounded_fatty_acyls
    session.ui["all_lipids_list"] = all_lipids_list
    session.ui["regulated_lipids_list"] = regulated_lipids_list
    session.ui["upregulated_lipids_list"] = upregulated_lipids_list
    session.ui["downregulated_lipids_list"] = downregulated_lipids_list
    session.ui["all_proteins_list"] = all_proteins_list
    session.ui["regulated_proteins_list"] = regulated_proteins_list
    session.ui["upregulated_proteins_list"] = upregulated_proteins_list
    session.ui["downregulated_proteins_list"] = downregulated_proteins_list
    session.ui["all_metabolites_list"] = all_metabolites_list
    session.ui["regulated_metabolites_list"] = regulated_metabolites_list
    session.ui["upregulated_metabolites_list"] = upregulated_metabolites_list
    session.ui["downregulated_metabolites_list"] = downregulated_metabolites_list
    session.ui["all_transcripts_list"] = all_transcripts_list
    session.ui["regulated_transcripts_list"] = regulated_transcripts_list
    session.ui["upregulated_transcripts_list"] = upregulated_transcripts_list
    session.ui["downregulated_transcripts_list"] = downregulated_transcripts_list
    session.ui["with_lipids"] = with_lipids
    session.ui["with_proteins"] = with_proteins
    session.ui["with_metabolites"] = with_metabolites
    session.ui["with_transcripts"] = with_transcripts
    session.ui["select_organism"] = select_organism
    session.ui["select_molecule_handling"] = select_molecule_handling
    session.ui["select_regulated_molecule_handling"] = select_regulated_molecule_handling
    session.ui["multiselect_filter_molecules"] = multiselect_filter_molecules
    session.ui["select_term_representation"] = select_term_representation
    session.ui["select_test_method"] = select_test_method
    session.ui["select_domains"] = select_domains
    session.ui["separate_updown_switch"] = separate_updown_switch

    return (
        no_update,
        no_update,
        f"Entries: {num_all_lipids}",
        f"Entries: {num_regulated_lipids}",
        f"Entries: {num_upregulated_lipids}",
        f"Entries: {num_downregulated_lipids}",
        f"Entries: {num_all_proteins}",
        f"Entries: {num_regulated_proteins}",
        f"Entries: {num_upregulated_proteins}",
        f"Entries: {num_downregulated_proteins}",
        f"Entries: {num_all_metabolites}",
        f"Entries: {num_regulated_metabolites}",
        f"Entries: {num_upregulated_metabolites}",
        f"Entries: {num_downregulated_metabolites}",
        f"Entries: {num_all_transcripts}",
        f"Entries: {num_regulated_transcripts}",
        f"Entries: {num_upregulated_transcripts}",
        f"Entries: {num_downregulated_transcripts}",
    )



@callback(
    Output("loading_output", "children", allow_duplicate = True),
    Output("download_data", "data", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("analytics", "children", allow_duplicate = True),
    Input("icon_download_results", "n_clicks"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def download_results_table(
    _,
    graph_enrichment_results,
    selected_rows,
    session_id,
):
    if session_id not in sessions:
        return (
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            no_update,
        )

    session = sessions[session_id]
    session.time = time.time()
    domains = []
    term_ids = []
    terms = []
    counts = []
    pvalues = []
    lors = []

    selected_term_ids = {row["termid"] for row in selected_rows}
    if len(graph_enrichment_results) == 0 or len(selected_term_ids) == 0:
        raise exceptions.PreventUpdate

    for row in graph_enrichment_results:
        term_id = row["termid"]
        if term_id not in selected_term_ids: continue
        domains.append(row["domain"])
        term_ids.append(term_id)
        terms.append(row["term"])
        counts.append(row["count"])
        pvalues.append(row["pvalue"])
        lors.append(row["log_odds_ratio"])

    df = pd.DataFrame({"Domain": domains, "Term ID": term_ids, "Term": terms, "Count": counts, "p-value": pvalues, "log-odds-ration": lors})
    data = session.data

    background_lipids = session.background_lipids
    regulated_lipid_set = session.regulated_lipids
    background_proteins = session.background_proteins
    regulated_protein_set = session.regulated_proteins
    background_metabolites = session.background_metabolites
    regulated_metabolite_set = session.regulated_metabolites
    background_transcripts = session.background_transcripts
    regulated_transcript_set = session.regulated_transcripts

    with_lipids = background_lipids != None
    with_proteins = background_proteins != None
    with_metabolites = background_metabolites != None
    with_transcripts = background_transcripts != None

    if with_lipids:
        background_lipids = sorted(list(background_lipids.keys()))
        regulated_lipids = sorted(list(regulated_lipid_set))
        associated_lipids = {}

    if with_proteins:
        background_proteins = sorted(list(background_proteins))
        regulated_proteins = sorted(list(regulated_protein_set))
        associated_proteins = {}

    if with_metabolites:
        background_metabolites = sorted(list(background_metabolites))
        regulated_metabolites = sorted(list(regulated_metabolite_set))
        associated_metabolites = {}

    if with_transcripts:
        background_transcripts = sorted(list(background_transcripts))
        regulated_transcripts = sorted(list(regulated_transcript_set))
        associated_transcripts = {}

    for term_id, term in zip(term_ids, terms):
        term_id = term_id.split("|")[0]
        if term_id not in data: continue
        molecules = data[term_id].source_terms

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

        if with_transcripts:
            transcripts = sorted(list(molecules & set(background_transcripts)))
            associated_transcripts[term] = [transcript for transcript in transcripts]
            associated_transcripts[f"{term}, regulated"] = ["X" if transcript in regulated_transcript_set else "" for transcript in transcripts]

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

    if with_transcripts:
        background_transcripts = [transcript for transcript in background_transcripts]
        regulated_transcripts = [transcript for transcript in regulated_transcripts]
        pd.DataFrame.from_dict(associated_transcripts, orient = "index").transpose().to_excel(writer, sheet_name = "Associated transcripts", index = False)
        pd.DataFrame({"Ensembl": background_transcripts}).to_excel(writer, sheet_name = "Background transcripts", index = False)
        pd.DataFrame({"Ensembl": regulated_transcripts}).to_excel(writer, sheet_name = "Regulated transcripts", index = False)

    writer._save()

    return (
        "",
        dcc.send_bytes(output.getvalue(), "GO_multiomics_results.xlsx"),
        no_update,
        no_update,
        "download_results",
    )



@callback(
    Output("loading_output", "children", allow_duplicate = True),
    Output("download_data", "data", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("analytics", "children", allow_duplicate = True),
    Input("icon_download_sankey_entries", "n_clicks"),
    State("sankey_entries", "virtualRowData"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def download_entries_table(_, sankey_entries_results, session_id):
    if session_id not in sessions:
        return (
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            no_update,
        )

    session = sessions[session_id]
    session.time = time.time()
    ontology = session.ontology
    domains = []
    term_ids = []
    terms = []
    counts = []
    pvalues = []
    lors = []

    if len(sankey_entries_results) == 0:
        raise exceptions.PreventUpdate

    for entry in sankey_entries_results:
        term_id = entry["termid"]
        entry["term name"] = "" if term_id not in ontology.ontology_terms else ontology.ontology_terms[term_id].name

    df = pd.DataFrame(sankey_entries_results)
    output = io.BytesIO()
    writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
    df.to_excel(writer, sheet_name = "Sankey entries", index = False)
    writer._save()

    return (
        "",
        dcc.send_bytes(output.getvalue(), "Sankey_entries.xlsx"),
        no_update,
        no_update,
        "download_sankey_entries",
    )



@callback(
    Output("download_data", "data", allow_duplicate = True),
    Output("html_download_trigger", "value", allow_duplicate = True),
    Input("html_download_trigger", "value"),
    prevent_initial_call = True
)
def download_html_figure(json_figure):
    if json_figure == None or len(json_figure) == 0:
        raise exceptions.PreventUpdate
    figure = go.Figure(json.loads(json_figure))
    buffer = io.StringIO()
    figure.write_html(buffer)
    html_code = buffer.getvalue()
    return dict(content = html_code, filename = "enrichment_figure.html"), ""



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
    Output("textarea_upregulated_lipids", "value", allow_duplicate = True),
    Output("textarea_downregulated_lipids", "value", allow_duplicate = True),
    Output("textarea_all_proteins", "value", allow_duplicate = True),
    Output("textarea_regulated_proteins", "value", allow_duplicate = True),
    Output("textarea_upregulated_proteins", "value", allow_duplicate = True),
    Output("textarea_downregulated_proteins", "value", allow_duplicate = True),
    Output("textarea_all_metabolites", "value", allow_duplicate = True),
    Output("textarea_regulated_metabolites", "value", allow_duplicate = True),
    Output("textarea_upregulated_metabolites", "value", allow_duplicate = True),
    Output("textarea_downregulated_metabolites", "value", allow_duplicate = True),
    Output("textarea_all_transcripts", "value", allow_duplicate = True),
    Output("textarea_regulated_transcripts", "value", allow_duplicate = True),
    Output("textarea_upregulated_transcripts", "value", allow_duplicate = True),
    Output("textarea_downregulated_transcripts", "value", allow_duplicate = True),
    Output("select_organism", "value", allow_duplicate = True),
    Output("checkbox_use_lipids", "checked", allow_duplicate = True),
    Output("checkbox_use_proteins", "checked", allow_duplicate = True),
    Output("checkbox_use_metabolites", "checked", allow_duplicate = True),
    Output("checkbox_use_transcripts", "checked", allow_duplicate = True),
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
        "\n".join(examples[index]["upregl"]),
        "\n".join(examples[index]["downregl"]),
        "\n".join(examples[index]["bgp"]),
        "\n".join(examples[index]["regp"]),
        "\n".join(examples[index]["upregp"]),
        "\n".join(examples[index]["downregp"]),
        "\n".join(examples[index]["bgm"]),
        "\n".join(examples[index]["regm"]),
        "\n".join(examples[index]["upregm"]),
        "\n".join(examples[index]["downregm"]),
        "\n".join(examples[index]["bgt"]),
        "\n".join(examples[index]["regt"]),
        "\n".join(examples[index]["upregt"]),
        "\n".join(examples[index]["downregt"]),
        examples[index]["org"],
        len(examples[index]["bgl"]) > 0 or len(examples[index]["regl"]) > 0 or len(examples[index]["upregl"]) > 0 or len(examples[index]["downregl"]) > 0,
        len(examples[index]["bgp"]) > 0 or len(examples[index]["regp"]) > 0 or len(examples[index]["upregp"]) > 0 or len(examples[index]["downregp"]) > 0,
        len(examples[index]["bgm"]) > 0 or len(examples[index]["regm"]) > 0 or len(examples[index]["upregm"]) > 0 or len(examples[index]["downregm"]) > 0,
        len(examples[index]["bgt"]) > 0 or len(examples[index]["regt"]) > 0 or len(examples[index]["upregt"]) > 0 or len(examples[index]["downregt"]) > 0,
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
    Output("sunburst_results", "disabled", allow_duplicate = True),
    Output("sankey_results", "disabled", allow_duplicate = True),
    Output("icon_download_results", "disabled", allow_duplicate = True),
    Input("graph_enrichment_results", "selectedRows"),
    prevent_initial_call = True,
)
def update_action_icons(selected_rows):
    return (len(selected_rows) == 0, ) * 4



@callback(
    Output("term_molecules_modal", "opened", allow_duplicate = True),
    Output("term_molecules_modal", "title", allow_duplicate = True),
    Output("term_lipids_modal_grid", "rowData", allow_duplicate = True),
    Output("term_proteins_modal_grid", "rowData", allow_duplicate = True),
    Output("term_metabolites_modal_grid", "rowData", allow_duplicate = True),
    Output("term_transcripts_modal_grid", "rowData", allow_duplicate = True),
    Output("term_molecules_modal_id", "children", allow_duplicate = True),
    Output("lipid_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("protein_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("metabolite_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("transcript_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("term_molecules_modal_tab", "value", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("term_path_area", "children", allow_duplicate = True),
    Output("contingency_table", "children", allow_duplicate = True),
    Input("graph_enrichment_results", "cellRendererData"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def open_term_window(
    row_data,
    session_id,
):
    if session_id not in sessions or "rowId" not in row_data:
        return (
            no_update,
            no_update,
            no_update,
            no_update,
            no_update,
            no_update,
            "",
            no_update,
            no_update,
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            "",
            "",
        )
    session = sessions[session_id]
    ontology = session.ontology
    session.time = time.time()
    if "rowId" not in row_data:
        raise exceptions.PreventUpdate

    term_id = row_data["rowId"].split("|")[0]
    if term_id not in session.data:
        raise exceptions.PreventUpdate

    background_lipids = session.background_lipids
    regulated_lipids = session.regulated_lipids
    upregulated_lipids = session.upregulated_lipids
    downregulated_lipids = session.downregulated_lipids
    background_proteins = session.background_proteins
    regulated_proteins = session.regulated_proteins
    upregulated_proteins = session.upregulated_proteins
    downregulated_proteins = session.downregulated_proteins
    background_metabolites = session.background_metabolites
    regulated_metabolites = session.regulated_metabolites
    upregulated_metabolites = session.upregulated_metabolites
    downregulated_metabolites = session.downregulated_metabolites
    background_transcripts = session.background_transcripts
    regulated_transcripts = session.regulated_transcripts
    upregulated_transcripts = session.upregulated_transcripts
    downregulated_transcripts = session.downregulated_transcripts

    with_lipids = background_lipids != None
    with_proteins = background_proteins != None
    with_metabolites = background_metabolites != None
    with_transcripts = background_transcripts != None

    lipid_table = []
    protein_table = []
    metabolite_table = []
    transcript_table = []

    result = sessions[session_id].data[term_id]
    molecules = sorted(list(result.source_terms))

    def molecule_regulated(molecule, m_type):
        match m_type:
            case "l":
                regulated_molecules = regulated_lipids
                upregulated_molecules = upregulated_lipids
                downregulated_molecules = downregulated_lipids
            case "p":
                regulated_molecules = regulated_proteins
                upregulated_molecules = upregulated_proteins
                downregulated_molecules = downregulated_proteins
            case "m":
                regulated_molecules = regulated_metabolites
                upregulated_molecules = upregulated_metabolites
                downregulated_molecules = downregulated_metabolites
            case "t":
                regulated_molecules = regulated_transcripts
                upregulated_molecules = upregulated_transcripts
                downregulated_molecules = downregulated_transcripts

        if session.separate_updown_switch:
            if molecule in upregulated_molecules: return "↑"
            if molecule in downregulated_molecules: return "↓"
            return ""
        else:
            if molecule in regulated_molecules: return "x"
            return ""

    if with_lipids:
        background_lipids = set(background_lipids.keys())
        lipid_table = [{"molecule_id": molecule, "molecule": molecule, "regulated": molecule_regulated(molecule, "l")} for molecule in molecules if molecule in background_lipids]

    if with_proteins:
        protein_table = [
            {
                "molecule_id": molecule,
                "molecule": molecule.replace("UNIPROT:", "") + f" ({ontology.proteins[molecule].name})" if molecule in ontology.proteins else "",
                "regulated": molecule_regulated(molecule, "p")
            } for molecule in molecules if molecule in background_proteins
        ]

    if with_metabolites:
        metabolite_table = [{"molecule_id": molecule, "molecule": molecule, "regulated": molecule_regulated(molecule, "m")} for molecule in molecules if molecule in background_metabolites]

    if with_transcripts:
        transcript_table = [
            {
                "molecule_id": molecule,
                "molecule": molecule + f" ({ontology.transcripts[m].name})" if m in ontology.transcripts else "",
                "regulated": molecule_regulated(molecule, "t")
            } for molecule in molecules if molecule in background_transcripts and (m := molecule.split(".")[0])
        ]

    tab_value = "lipid_tab_modal"
    if not with_lipids:
        tab_value = "protein_tab_modal"
        if not with_proteins:
            tab_value = "metabolite_tab_modal"
            if not with_metabolites:
                tab_value = "transcript_tab_modal"

    col_1 = ""
    col_2 = "Molecules associated with Term"
    col_3 = "Molecules not associated with Term"
    col_4 = "Sum"

    if session.separate_updown_switch:
        fisher_data_up = result.fisher_data[0]
        fisher_data_down = result.fisher_data[1]
        data = {
            col_1: ["Up-regulated", "Down-regulated", "Not regulated", "Sum"],
            col_2: [
                fisher_data_up[0],
                fisher_data_down[0],
                fisher_data_up[1] - fisher_data_down[0],
                fisher_data_up[0] + fisher_data_up[1],
            ],
            col_3: [
                fisher_data_up[2],
                fisher_data_down[2],
                fisher_data_up[3] - fisher_data_down[2],
                fisher_data_up[2] + fisher_data_up[3],
            ],
            col_4: [
                fisher_data_up[0] + fisher_data_up[2],
                fisher_data_down[0] + fisher_data_down[2],
                (fisher_data_up[1] - fisher_data_down[0]) + (fisher_data_up[3] - fisher_data_down[2]),
                sum(fisher_data_up),
            ]
        }

    else:
        fisher_data = result.fisher_data
        data = {
            col_1: ["Regulated molecules", "Not regulated Molecules", "Sum"],
            col_2: [fisher_data[0], fisher_data[1], fisher_data[0] + fisher_data[1]],
            col_3: [fisher_data[2], fisher_data[3], fisher_data[2] + fisher_data[3]],
            col_4: [fisher_data[0] + fisher_data[2], fisher_data[1] + fisher_data[3], fisher_data[0] + fisher_data[1] + fisher_data[2] + fisher_data[3]]
        }
    dfc = pd.DataFrame(data)

    contingency_table = dash_table.DataTable(
        columns=[{"name": col, "id": col} for col in dfc.columns],
        data = dfc.to_dict("records"),
        style_table={'width': '80%', 'margin': 'auto'},

        style_cell_conditional=[
            {'if': {'column_id': col_2}, 'width': '30%'},
            {'if': {'column_id': col_3}, 'width': '30%'},
            {'if': {'column_id': col_4}, 'width': '10%'}
        ],

        # Style header row
        style_header={'backgroundColor': '#f8f8f8', 'fontWeight': 'bold'},

        style_header_conditional=[
            {'if': {'column_id': col_3}, 'borderRight': '1px solid black'}
        ],

        # Style first column to match the header row
        style_data_conditional=[
            {'if': {'column_id': col_1}, 'backgroundColor': '#f8f8f8'},
            {'if': {'row_index': 2 if session.separate_updown_switch else 1}, 'borderBottom': '1px solid black'},
            {'if': {'column_id': col_3}, 'borderRight': '1px solid black'},
            {'if': {'column_id': col_1}, 'fontWeight': 'bold'},
        ]
    )

    return (
        True,
        f"Molecules for '{result.term.name}'",
        lipid_table,
        protein_table,
        metabolite_table,
        transcript_table,
        list(result.term.term_id)[0],
        not with_lipids,
        not with_proteins,
        not with_metabolites,
        not with_transcripts,
        tab_value,
        no_update,
        no_update,
        "",
        contingency_table,
    )



@callback(
    Output("barplot_terms_modal", "opened", allow_duplicate = True),
    Output("barplot_terms", "figure", allow_duplicate = True),
    Output("barplot_terms_modal", "size", allow_duplicate = True),
    Output("barplot_terms_wrapper", "style", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("barplot_controls", "style", allow_duplicate = True),
    Output("sunburst_controls", "style", allow_duplicate = True),
    Output("barplot_numberinput_max_pvalue", "value", allow_duplicate = True),
    Output("barplot_numberinput_min_pvalue", "value", allow_duplicate = True),
    Input("sunburst_results", "n_clicks"),
    Input("barplot_numberinput_max_pvalue", "value"),
    Input("barplot_numberinput_min_pvalue", "value"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    State("sessionid", "children"),
    State("barplot_controls", "style"),
    State("sunburst_controls", "style"),
    State("barplot_terms_wrapper", "style"),
    prevent_initial_call = True,
)
def open_sunburstplot(
    n_clicks,
    pval_max,
    pval_min,
    row_data,
    selected_rows,
    session_id,
    barplot_controls_style,
    sunburst_controls_style,
    barplot_terms_wrapper_style,
):
    if session_id == None or n_clicks == None:
        raise exceptions.PreventUpdate

    if session_id not in sessions:
        return (
            no_update,
            no_update,
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            no_update,
            no_update,
            no_update,
            no_update,
        )

    session = sessions[session_id]

    pval_max_float, pval_min_float = pval_max, pval_min

    try:
        pval_max_float = float(pval_max_float)
    except Exception as e:
        return tuple([no_update] * 8 + [session.max_pvalue, session.min_pvalue])

    try:
        pval_min_float = float(pval_min_float)
    except Exception as e:
        return tuple([no_update] * 8 + [session.max_pvalue, session.min_pvalue])

    if not (0 <= pval_max_float <= 1):
        return tuple([no_update] * 8 + [session.max_pvalue, session.min_pvalue])

    if not (0 <= pval_min_float <= 1):
        return tuple([no_update] * 8 + [session.max_pvalue, session.min_pvalue])

    if pval_max_float < pval_min_float:
        return tuple([no_update] * 8 + [session.max_pvalue, session.min_pvalue])

    session.min_pvalue = pval_min
    session.max_pvalue = pval_max

    barplot_controls_style["display"] = "none"
    sunburst_controls_style["display"] = "block"
    fig = go.Figure()
    ontology = session.ontology
    domains = session.domains

    session_data = session.data
    selected_term_ids = [row["termid"] for row in selected_rows]
    terms = ontology.ontology_terms

    pval_max_float = -np.log10(pval_max_float)
    pval_min_float = -np.log10(pval_min_float)

    sunburst_terms = {}
    for term_id in selected_term_ids:
        queue = [term_id]
        term = terms[term_id.split("|")[0]]
        term_prefix = [t.split(":")[0] if t.find(":") > -1 else "R-" for t in term_id.split("|")]
        term_prefix = [t for t in term_prefix if ((t in LIBRARY_TO_DOMAINS) and (len(LIBRARY_TO_DOMAINS[t] & domains) > 0))]
        if term in sunburst_terms:
            sunburst_terms[term].entry_point = True
            continue
        if len(term_prefix) == 0: continue
        term_prefix = term_prefix[0]

        while len(queue) > 0:
            child_term_id = queue.pop().split("|")[0]
            if child_term_id not in terms: continue
            child_term = terms[child_term_id]
            if child_term in sunburst_terms: continue
            if session.separate_updown_switch:
                pvalue = session_data[child_term_id].pvalue_corrected[0] if child_term_id in session_data else -1
            else:
                pvalue = session_data[child_term_id].pvalue_corrected if child_term_id in session_data else -1
            sunburst_terms[child_term] = SunburstTerm(child_term, term == child_term, pvalue)

            for parent_term in child_term.relations:
                for parent_term_id in parent_term.term_id:
                    parent_prefix = parent_term_id.split(":")[0] if parent_term_id.find(":") > -1 else "R-"
                    if parent_prefix not in LIBRARY_TO_DOMAINS \
                        or len(LIBRARY_TO_DOMAINS[parent_prefix] & domains) == 0 \
                        or term_prefix != parent_prefix \
                        or parent_term_id not in terms: continue

                    parent_term = terms[parent_term_id]
                    sunburst_terms[child_term].parent = parent_term.term_id_str
                    queue.append(sunburst_terms[child_term].parent)
                    break

    labels, parents, colors, values, names, custom_data = [], [], [], [], [], []
    num_letters = 10 + int(80 / len(selected_term_ids))

    for term, sunburst_term in sunburst_terms.items():
        labels.append(term.term_id_str)
        parents.append(sunburst_term.parent)
        values.append(int(sunburst_term.entry_point))
        color = "#cccccc"
        names.append(shorten_label(term.name, num_letters) if sunburst_term.entry_point else "")
        custom_text = f"Id: <b>{labels[-1]}</b><br>Name: {term.name}"
        if sunburst_term.pvalue > 0:
            pvalue = max(min(-np.log10(sunburst_term.pvalue), pval_min_float), pval_max_float)
            ratio = (pvalue - pval_max_float) / (pval_min_float - pval_max_float)
            term_session_data = session_data[list(term.term_id)[0]]
            fisher_data = term_session_data.fisher_data

            if session.separate_updown_switch:
                pvalues = session_data[list(term.term_id)[0]].pvalue_corrected
                n_up, K_up, N_up = fisher_data[0][0] + fisher_data[0][2], fisher_data[0][0] + fisher_data[0][1], sum(fisher_data[0])
                n_down, K_down, N_down = fisher_data[1][0] + fisher_data[1][2], fisher_data[1][0] + fisher_data[1][1], sum(fisher_data[1])
                overrepresented_up = fisher_data[0][0] >= ceil((n_up + 1) * (K_up + 1) / (N_up + 2))
                overrepresented_down = fisher_data[1][0] >= ceil((n_down + 1) * (K_down + 1) / (N_down + 2))
                if (
                    0 < pvalues[1] <= pval_max_float
                    and overrepresented_up
                    and 0 < pvalues[2] <= pval_max_float
                    and overrepresented_down
                ):
                    h, s, l = 308, int(ratio * 100), 75 - int(ratio * 39) # 39 = 75 - 36
                    regulated_metabolites = fisher_data[0][0] + fisher_data[1][0]
                    all_metabolites = regulated_metabolites + fisher_data[0][1] + fisher_data[1][1]
                    custom_text += f"<br>p-value: {sunburst_term.pvalue}<br>Up/Down-regulated molecules: {regulated_metabolites} / {all_metabolites}"

                elif 0 < pvalues[1] <= pval_max_float and overrepresented_up:
                    h, s, l = 3, int(ratio * 100), 75 - int(ratio * 39) # 39 = 75 - 36
                    regulated_metabolites = fisher_data[0][0]
                    all_metabolites = regulated_metabolites + fisher_data[0][1]
                    custom_text += f"<br>p-value: {sunburst_term.pvalue}<br>Up-regulated molecules: {regulated_metabolites} / {all_metabolites}"

                elif 0 < pvalues[2] <= pval_max_float and overrepresented_down:
                    h, s, l = 207, int(ratio * 100), 75 - int(ratio * 39) # 39 = 75 - 36
                    regulated_metabolites = fisher_data[1][0]
                    all_metabolites = regulated_metabolites + fisher_data[1][1]
                    custom_text += f"<br>p-value: {sunburst_term.pvalue}<br>Down-regulated molecules: {regulated_metabolites} / {all_metabolites}"

                else:
                    h, s, l = 0, 0, 80
                    custom_text += f"<br>p-value: {sunburst_term.pvalue}"
            else:
                n, K, N = fisher_data[0] + fisher_data[2], fisher_data[0] + fisher_data[1], sum(fisher_data)
                overrepresented = fisher_data[0] >= ceil((n + 1) * (K + 1) / (N + 2))
                if overrepresented:
                    h, s, l = 3, int(ratio * 100), 75 - int(ratio * 39) # 39 = 75 - 36
                    regulated_metabolites = fisher_data[0]
                    all_metabolites = regulated_metabolites + fisher_data[1]
                    custom_text += f"<br>p-value: {sunburst_term.pvalue}<br>Regulated molecules: {regulated_metabolites} / {all_metabolites}"

                else:
                    h, s, l = 0, 0, 80
                    custom_text += f"<br>p-value: {sunburst_term.pvalue}"

            color = f"hsl({h}, {s}%, {l}%)"


        colors.append(color)
        custom_data.append(custom_text)


    fig.add_trace(go.Sunburst(
        ids = labels,
        parents = parents,
        labels = names,
        values = values,
        marker = dict(
            colors = colors,
            colorscale = None,
            cmin = 0,
            cmax = 1,
            showscale = False,
        ),
        #textinfo = 'none',
        customdata = custom_data,
        hovertemplate="%{customdata}<extra></extra>",
    ))
    fig.update_layout(
        margin = dict(t = 5, l = 5, r = 5, b = 5),
    )
    barplot_terms_wrapper_style["height"] = "70vh"

    return (
        True,
        fig,
        "70%",
        barplot_terms_wrapper_style,
        no_update,
        no_update,
        barplot_controls_style,
        sunburst_controls_style,
        no_update,
        no_update,
    )



@callback(
    Output("sankey_modal", "opened", allow_duplicate = True),
    Output("sankey_graph", "figure", allow_duplicate = True),
    Output("sankey_graph_wrapper", "style", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("sankey_entries", "rowData", allow_duplicate = True),
    Output("icon_download_sankey_entries", "disabled", allow_duplicate = True),
    Output("radiogroup_sankey", "value", allow_duplicate = True),
    Input("sankey_results", "n_clicks"),
    Input("radiogroup_sankey", "value"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    State("sessionid", "children"),
    State("barplot_controls", "style"),
    State("sunburst_controls", "style"),
    State("barplot_terms_wrapper", "style"),
    prevent_initial_call = True,
)
def open_sankeyplot(
    n_clicks,
    radiogroup_sankey,
    row_data,
    selected_rows,
    session_id,
    barplot_controls_style,
    sunburst_controls_style,
    barplot_terms_wrapper_style,
):
    if session_id == None or n_clicks == None or radiogroup_sankey == None:
        raise exceptions.PreventUpdate

    if session_id not in sessions:
        return (
            no_update,
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            no_update,
            no_update,
            no_update,
        )
    session = sessions[session_id]
    selected_term_ids = [row["termid"] for row in selected_rows]
    ontology = session.ontology
    domains = session.domains
    terms = ontology.ontology_terms

    barplot_controls_style["display"] = "none"
    sunburst_controls_style["display"] = "none"

    if ctx.triggered_id != "radiogroup_sankey":
        radiogroup_sankey = "sankey_all"

    background_lipids = session.background_lipids
    regulated_lipids = session.regulated_lipids
    upregulated_lipids = session.upregulated_lipids
    downregulated_lipids = session.downregulated_lipids
    background_proteins = session.background_proteins
    regulated_proteins = session.regulated_proteins
    upregulated_proteins = session.upregulated_proteins
    downregulated_proteins = session.downregulated_proteins
    background_metabolites = session.background_metabolites
    regulated_metabolites = session.regulated_metabolites
    upregulated_metabolites = session.upregulated_metabolites
    downregulated_metabolites = session.downregulated_metabolites
    background_transcripts = session.background_transcripts
    regulated_transcripts = session.regulated_transcripts
    upregulated_transcripts = session.upregulated_transcripts
    downregulated_transcripts = session.downregulated_transcripts

    if session.separate_updown_switch:
        upregulated_molecules = (
            (upregulated_lipids if upregulated_lipids else set()) |
            (upregulated_proteins if upregulated_proteins else set()) |
            (upregulated_metabolites if upregulated_metabolites else set()) |
            (upregulated_transcripts if upregulated_transcripts else set())
        )
        downregulated_molecules = (
            (downregulated_lipids if downregulated_lipids else set()) |
            (downregulated_proteins if downregulated_proteins else set()) |
            (downregulated_metabolites if downregulated_metabolites else set()) |
            (downregulated_transcripts if downregulated_transcripts else set())
        )

        regulated_molecules = upregulated_molecules | downregulated_molecules

    else:
        regulated_molecules = (
            (regulated_lipids if regulated_lipids else set()) |
            (regulated_proteins if regulated_proteins else set()) |
            (regulated_metabolites if regulated_metabolites else set()) |
            (regulated_transcripts if regulated_transcripts else set())
        )
        upregulated_molecules, downregulated_molecules = set(), set()

    pastel_colors = ["#FFD1DC", "#AAF0D1", "#FFB347", "#B5EAAA", "#CBAACB", "#FDFD96", "#CFCFC4", "#E3E4FA", "#AEC6CF", "#FF6961", "#FFE5B4", "#77DD77"]
    def compute_layers(n_nodes, source, target):
        # Build graph
        incoming = defaultdict(set)
        outgoing = defaultdict(set)

        for s, t in zip(source, target):
            outgoing[s].add(t)
            incoming[t].add(s)

        # Kahn's algorithm (topological layering)
        layer = {i: 0 for i in range(n_nodes)}
        queue = deque()

        # Start with nodes with no incoming edges
        for i in range(n_nodes):
            if not incoming[i]:
                queue.append(i)

        while queue:
            node = queue.popleft()
            for child in outgoing[node]:
                layer[child] = max(layer[child], layer[node] + 1)
                incoming[child].remove(node)
                if not incoming[child]:
                    queue.append(child)

        return layer

    sankey_data = {}
    sankey_data_ids = {}
    for target_term_id in selected_term_ids:
        target_term_id_single = target_term_id.split("|")[0]
        target_term = terms[target_term_id_single]

        for molecule in session.data[target_term_id_single].source_terms:
            if (
                (radiogroup_sankey == "sankey_non" and molecule in regulated_molecules)
                or (radiogroup_sankey == "sankey_regulated" and molecule not in regulated_molecules)
                or (radiogroup_sankey == "sankey_upregulated" and molecule not in upregulated_molecules)
                or (radiogroup_sankey == "sankey_downregulated" and molecule not in downregulated_molecules)
            ): continue

            if molecule in background_lipids.keys():
                input_molecule = "Input lipid"
            elif molecule in background_proteins:
                input_molecule = "Input protein"
            elif molecule in background_metabolites:
                input_molecule = "Input metabolite"
            elif molecule in background_transcripts:
                input_molecule = "Input transcript"

            path_layers = [[TermType.INPUT_TERM, [input_molecule], ItemCounter(molecule)]]

            for i, term in enumerate(get_path(sessions[session_id].all_parent_nodes[molecule], target_term)):
                term_id = term if type(term) == str else term.term_id[0]

                if not term_id in terms:
                    if term_id in background_lipids and type(background_lipids[term_id]) == LipidAdduct:
                        lipid_category = background_lipids[term_id].get_lipid_string(LipidLevel.CATEGORY)
                        if path_layers[-1][0] != TermType.LIPID_SPECIES:
                            path_layers.append([TermType.LIPID_SPECIES, [lipid_category], ItemCounter(term_id)])
                        else:
                            path_layers[-1][2].add(term_id)

                else:
                    term_type = term.term_type
                    if term_type == TermType.LIPID_CLASS: term_type = TermType.LIPID_SPECIES
                    elif term_type in {TermType.UNREVIEWED_PROTEIN, TermType.ENSEMBLE_PROTEIN}: term_type = TermType.REVIEWED_PROTEIN
                    elif term_type == TermType.ENSEMBLE_GENE: term_type = TermType.GENE
                    elif term_type == TermType.METABOLITE: term_type = TermType.LIPID_SPECIES

                    if term_type != TermType.UNCLASSIFIED_TERM:
                        if not path_layers or path_layers[-1][0] != term_type:
                            if len(term.categories) > 0:
                                path_layers.append([term_type, term.categories, ItemCounter(term_id)])
                            else:
                                path_layers.append([term_type, [term.name], ItemCounter(term_id)])
                        else:
                            path_layers[-1][2].add(term_id)


                    else:
                        if path_layers[-1][0] != TermType.UNCLASSIFIED_TERM:
                            path_layers.append([TermType.UNCLASSIFIED_TERM, [term.name], ItemCounter(term_id)])
                        else:
                            path_layers[-1][1] = [term.name]
                            path_layers[-1][2].add(term_id)

            if len(path_layers) < 2: continue
            for i, path_layer in enumerate(path_layers):
                layer_term, layer_categories, term_ids = path_layer
                for layer_category in layer_categories:
                    if (sankey_key := (layer_term, layer_category)) not in sankey_data:
                        sankey_data[sankey_key] = {}
                        item_counter = ItemCounter()
                        item_counter.merge(term_ids)
                        sankey_data_ids[sankey_key] = [len(sankey_data_ids), item_counter]
                    else:
                        sankey_data_ids[sankey_key][1].merge(term_ids)

                    if i == 0: continue
                    layer_term_prev = path_layers[i - 1][0]
                    for layer_category_prev in path_layers[i - 1][1]:
                        sankey_key_prev = (layer_term_prev, layer_category_prev)
                        if sankey_key not in sankey_data[sankey_key_prev]:
                            sankey_data[sankey_key_prev][sankey_key] = 0
                        else: sankey_data[sankey_key_prev][sankey_key] += 1


    fig_labels = []
    fig_source = []
    fig_target = []
    fig_value =  []
    session.sankey_data = []
    for source, targets in sankey_data.items():
        (source_term, source_category) = source
        fig_labels.append(source_category)
        session.sankey_data.append(sankey_data_ids[source][1])
        source_id = sankey_data_ids[source][0]
        for target, count in targets.items():
            fig_source.append(source_id)
            fig_target.append(sankey_data_ids[target][0])
            fig_value.append(count)

    n_nodes = len(set(fig_source) | set(fig_target))
    layers = compute_layers(n_nodes, fig_source, fig_target)
    max_layer = max(layers.values())
    # Normalize x positions between 0 and 1
    node_x = [layers[i] / max_layer for i in range(n_nodes)]
    node_colors = [pastel_colors[i % len(pastel_colors)] for i in range(n_nodes)]
    edge_colors = [pastel_colors[s % len(pastel_colors)] for s in fig_source]

    fig = go.Figure(
        go.Sankey(
            node = dict(
                label = [(ontology.categories[c] if type(c) == int else c) for c in fig_labels],
                x = node_x,
                color = node_colors,
                pad = 20,
                thickness = 40,
            ),
            link = dict(
                arrowlen = 15,
                source = fig_source,
                target = fig_target,
                value = fig_value,
                color = edge_colors,
            )
        )
    )
    fig.update_layout(
        margin = dict(t = 5, l = 5, r = 5, b = 5),
        autosize = True,
        height = None,
    )
    barplot_terms_wrapper_style["height"] = "80vh"

    return (
        True,
        fig,
        barplot_terms_wrapper_style,
        no_update,
        no_update,
        [],
        True,
        "sankey_all" if ctx.triggered_id != "radiogroup_sankey" else no_update,
    )



@app.callback(
    Output("sankey_entries", "rowData"),
    Output("icon_download_sankey_entries", "disabled", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Input("sankey_graph", "clickData"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def sankey_node_clicked(clickData, session_id):
    if session_id == None or not clickData:
        raise exceptions.PreventUpdate

    if session_id not in sessions:
        return (
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
        )
    session = sessions[session_id]

    point = clickData["points"][0]
    entries_list = []
    if "pointNumber" in point:
        if "group" not in point:
            raise exceptions.PreventUpdate

        idx = point["pointNumber"]
        if session.sankey_data != None and idx < len(session.sankey_data):
            entries_list = sorted([
                {
                    "termid": k,
                    "count": v,
                } for k, v in session.sankey_data[idx].counter.items()
            ], key = lambda x: x["termid"])

    return entries_list, False, no_update, no_update



@callback(
    Output("barplot_terms_modal", "opened", allow_duplicate = True),
    Output("barplot_terms", "figure", allow_duplicate = True),
    Output("barplot_terms_modal", "size", allow_duplicate = True),
    Output("barplot_terms_wrapper", "style", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("barplot_controls", "style", allow_duplicate = True),
    Output("sunburst_controls", "style", allow_duplicate = True),
    Input("chart_results", "n_clicks"),
    Input("barplot_numberinput_connect_ths", "value"),
    Input("barplot_numberinput_font_size", "value"),
    Input("barplot_select_name", "value"),
    Input("barplot_select_sorting", "value"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    State("sessionid", "children"),
    State("barplot_controls", "style"),
    State("sunburst_controls", "style"),
    State("barplot_terms_wrapper", "style"),
    prevent_initial_call = True,
)
def open_barplot(
    n_clicks,
    jaccard_ths,
    font_size,
    bar_label,
    bar_sorting,
    row_data,
    selected_rows,
    session_id,
    barplot_controls_style,
    sunburst_controls_style,
    barplot_terms_wrapper_style,
):
    if session_id == None or jaccard_ths == None or n_clicks == None or font_size == None or bar_label not in {"id", "name"} or bar_sorting not in {BAR_SORTING_PVALUE, BAR_SORTING_SIMILARITY} or type(jaccard_ths) not in {float, int}:
        raise exceptions.PreventUpdate

    if session_id not in sessions:
        return (
            no_update,
            no_update,
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            no_update,
            no_update,
        )

    session = sessions[session_id]

    background_lipids = session.background_lipids
    background_proteins = session.background_proteins
    background_metabolites = session.background_metabolites
    background_transcripts = session.background_transcripts

    with_lipids = background_lipids != None
    with_proteins = background_proteins != None
    with_metabolites = background_metabolites != None
    with_transcripts = background_transcripts != None

    multiomics = sum([with_lipids + with_proteins + with_metabolites + with_transcripts]) > 1
    barplot_controls_style["display"] = "block"
    sunburst_controls_style["display"] = "none"

    fig = go.Figure()
    id_position = {row["termid"]: i for i, row in enumerate(row_data)}
    session_data = session.data
    selected_term_ids = np.array([t.split("|")[0] for t in sorted([row["termid"] for row in selected_rows], key = lambda x: id_position[x])])
    n = len(selected_term_ids)

    source_terms = [
        set(session_data[term_id].source_terms) for term_id in selected_term_ids
    ]
    jaccard_values = np.ones((n , n))
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            i_terms = source_terms[i]
            j_terms = source_terms[j]
            jaccard_values[j, i] = jaccard_values[i, j] = len(i_terms & j_terms) / len(i_terms | j_terms)

    if bar_sorting == BAR_SORTING_SIMILARITY and n > 1:
        Z = linkage(jaccard_values**8, method = 'average', metric = "cosine")
        sort_order = leaves_list(Z)
        selected_term_ids = selected_term_ids[sort_order]
        jaccard_values = jaccard_values[sort_order]
        jaccard_values = jaccard_values[:, sort_order]


    if session.separate_updown_switch:
        pvalues = -np.log10([session_data[term_id].pvalue_corrected[0] for term_id in selected_term_ids])
        number_regulated_entities = np.array([session_data[term_id].fisher_data[0][0] + session_data[term_id].fisher_data[1][0] for term_id in selected_term_ids])
    else:
        pvalues = -np.log10([session_data[term_id].pvalue_corrected for term_id in selected_term_ids])
        number_regulated_entities = np.array([session_data[term_id].fisher_data[0] for term_id in selected_term_ids])
    number_entities = np.array([len(session_data[term_id].source_terms) for term_id in selected_term_ids])

    domains = ["|".join(session_data[term_id].term.domain) for term_id in selected_term_ids]
    term_names = [session_data[term_id].term.name for term_id in selected_term_ids]
    jaccard_ths /= 100

    custom_data = [[s, row["term"], row["pvalue"], r, n, d] for s, row, r, n, d in zip(selected_term_ids, selected_rows, number_regulated_entities, number_entities, domains)]

    angle_gap = 360 / n
    angles = [i * angle_gap for i in range(n)]
    bar_width = angle_gap - 1
    arc_outer_radius = 100

    # Arc settings
    domains_arc_outer_radius = arc_outer_radius * 1
    domains_arc_inner_radius = arc_outer_radius * 0.84
    description_arc_outer_radius = arc_outer_radius * 0.8
    description_arc_inner_radius = arc_outer_radius * 0.74
    entities_arc_outer_radius = arc_outer_radius * 0.7
    entities_arc_inner_radius = arc_outer_radius * 0.64
    molecules_arc_inner_radius = arc_outer_radius * 0.6
    molecules_arc_outer_radius = arc_outer_radius * 0.54
    categories_arc_outer_radious = arc_outer_radius * (0.50 + (not multiomics) * 0.1)
    categories_arc_inner_radious = arc_outer_radius * (0.34 + (not multiomics) * 0.1)
    pvalue_arc_outer_radius = arc_outer_radius * (0.30 + (not multiomics) * 0.1)
    pvalue_arc_inner_radius = arc_outer_radius * (0.25 + (not multiomics) * 0.1)

    number_entities_sizes = number_entities + 5
    number_regulated_entities_sizes = number_regulated_entities + 5

    annotations = []
    mnes = np.max(number_entities_sizes)
    number_entities_sizes = number_entities_sizes / mnes
    number_regulated_entities_sizes = number_regulated_entities_sizes / mnes
    fig = go.Figure()

    if with_lipids:
        background_lipids = set(background_lipids.keys())

    # Add thick arcs above bars
    annotation_label = selected_term_ids if bar_label == "id" else term_names
    ontology = session.ontology
    sorted_domain_names = sorted(list(domain_colors.keys()))

    max_domain_num = 0
    max_categories_num = 0
    remaining_domains_categories = [[0, {dom: 0 for dom in domain_colors.keys()}, {}] for i in range(n)]
    domain_bar = bar_width * 0.95 / len(sorted_domain_names)
    category_bar = bar_width * 0.95 / 5

    for i in range(n):
        result = session_data[selected_term_ids[i]]
        center_angle = angles[i]

        description_start_angle = center_angle + bar_width / 2
        description_end_angle = center_angle - bar_width / 2

        # draw grey background
        add_arc(
            fig,
            description_start_angle,
            description_end_angle,
            0,
            arc_outer_radius,
            "#f7f7f7",
        )

        # draw the term description as horizontal radial bar
        custom = custom_data[i]
        add_arc(
            fig,
            description_start_angle,
            description_end_angle,
            description_arc_inner_radius,
            description_arc_outer_radius,
            domain_colors[list(result.term.domain)[0]][0],
            hoverinfo = f"Term ID: {custom[0]}<br />Term: {custom[1]}<br />p-value: {custom[2]}<br />Domain: {custom[5]}",
        )
        arc_mid_radius = (description_arc_inner_radius + description_arc_outer_radius) / 2
        annotate_arc(fig, annotations, annotation_label[i], angles[i], bar_width, arc_mid_radius, font_size = font_size, max_distance = arc_outer_radius, shorten = True)

        if pvalues[i] > 0:
            # draw the p value as horizontal radial bar
            pvalue_start_angle = description_start_angle
            pvalue_end_angle = pvalue_start_angle - bar_width / max(pvalues) * pvalues[i]
            pvalue = result.pvalue_corrected[0] if session.separate_updown_switch else result.pvalue_corrected
            add_arc(
                fig,
                pvalue_start_angle,
                pvalue_end_angle,
                pvalue_arc_inner_radius,
                pvalue_arc_outer_radius,
                domain_colors[list(result.term.domain)[0]][1],
                hoverinfo = f"p-value: {pvalue:.6g}"
            )

        len_lipid_table = 0
        len_protein_table = 0
        len_metabolite_table = 0
        len_transcript_table = 0
        molecules = result.source_terms
        if with_lipids: len_lipid_table = len(molecules & background_lipids)
        if with_proteins: len_protein_table = len(molecules & background_proteins)
        if with_metabolites: len_metabolite_table = len(molecules & background_metabolites)
        if with_transcripts: len_transcript_table = len(molecules & background_transcripts)

        entities_start_angle = center_angle + bar_width / 2
        entities_mid_angle = entities_start_angle
        hoverinfo = f"Regulated molecules: {number_regulated_entities[i]}<br />Associated molecules: {number_entities[i]}"
        if number_regulated_entities[i] > 0:
            entities_end_angle = entities_start_angle - number_regulated_entities_sizes[i] * bar_width
            add_arc(
                fig,
                entities_start_angle,
                entities_end_angle,
                entities_arc_inner_radius,
                entities_arc_outer_radius,
                REGULATED_COLOR,
                hoverinfo,
                degree_factor = 3,
            )
            entities_mid_angle = entities_end_angle

        if number_entities[i] - number_regulated_entities[i] > 1e-60:
            entities_end_angle = entities_mid_angle - (number_entities_sizes[i] - number_regulated_entities_sizes[i]) * bar_width
            add_arc(
                fig,
                entities_mid_angle,
                entities_end_angle,
                entities_arc_inner_radius,
                entities_arc_outer_radius,
                ASSOCIATED_COLOR,
                hoverinfo,
                degree_factor = 3,
            )

        entities_arc_mid_radius = (entities_arc_inner_radius + entities_arc_outer_radius) / 2
        arc_angle = (entities_start_angle + entities_end_angle) / 2
        entities_number_text = f"{str(int(number_regulated_entities[i]))} / {str(int(number_entities[i]))}"
        annotate_arc(fig, annotations, entities_number_text, arc_angle, abs(entities_start_angle - entities_end_angle), entities_arc_mid_radius, font_size = font_size, max_distance = arc_outer_radius)

        molecules_label = ""
        molecule_normalizer = len_lipid_table + len_protein_table + len_metabolite_table + len_transcript_table
        if multiomics:
            hoverinfo = ""
            molecules_l_start_angle = center_angle + bar_width / 2
            molecules_p_start_angle = molecules_l_start_angle
            molecules_m_start_angle = molecules_l_start_angle
            molecules_t_start_angle = molecules_l_start_angle
            arcs = []
            if with_lipids:
                molecules_label = f"{len_lipid_table}"
                hoverinfo = f"Lipids: {len_lipid_table}"
                if len_lipid_table > 0:
                    molecules_l_end_angle = molecules_l_start_angle - len_lipid_table / molecule_normalizer * bar_width
                    arcs.append([
                        molecules_l_start_angle,
                        molecules_l_end_angle,
                        molecules_arc_inner_radius,
                        molecules_arc_outer_radius,
                        LIPIDS_COLOR,
                    ])
                    molecules_p_start_angle = molecules_l_end_angle
                    molecules_m_start_angle = molecules_l_end_angle
                    molecules_t_start_angle = molecules_l_end_angle

            if with_proteins:
                molecules_label += (" / " if len(molecules_label) > 0 else "") + f"{len_protein_table}"
                hoverinfo += ("<br />" if len(hoverinfo) > 0 else "") + f"Proteins: {len_protein_table}"
                if len_protein_table > 0:
                    molecules_p_end_angle = molecules_p_start_angle - len_protein_table / molecule_normalizer * bar_width
                    arcs.append([
                        molecules_p_start_angle,
                        molecules_p_end_angle,
                        molecules_arc_inner_radius,
                        molecules_arc_outer_radius,
                        PROTEINS_COLOR,
                    ])
                    molecules_m_start_angle = molecules_p_end_angle
                    molecules_t_start_angle = molecules_p_end_angle

            if with_metabolites:
                molecules_label += (" / " if len(molecules_label) > 0 else "") + f"{len_metabolite_table}"
                hoverinfo += ("<br />" if len(hoverinfo) > 0 else "") + f"Metabolites: {len_metabolite_table}"
                if len_metabolite_table > 0:
                    molecules_m_end_angle = molecules_m_start_angle - len_metabolite_table / molecule_normalizer * bar_width
                    arcs.append([
                        molecules_m_start_angle,
                        molecules_m_end_angle,
                        molecules_arc_inner_radius,
                        molecules_arc_outer_radius,
                        METABOLITES_COLOR,
                    ])
                    molecules_t_start_angle = molecules_m_end_angle

            if with_transcripts:
                molecules_label += (" / " if len(molecules_label) > 0 else "") + f"{len_transcript_table}"
                hoverinfo += ("<br />" if len(hoverinfo) > 0 else "") + f"Transcripts: {len_transcript_table}"
                if len_transcript_table > 0:
                    molecules_t_end_angle = molecules_t_start_angle - len_transcript_table / molecule_normalizer * bar_width
                    arcs.append([
                        molecules_t_start_angle,
                        molecules_t_end_angle,
                        molecules_arc_inner_radius,
                        molecules_arc_outer_radius,
                        TRANSCRIPTS_COLOR,
                    ])
            for arc in arcs: add_arc(fig, *arc, hoverinfo)

        molecule_arc_mid_radius = (molecules_arc_inner_radius + molecules_arc_outer_radius) / 2
        annotate_arc(fig, annotations, molecules_label, angles[i], bar_width, molecule_arc_mid_radius, font_size = font_size, max_distance = arc_outer_radius)

        # Add related domain bars
        len_source_terms = len(molecules)
        for foreign_term, term_input_molecules in session.search_terms.items():
            if min(len_source_terms, len(term_input_molecules)) / max(len_source_terms, len(term_input_molecules)) < jaccard_ths: continue

            for dom in foreign_term.domain:
                remaining_domains_categories[i][1][dom] += 1
                max_domain_num = max(max_domain_num, remaining_domains_categories[i][1][dom])
        remaining_domains_categories[i][0] = description_start_angle - bar_width * 0.025

        for molecule in molecules:
            for term in get_path(session.all_parent_nodes[molecule], result.term):
                if type(term) == OntologyTerm: break

            for category in term.categories:
                category = ontology.categories[category]
                if category not in remaining_domains_categories[i][2]:
                    remaining_domains_categories[i][2][category] = 0
                remaining_domains_categories[i][2][category] += 1
                max_categories_num = max(max_categories_num, remaining_domains_categories[i][2][category])

    for i in range(n):
        domains_position = remaining_domains_categories[i][0]
        for domain_name in sorted_domain_names:
            outer_radius = (domains_arc_outer_radius - domains_arc_inner_radius) / max_domain_num * remaining_domains_categories[i][1][domain_name]
            if remaining_domains_categories[i][1][domain_name] > 0:
                add_arc(
                    fig,
                    domains_position,
                    domains_position - domain_bar,
                    domains_arc_inner_radius,
                    domains_arc_inner_radius + outer_radius,
                    domain_colors[domain_name][0],
                    hoverinfo = f"{domain_name}: {remaining_domains_categories[i][1][domain_name]}",
                    degree_factor = 5,
                )
            domains_position -= domain_bar

        categories = sorted([(k, v) for k, v in remaining_domains_categories[i][2].items()], key = lambda x: x[1], reverse = True)

        category_position = remaining_domains_categories[i][0]
        for i, (category, num) in enumerate(categories[:5]):
            outer_radius = (categories_arc_outer_radious - categories_arc_inner_radious) / max_categories_num * num
            add_arc(
                fig,
                category_position,
                category_position - category_bar,
                categories_arc_inner_radious,
                categories_arc_inner_radious + outer_radius,
                f"hsl(207, 100%, {int(36 + i * 10)}%)",
                hoverinfo = f"{category}: {num}",
                degree_factor = 5,
            )
            category_position -= category_bar


    # Add white donut hole
    outer_inner_radius = 25 + (not multiomics) * 10
    add_arc(fig, 0, 360, 0, outer_inner_radius, "#ffffff")

    if bar_sorting == BAR_SORTING_SIMILARITY and n > 1:
        add_arc(fig, 0, 360, 0, outer_inner_radius * 0.05, "#000000")

        dendrogram_points = np.zeros(((n - 1) * 4, 2))
        ddata = dendrogram(Z, no_plot = True)
        icoord = ddata['icoord']  # x-coordinates for each link
        dcoord = ddata['dcoord']  # y-coordinates (heights) for each link
        i = 0
        for xs, ys in zip(icoord, dcoord):
            for x, y in zip(xs, ys):
                dendrogram_points[i, 0] = x
                dendrogram_points[i, 1] = y
                i += 1

        dendrogram_points[:, 1] = outer_inner_radius * (1 - dendrogram_points[:,1] / max(dendrogram_points[:,1]))

        min_angle, max_angle = angles[0], angles[-1]
        min_left, max_left = min(dendrogram_points[:,0]), max(dendrogram_points[:,0])

        dendrogram_points[:, 0] = (dendrogram_points[:, 0] - min_left) / (max_left - min_left) * (max_angle - min_angle) + min_angle

        def draw_arc(figure, theta_start, theta_end, radius):
            if radius <= 0: return
            theta = np.linspace(np.radians(theta_start), np.radians(theta_end), int(abs(theta_end - theta_start)))

            # Parametric arc coordinates
            x = radius * np.cos(theta)
            y = radius * np.sin(theta)

            # Add the arc as a line
            figure.add_trace(go.Scatter(
                x = x,
                y = y,
                mode = 'lines',
                line = dict(color = 'black', width = 2),
                hoverinfo = 'skip',
                showlegend = False
            ))

        def draw_line(figure, theta, radius_1, radius_2):
            if radius_1 < 0 or radius_2 < 0: return
            theta = np.radians(theta)

            # Parametric arc coordinates
            ctheta, stheta = np.cos(theta), np.sin(theta)
            x1 = radius_1 * ctheta
            y1 = radius_1 * stheta
            x2 = radius_2 * ctheta
            y2 = radius_2 * stheta

            # Add the arc as a line
            figure.add_trace(go.Scatter(
                x = [x1, x2],
                y = [y1, y2],
                mode = 'lines',
                line = dict(color = 'black', width = 2),
                hoverinfo = 'skip',
                showlegend = False
            ))

        for i in range(n - 1):
            i *= 4
            draw_line(fig, dendrogram_points[i, 0], dendrogram_points[i, 1], dendrogram_points[i + 1, 1])
            draw_arc(fig, dendrogram_points[i + 1, 0], dendrogram_points[i + 2, 0], dendrogram_points[i + 1, 1])
            draw_line(fig, dendrogram_points[i + 2, 0], dendrogram_points[i + 2, 1], dendrogram_points[i + 3, 1])

    elif bar_sorting == BAR_SORTING_PVALUE:
        jaccard_max_saturation = 90
        for i in range(0, n - 1):
            for j in range(i + 1, n):
                jaccard = jaccard_values[i, j]
                if jaccard < jaccard_ths: continue
                if jaccard_ths < 1:
                    jaccard = int(jaccard_max_saturation - jaccard_max_saturation * (jaccard - jaccard_ths) / (1 - jaccard_ths))
                else:
                    jaccard = jaccard_max_saturation

                x0 = outer_inner_radius * np.cos(np.radians(angles[i]))
                y0 = outer_inner_radius * np.sin(np.radians(angles[i]))
                x2 = outer_inner_radius * np.cos(np.radians(angles[j]))
                y2 = outer_inner_radius * np.sin(np.radians(angles[j]))
                x1 = (x0 + x2) * 0.25
                y1 = (y0 + y2) * 0.25

                P0 = np.array([x0, y0])
                P1 = np.array([x1, y1])
                P2 = np.array([x2, y2])

                # Generate Bezier curve points
                t = np.linspace(0, 1, 10)
                bezier_points = (1 - t)[:, None]**2 * P0 + \
                                2 * (1 - t)[:, None] * t[:, None] * P1 + \
                                t[:, None]**2 * P2

                fig.add_trace(go.Scatter(
                    x = bezier_points[:, 0],
                    y = bezier_points[:, 1],
                    mode = 'lines',
                    line_shape = 'spline',
                    line = dict(color = f'hsv(0,0,{jaccard})', width = 2),
                    hoverinfo = 'skip',
                    showlegend = False
                ))

    # Layout settings with visible radial axis (height of bars)
    fig.update_layout(
        annotations = annotations,
        xaxis = dict(
            domain=[0, 1],
            range = [-arc_outer_radius, arc_outer_radius],
            showgrid = False,
            zeroline = True,
            showticklabels = False,  # Hide the tick labels on the x-axis
        ),
        yaxis = dict(
            domain = [0, 1],
            range = [-arc_outer_radius, arc_outer_radius],
            showgrid = False,
            zeroline = True,
            showticklabels = False,  # Hide the tick labels on the x-axis
        ),
        showlegend = False,
        paper_bgcolor = 'white',
        plot_bgcolor = 'white',
        dragmode = False,
        width = CIRCLE_WIDTH * 2,
        height = CIRCLE_WIDTH * 2,
        margin = dict(t = 5, l = 5, r = 5, b = 5),
    )
    barplot_terms_wrapper_style["height"] = "50%"

    return (
        True,
        fig,
        CIRCLE_WIDTH * 2 + 50,
        barplot_terms_wrapper_style,
        no_update,
        no_update,
        barplot_controls_style,
        sunburst_controls_style,
    )



@callback(
    Output("barplot_terms_modal", "opened", allow_duplicate = True),
    Output("barplot_terms", "figure", allow_duplicate = True),
    Output("barplot_terms_modal", "size", allow_duplicate = True),
    Output("barplot_terms_wrapper", "style", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("barplot_controls", "style", allow_duplicate = True),
    Output("sunburst_controls", "style", allow_duplicate = True),
    Input("histogram_results", "n_clicks"),
    State("sessionid", "children"),
    State("barplot_controls", "style"),
    State("sunburst_controls", "style"),
    State("barplot_terms_wrapper", "style"),
    prevent_initial_call = True,
)
def open_histogram(
    n_clicks,
    session_id,
    barplot_controls_style,
    sunburst_controls_style,
    barplot_terms_wrapper_style,
):

    if session_id not in sessions:
        return (
            False,
            no_update,
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            no_update,
            no_update,
        )

    sessions[session_id].time = time.time()
    results = sessions[session_id].result
    barplot_controls_style["display"] = "none"
    sunburst_controls_style["display"] = "none"

    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x = [r.pvalue for r in results],
            nbinsx = 20,
            marker = dict(color = "#1f77b4"),
        )
    )
    fig.update_layout(
        xaxis_title = 'Uncorrected p-values',
        yaxis_title = 'Count',
        margin = dict(
            l = 0,
            r = 0,
            t = 0,
            b = 0
        ),
    )
    barplot_terms_wrapper_style["height"] = "80%"
    return (
        True,
        fig,
        "90%",
        barplot_terms_wrapper_style,
        no_update,
        no_update,
        barplot_controls_style,
        sunburst_controls_style,
    )


def graph_config_code(graph_type):
    return """
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
            Plotly.downloadImage(gd, {format: "png", filename: "%s"});
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
            Plotly.downloadImage(gd, {format: "svg", filename: "%s"});
        }
    });

    graph_config.modeBarButtonsToAdd.push({
        name: 'download_html',
        title: 'Download as standalone interactive html',
        icon: {
            width: 24,
            height: 24,
            path: 'M22.71 6.29a1 1 0 0 0-1.42 0L20 7.59V2a1 1 0 0 0-2 0v5.59l-1.29-1.3a1 1 0 0 0-1.42 1.42l3 3a1 1 0 0 0 .33.21a.94.94 0 0 0 .76 0a1 1 0 0 0 .33-.21l3-3a1 1 0 0 0 0-1.42M19 13a1 1 0 0 0-1 1v.38l-1.48-1.48a2.79 2.79 0 0 0-3.93 0l-.7.7l-2.48-2.48a2.85 2.85 0 0 0-3.93 0L4 12.6V7a1 1 0 0 1 1-1h8a1 1 0 0 0 0-2H5a3 3 0 0 0-3 3v12a3 3 0 0 0 3 3h12a3 3 0 0 0 3-3v-5a1 1 0 0 0-1-1M5 20a1 1 0 0 1-1-1v-3.57l2.9-2.9a.79.79 0 0 1 1.09 0l3.17 3.17l4.3 4.3Zm13-1a.9.9 0 0 1-.18.53L13.31 15l.7-.7a.77.77 0 0 1 1.1 0L18 17.21Z',
        },
        click: function(gd) {
            var gd_figure = {
                data: gd.data,
                layout: gd.layout
            };
            var html_download_trigger = document.getElementById('html_download_trigger');
            var gd_str = JSON.stringify(gd_figure);
            var _setter = Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, "value").set;
            _setter.call(html_download_trigger, gd_str);
            var ev = new Event('input', { bubbles: true });
            html_download_trigger.dispatchEvent(ev);
        }
    });
    return window.dash_clientside.no_update;
}
""" % (graph_type, graph_type)


clientside_callback(
    graph_config_code("barplot"),
    Output("barplot_terms", "config"),
    Input("barplot_terms", "config"),
)


clientside_callback(
    graph_config_code("sankey"),
    Output("sankey_graph", "config"),
    Input("sankey_graph", "config"),
)


clientside_callback(
    """
function (recorded_action) {
    var url = "https://lifs-tools.org/matomo/matomo.php?idsite=17&rec=1&e_c=MAPtoGO-%s&e_a=" + recorded_action;
    console.log(url);
    var xhr = new XMLHttpRequest();
    xhr.open('GET', url, true);
    xhr.onreadystatechange = function () {
        if (xhr.readyState === 4) {
            if (xhr.status === 200) {
                console.log("ok");
            }
        }
    };
    xhr.send();
    return "";
}
    """ % ini_config["git"]["version"],
    Output("analytics", "children", allow_duplicate = True),
    Input("analytics", "children"),
    prevent_initial_call = True,
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
    Input("button_load_reviewed_proteins", "n_clicks"),
    Input("button_load_unreviewed_proteins", "n_clicks"),
    Input("button_load_all_proteins", "n_clicks"),
    State("select_organism", "value"),
    prevent_initial_call = True
)
def select_predefined_modal(
    reviewed_clicks,
    unreviewed_clicks,
    all_clicks,
    selected_organism,
):
    if ctx.triggered_id not in {
        "button_load_reviewed_proteins",
        "button_load_unreviewed_proteins",
        "button_load_all_proteins",
    } or selected_organism not in enrichment_ontologies:
        raise exceptions.PreventUpdate

    ontology = enrichment_ontologies[selected_organism]

    if ctx.triggered_id == "button_load_reviewed_proteins":
        protein_set = ontology.reviewed_proteins
    elif ctx.triggered_id == "button_load_unreviewed_proteins":
        protein_set = set(ontology.proteins.keys()) - ontology.reviewed_proteins
    else:
        protein_set = ontology.proteins.keys()

    return "\n".join([protein.replace("UNIPROT:", "") for protein in protein_set])



@callback(
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("term_path_area", "children", allow_duplicate = True),
    Input("term_lipids_modal_grid", "cellRendererData"),
    Input("term_proteins_modal_grid", "cellRendererData"),
    Input("term_metabolites_modal_grid", "cellRendererData"),
    Input("term_transcripts_modal_grid", "cellRendererData"),
    State("sessionid", "children"),
    State("term_molecules_modal_id", "children"),
    State("term_lipids_modal_grid", "rowData"),
    State("term_proteins_modal_grid", "rowData"),
    State("term_metabolites_modal_grid", "rowData"),
    State("term_transcripts_modal_grid", "rowData"),
    State("select_organism", "value"),
    prevent_initial_call = True,
)
def show_molecule_term_path(
    renderer_data_lipids,
    renderer_data_proteins,
    renderer_data_metabolites,
    renderer_data_transcripts,
    session_id,
    target_term_id,
    row_data_lipids,
    row_data_proteins,
    row_data_metabolites,
    row_data_transcripts,
    organism,
):
    if session_id not in sessions:
        return True, "Your session has expired. Please refresh the website.", ""

    ontology = enrichment_ontologies[organism]
    if target_term_id not in ontology.ontology_terms:
        return True, "Your session has expired. Please refresh the website.", ""

    target_term = ontology.ontology_terms[target_term_id]
    if target_term not in sessions[session_id].search_terms:
        return True, "Your session has expired. Please refresh the website.", ""

    trigger = callback_context.triggered[0]["prop_id"].split(".")[0]
    if trigger == "term_lipids_modal_grid":
        row_index = int(renderer_data_lipids["rowId"])
        if len(row_data_lipids) <= row_index:
            return True, "Your session has expired. Please refresh the website.", ""
        molecule = row_data_lipids[row_index]["molecule_id"]

    elif trigger == "term_proteins_modal_grid":
        row_index = int(renderer_data_proteins["rowId"])
        if len(row_data_proteins) <= row_index:
            return True, "Your session has expired. Please refresh the website.", ""
        molecule = row_data_proteins[row_index]["molecule_id"]

    elif trigger == "term_metabolites_modal_grid":
        row_index = int(renderer_data_metabolites["rowId"])
        if len(row_data_metabolites) <= row_index:
            return True, "Your session has expired. Please refresh the website.", ""
        molecule = row_data_metabolites[row_index]["molecule_id"]

    elif trigger == "term_transcripts_modal_grid":
        row_index = int(renderer_data_transcripts["rowId"])
        if len(row_data_transcripts) <= row_index:
            return True, "Your session has expired. Please refresh the website.", ""
        molecule = row_data_transcripts[row_index]["molecule_id"]

    term_path = []
    for i, term in enumerate(get_path(sessions[session_id].all_parent_nodes[molecule], target_term)):
        if type(term) == str: term_id = term
        else: term_id = list(term.term_id)[0]
        if i > 0: term_path.append(dmc.Text("▼", style = {"textAlign": "center"}))

        href = get_term_link(term_id)
        if term_id in ontology.ontology_terms:
            term_name = term.name if type(term) == OntologyTerm else ontology.ontology_terms[term_id].name
        else:
            term_name = term_id

        term_path.append(
            dmc.Paper(
                html.A(
                    dmc.Text(term_name),
                    href = href,
                    target = "_blank",
                    style = {
                        "color": LINK_COLOR,
                    },
                ) if href != "." else dmc.Text(term_name),
                shadow = "xs",
                style = {
                    "padding": "10px",
                    "width": "100%",
                    "marginBottom": "6px",
                    "marginTop": "6px",
                    "textAlign": "center",
                },
            ),
        )
    return no_update, no_update, term_path



ns = api.namespace("api", description = "Operations")

# Model for GO analysis
enrichment_model = api.model("Enrichment", {
    "background_lipids": fields.List(
        fields.String,
        description = "All lipid names in experiment (background, required)",
        example = ["DG 30:0", "DG 32:0", "DG 32:1", "DG 34:0", "DG 34:1", "DG 34:2", "DG 36:1", "DG 36:2", "DG 36:4", "DG 38:5", "DG 43:6", "DG 46:1", "DG 48:1", "DG 50:1", "LPA 16:0", "LPA 18:3", "LPA 18:0", "LPC 14:0", "LPC 15:0", "LPC 16:0", "LPC 16:1", "LPC 17:0", "LPC 18:0", "LPC 18:1", "LPC 18:2", "LPC 20:0", "LPC 20:1", "LPC 20:3", "LPC 20:4", "LPE 16:0", "LPE 16:1", "LPE 17:0", "LPE 18:0", "LPE 18:1", "LPE 20:4", "LPE 20:5", "LPE 22:5", "LPE 22:6", "PA 30:0", "PA 32:0", "PA 32:1", "PA 34:0", "PA 34:1", "PA 34:2", "PA 35:2", "PA 36:2", "PA 38:4", "PC O-32:1", "PC O-34:0", "PC O-34:1", "PC O-36:3", "PC O-38:3", "PC O-38:4", "PC 32:0", "PC 32:1", "PC 32:2", "PC 33:1", "PC 33:2", "PC 34:1", "PC 34:2", "PC 34:3", "PC 35:1", "PC 35:2", "PC 36:1", "PC 36:2", "PC 36:3", "PC 36:4", "PC 38:4", "PC 38:5", "PC 40:5", "PC 40:6", "PE O-32:2", "PE O-34:1", "PE O-34:2", "PE O-34:3", "PE O-35:2", "PE O-35:4", "PE O-36:1", "PE O-36:2", "PE O-36:3", "PE O-36:4", "PE O-36:5", "PE O-36:6", "PE O-37:4", "PE O-37:5", "PE O-37:6", "PE O-38:1", "PE O-38:2", "PE O-38:3", "PE O-38:4", "PE O-38:5", "PE O-38:7", "PE O-40:2", "PE O-40:3", "PE O-40:5", "PE 32:1", "PE 32:2", "PE 33:1", "PE 34:1", "PE 34:2", "PE 35:1", "PE 35:2", "PE 36:1", "PE 36:2", "PE 36:3", "PE 36:4", "PE 37:2", "PE 37:3", "PE 37:4", "PE 38:2", "PE 38:3", "PE 38:4", "PE 38:5", "PE 38:6", "PE 38:7", "PE 39:5", "PE 39:6", "PE 40:1", "PE 40:2", "PE 40:4", "PE 40:5", "PE 40:6", "PE 40:7", "PE 42:7", "PG 32:0", "PG 32:1", "PG 34:1", "PG 34:2", "PG 34:3", "PG 34:4", "PG 36:1", "PG 36:3", "PG 36:4", "PG 38:1", "PG 38:2", "PG 38:3", "PG 38:4", "PG 38:5", "PG 38:6", "PG 40:3", "PG 40:4", "PG 40:5", "PG 40:6", "PG 40:7", "PG 42:7", "PI 32:1", "PI 34:1", "PI 34:2", "PI 36:2", "PI 36:3", "PI 36:4", "PI 36:5", "PI 37:2", "PI 37:3", "PI 37:4", "PI 38:3", "PI 38:4", "PI 38:5", "PI 38:6", "PI 39:5", "PI 40:5", "PI 40:6", "PI 40:7", "PS 32:1", "PS 34:1", "PS 34:2", "PS 36:1", "PS 36:2", "PS 36:3", "PS 36:4", "PS 37:1", "PS 38:1", "PS 38:2", "PS 38:3", "PS 38:4", "PS 38:5", "PS 39:3", "PS 40:1", "PS 40:2", "PS 40:3", "PS 40:4", "PS 40:5", "PS 40:6", "PS 40:7", "Cer 34:1", "Cer 34:2", "Cer 42:1", "Cer 42:2", "SM 34:0", "SM 34:1", "SM 35:1", "SM 32:1", "SM 33:1", "SM 36:1", "SM 38:1", "SM 38:2", "SM 39:1", "SM 40:1", "SM 40:2", "SM 41:1", "SM 42:1", "SM 42:2", "TAG 44:1", "TAG 46:2", "TAG 47:2", "TAG 48:1", "TAG 48:2", "TAG 48:3", "TAG 49:2", "TAG 49:3", "TAG 50:2", "TAG 50:3", "TAG 50:4", "TAG 51:2", "TAG 51:3", "TAG 52:3", "TAG 52:4", "TAG 52:5", "TAG 53:3", "TAG 53:4", "TAG 54:3", "TAG 54:4", "TAG 54:5", "TAG 54:6", "TAG 54:7", "TAG 56:2", "TAG 56:3", "TAG 56:4", "TAG 56:5", "TAG 56:6", "TAG 56:7", "TAG 56:8", "TAG 58:2", "TAG 58:3", "TAG 58:4", "TAG 58:5", "TAG 58:7", "TAG 58:8", "TAG 58:9", "TAG 60:3", "TAG 60:4"],

    ),
    "regulated_lipids": fields.List(
        fields.String,
        description = "All regulated lipid names in experiment",
        example = ["DG 30:0", "DG 32:0", "DG 32:1", "DG 34:0", "DG 46:1", "DG 48:1", "LPA 16:0", "LPC 14:0", "LPC 15:0", "LPC 18:2", "LPC 20:1", "LPE 16:0", "LPE 16:1", "LPE 17:0", "LPE 18:0", "LPE 20:5", "LPE 22:5", "LPE 22:6", "PA 30:0", "PA 32:0", "PA 32:1", "PA 34:0", "PA 34:1", "PA 34:2", "PA 35:2", "PA 38:4", "PC O-34:1", "PC O-38:3", "PC O-38:4", "PC 32:0", "PC 32:1", "PC 33:1", "PC 33:2", "PC 35:1", "PC 40:6", "PE O-32:2", "PE O-34:1", "PE O-35:2", "PE O-35:4", "PE O-36:1", "PE O-36:2", "PE O-36:6", "PE O-37:6", "PE O-38:1", "PE O-38:2", "PE O-38:3", "PE O-40:2", "PE O-40:3", "PE O-40:5", "PE 32:1", "PE 32:2", "PE 33:1", "PE 34:1", "PE 34:2", "PE 36:4", "PE 37:2", "PE 38:3", "PE 38:7", "PE 39:5", "PE 40:1", "PE 40:2", "PE 42:7", "PG 32:0", "PG 32:1", "PG 34:3", "PG 34:4", "PG 36:3", "PG 36:4", "PG 38:1", "PG 38:2", "PG 38:3", "PG 38:4", "PG 38:5", "PG 40:3", "PG 40:4", "PG 40:5", "PG 40:6", "PG 40:7", "PG 42:7", "PI 32:1", "PI 34:1", "PI 34:2", "PI 36:2", "PI 36:4", "PI 36:5", "PI 38:3", "PI 38:6", "PI 39:5", "PI 40:7", "PS 32:1", "PS 34:1", "PS 34:2", "PS 36:3", "PS 36:4", "PS 37:1", "PS 40:5", "PS 40:6", "SM 34:0", "TAG 44:1", "TAG 46:2", "TAG 47:2", "TAG 48:1", "TAG 48:2", "TAG 49:2", "TAG 49:3", "TAG 50:2", "TAG 50:3", "TAG 51:2", "TAG 51:3", "TAG 52:3", "TAG 52:5", "TAG 53:3", "TAG 53:4", "TAG 54:7", "TAG 56:2", "TAG 56:3", "TAG 56:4", "TAG 56:5", "TAG 56:8", "TAG 58:2", "TAG 58:5", "TAG 58:9", "TAG 60:3", "TAG 60:4"],
    ),
    "background_proteins": fields.List(
        fields.String,
        description = "All protein accessions in experiment (background)",
        example = ["P04117", "Q05816", "Q99LJ8", "O35206", "Q9Z2A7", "Q9CZU6", "Q7TN99", "P50544", "P41216", "Q8K1N1", "Q05920", "Q8R1X6", "Q80UU9", "Q62351", "Q9D6R2", "Q9QXS6", "Q62314", "P63038", "Q04857", "Q9QZA0", "Q8QZT1", "P59017", "Q99L13", "P56395", "Q9D9V3", "Q62188", "P28301", "P51881", "P32020", "O08756", "P16546", "Q9WV92", "P37040", "O35855", "P54310", "Q8R050", "Q62219", "Q8CAQ8", "Q9D281", "P13707", "Q5SSZ5", "P51660", "P97315", "P24270", "P42125", "P09103", "Q9CQN1", "P53395", "Q5U4D8", "Q921H8", "Q8K2B3", "Q5SW19", "Q61646", "Q9R0H0", "O54724", "Q8VDN2", "P37804", "Q8BFR5", "Q8C0N2", "Q9JKR6", "P61022", "P49717", "Q9DBW0", "Q8BMF4", "Q69ZN7", "P11087", "Q60953", "O08749", "Q8C1E7", "Q9DBL9", "Q9DB15", "P97449", "P42227", "Q9CQ62", "P31324", "Q64429", "Q99JR1", "Q8C129", "Q7TPR4", "Q9CRD2", "Q02788", "Q8BGS7", "Q9D819", "Q08857", "Q9R0E1", "Q99KQ4", "Q07797", "Q9CXW2", "Q8BMS1", "Q9D517", "Q02248", "Q9D024", "Q60994", "Q921G7", "Q9D051", "Q8CGN5", "Q9JLM8", "Q9DAW9", "E9Q634", "Q63918", "P56480", "Q99J47", "Q9DBN5", "P47757", "P63323", "P51912", "Q8K1N2", "P06801", "O35945", "Q9QUI0", "Q61391", "Q99JY0", "Q99L88", "Q8CGU1", "Q8BHS6", "Q923Z0", "Q8BWT1", "Q9JJE7", "Q8BY87", "Q9CXW4", "Q07417", "Q9CQ19", "O55239", "P53690", "P20108", "Q62523", "Q9JLJ2", "Q8C8U0", "Q8BVA4", "Q9R1Z8", "P28271", "Q9DB76", "P08122", "Q9EP72", "P11152", "Q91YH5", "P47934", "Q8K5B2", "Q5XG73", "P12382", "P62737", "Q8K3K7", "Q71LX4", "Q9CPQ8", "Q91W92", "D0QMC3", "P08249", "Q9Z0L0", "E9Q555", "P49817", "P63268", "P51174", "Q9D2R0", "Q9WTP6", "Q03265", "Q8K0C4", "Q9WVB0", "Q91VD9", "B2RXS4", "Q62137", "A2AAY5", "Q8BWW4", "O08573", "P97372", "O35639", "Q8K1Z0", "O08715", "P13595", "Q9QXS1", "Q61739", "Q8BX70", "Q91V12", "P35922", "Q80SW1", "Q6PDI5", "Q68FL4", "Q9QZD8", "Q58NB6", "Q9CZX7", "Q8BQS5", "P97807", "P52825", "Q9CZU4", "Q64282", "Q52KR3", "Q9CXE7", "Q9QYR9", "P18872", "P97314", "Q01149", "Q8BML1", "O08553", "Q9D880", "P63017", "Q61009", "P52196", "Q9DB20", "Q7TT50", "Q8C1B7", "Q9WVK4", "O08528", "Q920E5", "Q9D032", "Q76LS9", "Q8CC88", "Q9JHS4", "Q69ZK0", "Q91V64", "P02463", "Q8BZA9", "Q91VS7", "Q9DBJ1", "Q9QXZ0", "Q64516", "Q8CAY6", "O35405", "P54116", "Q8C838", "Q91WK5", "P99029", "Q8CI51", "Q8VDP4", "Q8R349", "O88492", "O70378", "Q9ET54", "Q62261", "P17809", "P41778", "P82347", "P49718", "P63328", "Q922J9", "Q62426", "Q9CR98", "Q3UJD6", "P61982", "Q9CZW5", "Q80UM7", "Q64433", "P11276", "O70503", "P45952", "Q8R550", "Q9CQ85", "P02469", "Q8BKZ9", "Q99KI3", "P68404", "Q6ZQ73", "Q93092", "Q8CGK3", "Q9DBC7", "P15261", "Q8JZK9", "Q8K0D5", "P38647", "Q2PZL6", "Q9QUN9", "Q71FD7", "Q8BH64", "P31786", "Q9EP71", "Q9Z2Z6", "P62908", "Q9R0Y5", "P97310", "Q5SWU9", "Q9JKB3", "P26645", "Q6NS46", "P26231", "Q8R059", "Q9CZ13", "Q8QZS1", "Q8K2Z4", "Q61881", "Q91YT0", "P19096", "Q99J87", "P97300", "Q91ZR2", "Q9Z1Z2", "Q91YM4", "P08121", "P48774", "Q9QZE5", "Q9Z0X1", "Q9CQC9", "Q3UMY5", "P97311", "P62075", "P10649", "P67778", "Q91YP2", "Q9CQ65", "P18242", "Q9CQF0", "P15626", "P09055", "P50136", "P50172", "Q60605", "Q9WV96", "P40142", "P40124", "Q61335", "P21981", "Q9D0A3", "Q99JY9", "Q9WVC3", "Q501J6", "Q9D7N9", "Q8BH95", "Q921Q7", "Q80X90", "Q8BX02", "Q64511", "Q9JIK5", "Q60710", "Q8VCN5", "Q99KI0", "Q99MQ1", "Q8CI94", "Q6PB66", "O88844", "Q62417", "Q6PIE5", "Q9Z2I8", "Q99J56", "O55143", "Q9CZE3", "Q5GH64", "Q5SNZ0", "Q9CXW3", "Q9QZU9", "Q99M71", "Q99LP6", "Q9QUJ7", "Q8VDD5", "Q60766", "Q9DBE0", "P46735", "Q9CQV8", "Q8R2Q4", "Q9CQX2", "Q8C7V3", "P01027", "P47856", "Q9DCW4", "Q8VEE4", "O70423", "P62331", "Q9D365", "P97465", "P60670", "Q99LC9", "A6H5Z3", "Q8R1A4", "Q8BGI4", "Q5F2E8", "Q9WVA4", "Q3TBT3", "Q60634", "O35129", "Q8VE99", "Q9R008", "Q9ER00", "P05622", "P00493", "F8VPU2", "Q9Z191", "Q3UIR3", "Q9Z0P4", "Q8R180", "Q8BTM8", "P70336", "Q9DBS9", "Q9QZH6", "Q9Z2V4", "Q64449", "Q61699", "Q3UGP9", "P98078", "O70583", "Q62130", "P50518", "P51150", "Q8JZQ2", "P54775", "Q91V41", "O35153", "P70404", "Q8BGH2", "Q8C6E0", "Q8C052", "Q8BKG3", "P25085", "P18654", "O08579", "Q9DBR7", "Q8JZV7", "P70168", "Q9WTP2", "Q99JW4", "Q9D1P0", "Q9DBX1", "Q9QXE0", "Q8BJ56", "Q8CG48", "Q6ZPJ0", "Q9QXG4", "P09671", "Q8CC35", "Q9R0E2", "Q60875", "Q9DBG5", "Q62356", "P62889", "Q9R1P0", "Q3UFY8", "Q8JZM0", "Q8K1E0", "Q6DID7", "Q9QVP9", "Q9Z1Y4", "Q8BMS9", "P55096", "Q8K124", "Q99MN9", "Q04899", "Q91YJ2", "Q9D7P6", "Q921N6", "P35486", "Q8BVI4", "Q60847", "Q9ES97", "Q9D0M5", "P81122", "P31428", "Q9CRD0", "Q9ER88", "Q8K2A1", "Q9CYV5", "Q8C3W1", "Q3TA38", "Q3U962", "P70677", "Q61425", "Q99LD8", "P08003", "Q8CJ26", "Q8BJS4", "Q9EPC1", "P11440", "P30999", "P27048", "Q3TDQ1", "Q8VCM8", "Q9D0K1", "P26041", "P28653", "P08113", "Q60767", "P17918", "P35700", "Q8K4Q8", "O09131", "Q9DAR7", "Q99LC5", "P35293", "Q99LJ0", "P14152", "Q9D1N9", "A2AJK6", "Q3U7R1", "P39749", "P62281", "Q9ERG0", "Q9QWY8", "P25206", "Q6ZQL4", "Q8R3V5", "Q01279", "O88746", "O55131", "Q91YX5", "P46638", "Q91VR2", "Q60992", "Q9CR62", "Q99KR7", "P13020", "Q9QZS8", "Q921T2", "Q9CZW4", "Q9CRB5", "Q9D8Z2", "A6X8Z5", "Q9WVA2", "Q8K297", "Q922R8", "A2AF47", "P47968", "O89110", "Q03173", "Q8R5K4", "O70401", "Q62425", "P55194", "Q8BHB4", "Q8CFI0", "Q8VDL4", "Q6P5E4", "P54227", "Q9D1A2", "P14901", "Q61553", "Q9D6Y7", "Q9Z1E4", "P20444", "P18828", "Q3U0B3", "Q3TJD7", "P58137", "Q8K2C7", "P10107", "Q3THG9", "Q9D6Z1", "Q9CQI3", "P98192", "P16332", "Q9ERT9", "Q9CX86", "Q61703", "O08547", "Q8C6B9", "P97864", "Q9D5V6", "Q8BG32", "P39098", "P58281", "Q9WV55", "Q9D154", "Q9JK81", "P70662", "Q01730", "Q8R404", "Q9DC69", "Q9DCD0", "O35379", "Q8BG51", "Q8R0N6", "P97379", "Q2TPA8", "Q8C156", "Q9WVL3", "Q8R0X7", "Q3TW96", "P05064", "Q60823", "Q8CCF0", "Q4LDD4", "Q62186", "Q9D1E6", "P18155", "P52432", "Q80Y17", "Q99K01", "Q6P3A8", "Q9QX11", "Q61107", "Q9R1J0", "Q8VHX6", "Q9D0G0", "Q925I1", "Q9CZ52", "Q99K70", "O70566", "Q8CCI5", "Q8C0J6", "Q04207", "Q9DCX2", "Q91ZA3", "Q9Z1Q2", "O55222", "Q3TN34", "Q6PHZ2", "P23492", "O35099", "Q91VU0", "Q3TDD9", "O35343", "Q8BVF2", "P00405", "P56135", "Q3V1G4", "Q8C079", "P19426", "Q8R1G6", "Q9D1L0", "Q9JLC4", "Q3UJB9", "Q9QZ23", "Q8K1I7", "Q8BPB5", "Q9D1L9", "P58771", "P23198", "Q8BWS5", "Q9R062", "P33434", "P48962", "Q99K51", "P70290", "Q91V01", "Q9JLN9", "Q78IK4", "Q71RI9", "Q9Z0H3", "Q8BUV3", "O35083", "P63037", "Q6ZWQ0", "P23475", "Q9WVL2", "Q8K479", "P56391", "Q9Z277", "Q99NH0", "O35874", "Q8C0L9", "P16460", "Q9Z0V7", "O35955", "O08807", "Q80WQ2", "Q8C6I2", "P42703", "Q78HU7", "P40224", "Q91VA6", "Q61245", "P56389", "Q99LX0", "Q9JKZ2", "P62320", "Q9CQQ7", "Q9CZR8", "O70252", "Q8CFX1", "Q91XY4", "Q80V53", "Q8BHJ5", "O88967", "Q08619", "Q3UVK0", "Q99MD9", "Q8C7X2", "O08795", "Q69Z37", "O88876", "Q9Z1Q9", "Q8C522", "E9Q9A9", "Q9WTP7", "Q9WUM4", "P70444", "Q9CWJ9", "A2A5R2", "Q60972", "Q8CBY0", "O35730", "P21126", "Q9Z2U1", "Q9Z1X4", "Q9DBH5", "O55023", "Q9Z2G6", "P62821", "Q8VCF0", "Q9CY64", "P08228", "P19783", "P16054", "O35215", "Q8VCI5", "Q920A7", "Q9D1X0", "Q9R0B9", "Q8BFW7", "Q6GQT9", "Q61490", "Q60930", "P10605", "P50580", "Q9CY50", "P61924", "Q8R2U6", "P15532", "Q7TNP2", "Q6R891", "P16110", "Q8BVE3", "Q8VE62", "Q8VHL1", "Q9QXL2", "Q91YR1", "Q99PG2", "Q9D0F9", "Q8CIE6", "Q8R317", "D3Z7P3", "Q9D023", "O70318", "Q9DCL9", "Q5SXY1", "P63087", "Q8R570", "O35927", "Q99LC8", "Q920B9", "Q9EQK5", "Q6P5B0", "Q925B0", "Q9DCJ1", "P70202", "Q9DCA2", "Q9DBE9", "Q8BIJ6", "Q8R5F7", "P53810", "Q8BQ30", "Q3UPL0", "Q9DB27", "P56394", "Q9Z0E6", "Q3U0J8", "Q8R138", "O35459", "Q99M87", "Q8VBZ3", "Q3U186", "Q14C51", "Q3TC93", "Q8BYK6", "Q9D7B7", "Q8BMD8", "Q8BGT5", "O88544", "Q8BWY9", "Q9ERS2", "Q8BP71", "Q99N89", "Q9DBL1", "Q922Q8", "Q80Y14", "Q8CFH6", "Q9QYC0", "Q60714", "Q8R3S6", "Q3TB82", "Q62087", "Q9EQ20", "Q9CQA3", "Q3TAS6", "Q3UIU2", "Q91V92", "Q91VE6", "Q8VDP3", "P27546", "Q8BJH1", "Q6P3Y5", "Q9D7G0", "P47809", "Q8BLN5", "P21550", "Q8VEE1", "Q9QZ85", "O54931", "P11928", "Q9Z0P5", "Q9CQ56", "Q99MR8", "Q6Q899", "Q91YD6", "Q8VDM6", "Q8BX65", "Q8K268", "P46664", "Q61686", "Q922F4", "Q91W50", "P97346", "O88384", "Q9DCS3", "Q60772", "Q9D898", "Q61879", "Q9QZK2", "Q9CPV2", "P11688", "Q9JKF1", "Q99LE6", "P70227", "Q8CDJ8", "Q9D1K2", "Q9EQC5", "Q9D7S7", "Q9WVL0", "P06151", "Q9JI39", "Q8R1U1", "Q9Z1Z0", "P61087", "P25799", "P45376", "P61290", "Q8VEA4", "Q99MU3", "Q8K2M0", "O88441", "Q91VJ2", "Q8R4X3", "Q9JHU4", "A2RTL5", "Q3THE2", "Q99J95", "Q925E0", "Q99JI4", "Q8BKI2", "Q8VE22", "Q8C033", "Q70E20", "Q9Z0H1", "P39053", "Q3UHJ0", "Q6PDQ2", "Q99LJ5", "Q8BWZ3", "Q5NBU8", "P24369", "Q640N1", "Q99ME9", "Q8R3B1", "B1AZP2", "Q91Z92", "Q80TA1", "Q9DBP5", "P41241", "O54890", "Q8BTT6", "Q9Z0R9", "Q60876", "Q9ER38", "P17047", "Q9CVB6", "Q9CPR5", "Q9CQV5", "Q99KV1", "Q3UBZ5", "Q8R0Y6", "Q7TSV4", "P56183", "Q8BR92", "Q00612", "Q9DBX2", "P26187", "Q2KN98", "O88379", "Q9D125", "P80316", "P48024", "P46471", "Q9D7V9", "Q0VGY8", "Q9JII6", "P54276", "Q9D0B6", "Q9EPQ7", "Q3TH73", "Q80ZJ1", "Q9QYX7", "P50171", "P83741", "Q3V1L4", "Q7TMY8", "Q9CZ30", "Q9CWK8", "Q9CY27", "Q9R1C7", "Q99P72", "Q8VDM4", "Q60960", "Q91ZW3", "P99026", "Q8CI95", "Q9R0M6", "Q91ZJ5", "P39447", "Q9D0M1", "Q9CU62", "Q9EPB5", "P46414", "Q9CQ60", "Q8BGS2", "Q62167", "Q9WVH4", "Q68FH4", "P29391", "Q8BPB0", "Q99JF5", "Q3TXS7", "Q8CCK0", "Q61398", "Q8BW96", "P15116", "P48428", "Q91WM2", "Q61941", "P14576", "P59325", "Q9Z2I9", "Q8K1A6", "Q60932", "P22682", "P05201", "P05533", "O09061", "Q61112", "Q3URD3", "O35623", "Q8C0I1", "P17141", "Q6ZWY3", "Q4VA53", "Q60817", "Q8C0C7", "Q5XG71", "P28660", "Q9ER72", "Q6NZC7", "Q9CQR2", "Q61810", "Q91VH6", "Q60838", "Q9DBL7", "Q99JP6", "Q63844", "Q9CWR0", "Q9JHJ0", "P15864", "Q9Z1F9", "Q6PHU5", "Q8K4J6", "Q9CQF9", "P99024", "O89112", "Q91VH2", "O35551", "P30681", "Q9DCR2", "Q9EPE9", "Q924Z4", "Q8R2G6", "P60824", "B7ZMP1", "Q60790", "P16330", "P33174", "Q80WJ7", "P09528", "Q9D8N0", "Q9DCN2", "P35564", "Q9DBZ5", "Q6P4T0", "Q99NB8", "Q6PAC3", "Q64345", "Q9DAK9", "Q61207", "Q3UR32", "Q8R310", "Q80TH2", "P97927", "Q8BH97", "P61161", "Q9WVA3", "Q7TMK9", "Q91V76", "P35441", "Q8R3F5", "Q9QY24", "Q64435", "Q9DBS1", "Q8C878", "Q9CPY7", "Q8K0Z7", "P15535", "Q61733", "Q8BSF4", "Q91YY4", "O35640", "P35276", "O88286", "P24288", "Q8BIW1", "Q3U0V2", "O88653", "Q9CXY9", "P10711", "Q9Z0S1", "Q921M7", "Q9JMH9", "P18572", "Q9D2G2", "Q6IRU2", "Q60739", "Q6X893", "O54754", "Q64133", "Q5I012", "O88271", "Q60902", "Q921M3", "Q9CXJ1", "Q00915", "Q9EPT5", "P70257", "Q5F2F2", "Q08093", "Q11136", "Q61235", "Q99JX4", "Q9DCH4", "P22907", "P62858", "O35857", "Q9DCV4", "Q9R1Q9", "Q5EBG8", "Q9JI13", "Q8K4F5", "Q6DFW4", "Q61598", "Q8BPU7", "P14873", "Q9Z0H4", "Q9Z247", "Q8VDJ3", "Q9JKV1", "P97360", "Q8R1K4", "Q9CWG9", "Q9R257", "P14733", "P25976", "Q08024", "Q64261", "Q8BHN3", "Q9CX97", "Q9D1J1", "O88532", "Q76MZ3", "P14602", "Q9DCF9", "Q8BNU0", "Q9Z0F8", "Q9Z2W0", "P27641", "Q569Z5", "Q8R3D1", "Q8BH79", "Q62000", "Q9CQZ1", "Q9DB77", "Q810A7", "Q9D6J6", "P97493", "Q62077", "Q8K1X1", "Q99K48", "Q99JZ7", "P24452", "P46935", "O08784", "Q9QYA2", "Q8BGD8", "Q3TBW2", "Q9D3L0", "Q6P9S7", "Q8BK72", "Q9DCC4", "P47738", "Q9CPQ3", "Q01705", "Q8BYU6", "Q9DCJ5", "B2RY56", "Q9D8C4", "Q91WC3", "Q06335", "P56873", "Q9WVG6", "Q8R1F1", "Q60631", "Q9JLI8", "Q8K370", "Q9EPJ9", "Q62418", "A6H5X4", "O08992", "Q9WUZ9", "Q8K157", "P62814", "P28063", "Q8K411", "Q68FD5", "Q8BTY1", "Q9CXY1", "P58044", "Q9Z2X1", "B1AR13", "Q9CQT5", "Q9JIM1", "P62715", "Q8C166", "Q9D3P8", "P35550", "Q9D0E1", "P10518", "P97452", "Q8K224", "Q6P6J9", "Q8C351", "O88630", "Q9CQ22", "Q8BK03", "Q9CYG7", "Q9D6S7", "Q9D1G2", "Q61510", "P97450", "P62754", "Q91X88", "Q32NZ6", "Q69ZW3", "Q9WTQ5", "P42932", "Q99K23", "Q9Z103", "P11862", "Q9Z280", "Q9WTM5", "Q9CQ91", "Q3TMX7", "Q8R081", "O55126", "Q8K298", "P62880", "Q91WQ3", "Q9CQE7", "Q3V3R1", "Q3TXU5", "Q9JIY5", "Q99L45", "P62192", "Q6PGL7", "Q80T69", "Q9Z2I0", "Q3TX08", "P62322", "Q9D358", "Q62318", "Q9WV98", "Q9ESX5", "Q3U3R4", "P40336", "Q61543", "Q8VCG3", "Q9CZL2", "O08539", "Q9D0K2", "Q6NZJ6", "Q78JW9", "O70589", "Q3UN04", "P09925", "Q3UTJ2", "Q9QYJ0", "Q811L6", "Q6PDG5", "Q9EPK6", "Q9CXY6", "P16125", "P61979", "Q3TVI8", "Q9WUH1", "Q99KU0", "Q8R035", "Q7TNG5", "Q00422", "P61588", "Q924C1", "Q9CR23", "O35309", "Q9WU56", "Q99M28", "P11352", "Q8CG47", "P51855", "Q9CQN3", "Q60854", "Q8BHK3", "Q80VQ0", "Q9DCS9", "O70433", "Q99N94", "P62897", "P53569", "Q921G6", "P09405", "P97429", "P63168", "Q9D1C9", "Q8R5J9", "P36916", "Q91ZP3", "P17710", "Q9CR16", "Q9Z1B3", "Q99NF3", "Q9D8S3", "Q8CAK1", "Q9D945", "P16045", "Q8BH48", "Q91YM2", "Q91VF6", "O88738", "Q9CQU0", "P57780", "P80318", "Q3ULJ0", "Q9Z2C4", "Q9CRC9", "Q7TMX5", "P33611", "Q9R112", "Q7TSQ8", "P70460", "P19536", "O54879", "Q9CQ92", "Q3UH60", "Q923D4", "P33610", "P05132", "O35245", "Q9WUM5", "Q8BGC0", "P58058", "Q9D0I9", "Q8BWW9", "Q9JJC6", "Q9R0P6", "Q9DCD5", "Q9JI48", "Q9Z1J3", "P60060", "Q8JZY4", "Q9EQ61", "Q9R0P5", "Q3KNM2", "Q8CHY6", "Q9Z2U0", "P17156", "Q922Q1", "Q6ZPZ3", "Q91Z53", "Q91WS0", "Q8R3U1", "P60710", "Q8BK67", "Q60648", "Q9QXY9", "Q99LH1", "Q9CX60", "P58774", "O88643", "Q6P542", "Q3TFD2", "P43274", "Q91YW3", "P46938", "Q8BWG8", "Q5EG47", "Q8VD75", "P43277", "Q9D7X3", "A2ASS6", "Q8CHH9", "P59999", "Q62241", "P17751", "P83887", "Q8BG95", "Q5FWI3", "Q9CWG8", "P10923", "O88487", "Q99PT1", "Q9ERR7", "Q99KB8", "Q64727", "Q3TMW1", "Q3U5Q7", "Q99M04", "P27612", "Q9CYA0", "O88322", "Q3ULF4", "Q9DCG2", "Q9D8U8", "Q8R361", "Q8BHG1", "Q80ZK0", "Q9QZH3", "Q61545", "Q8JZN5", "Q64339", "Q921F2", "P56399", "O70475", "P53026", "P62141", "Q8CHT3", "Q9CQJ8", "P97819", "A6H8H2", "E9Q557", "P60843", "P35951", "Q9DA69", "Q80VD1", "Q3UDR8", "Q8VCH8", "P68368", "P52912", "Q99JF8", "Q9D2C7", "Q9R1P3", "Q91YQ5", "Q6PHQ8", "Q920A5", "P84078", "P26039", "A3KGB4", "Q68FF6", "Q99MS7", "Q61550", "P50516", "Q3UMC0", "P63028", "P14142", "Q7TQI3", "P08752", "P25444", "Q6PCP5", "Q60692", "Q9WTQ8", "Q9JJ28", "Q9ESZ8", "Q06180", "P43276", "Q9CRA4", "Q9CYR6", "P03958", "Q8BJU0", "Q61081", "Q9D6Y9", "Q91ZE0", "Q9QZM4", "Q6PGH2", "Q01768", "P97384", "Q501J7", "Q80X85", "P01895", "Q61316", "P97742", "Q8C5K5", "Q64669", "Q6P5E6", "Q9D6J5", "Q9DCM0", "Q9WV84", "Q6P4S6", "Q9D0I8", "Q7TSI3", "Q61704", "A2AWA9", "Q9Z2L6", "P97371", "Q8BSY0", "Q921S7", "P53564", "Q3UH93", "P53702", "Q9QZ08", "Q8BV66", "Q8BRT1", "Q6PB93", "Q99KN9", "Q3TPE9", "Q9D9Z5", "Q9CXT8", "Q9DBR4", "Q9JJU8", "P63158", "P62245", "Q9JI10", "Q499X9", "Q9Z255", "P10493", "O35114", "Q8BL07", "Q8BNY6", "Q61656", "P62242", "Q78PY7", "P62996", "Q8C4J7", "Q9QYK7", "P36993", "P61166", "Q8R3N1", "P05213", "Q6A0A2", "Q9JM14", "O35973", "P55302", "Q3UFY7", "O54692", "Q6KCD5", "O88712", "Q922B2", "Q9JIK9", "O55013", "P10126", "P35285", "P70372", "Q9QYR6", "Q8CHP8", "P26040", "Q9DB60", "O70481", "Q6P6L0", "Q924L1", "Q812A2", "Q9CQT9", "Q921U8", "Q9JJI8", "Q9CWE0", "Q9JLZ3", "Q8C1A5", "P56959", "Q8BWM0", "P97376", "P97461", "Q6Y7W8", "P15092", "P35235", "Q9Z0E0", "Q3ULD5", "Q8BH59", "P34022", "P39689", "P62962", "Q11011", "O88342", "O08750", "Q8BG60", "O35710", "Q64521", "Q8K3C3", "Q922Q9", "Q60575", "P50285", "Q6PGB6", "Q91VN4", "Q8R151", "O35435", "Q9D6M3", "P41233", "Q8BUN5", "Q80X82", "Q6A058", "G5E8K5", "Q9D8Y0", "Q8BP48", "Q8VEB4", "Q62384", "Q9D1M0", "Q80UK0", "Q99PP9", "Q9CYL5", "Q8BFS6", "P97386", "P45878", "Q61211", "Q9D967", "P70441", "Q922W5", "Q9D379", "Q60967", "Q3UGC7", "Q9D0T1", "P20029", "Q9R0P4", "P20357", "Q9CR88", "P35278", "Q7TQ95", "Q8C0L6", "P61967", "O54825", "Q91XA2", "Q9R0Z9", "Q3TZZ7", "O70492", "Q9QZL0", "Q8K3X4", "Q922K7", "P62484", "P54822", "P07214", "P22892", "Q8C0D5", "O70310", "Q61206", "O88796", "P56380", "P55937", "Q99J99", "Q91YI0", "Q99KW3", "Q9CRB9", "Q8BGA5", "A2AGT5", "Q9DAI2", "P47754", "Q9Z2F2", "P97820", "Q99L28", "P14685", "Q8K371", "Q6ZQK5", "P14211", "P15105", "Q8R409", "Q9D8S9", "Q8VHE0", "Q8CJF9", "Q2EMV9", "Q8CGF1", "Q8BYA0", "Q99N93", "Q8BH04", "Q61576", "Q8C5Q4", "Q8VDS8", "Q9R0P3", "O35737", "Q9JHL1", "Q91Y74", "O35658", "P63276", "Q8BJZ4", "Q9CPU4", "Q8C5H8", "Q61824", "Q9ERG2", "P62259", "P63094", "Q05BC3", "Q9WV60", "Q68FL6", "Q3UDK1", "Q91VD1", "Q9JKF7", "P56546", "Q99JH1", "Q62074", "O88520", "P97377", "Q8BLY1", "Q9D1Q6", "Q4VBD2", "Q6P8X1", "P61961", "O55229", "Q8VEH6", "Q8QZR5", "Q60866", "Q9JJG9", "Q60759", "Q6PAR5", "Q61823", "P07091", "Q61768", "Q9R1Z7", "Q91YN9", "Q8BT60", "P08030", "Q8K2W3", "O88447", "Q8K0C9", "P06745", "Q9JM76", "P57746", "P80314", "Q8R105", "Q8CDN6", "P27659", "Q9D7E3", "P61025", "Q0GNC1", "Q9QXX4", "P16406", "Q9CZX9", "Q8JZZ7", "Q3UKJ7", "Q5BLK4", "Q8K0C8", "Q91YE6", "P17563", "P52795", "Q9JIZ9", "Q8BFY9", "Q9CZD3", "Q9CR68", "O70305", "Q9Z1Q5", "O35344", "Q9Z2A5", "Q99K85", "P04104", "Q9WV80", "Q64253", "P57784", "Q9DC33", "Q5RKZ7", "P27808", "P38060", "Q6ZWN5", "Q8BMS4", "O88428", "Q9D906", "Q8R311", "Q9D6N1", "P12787", "Q8BHD7", "O55028", "P49070", "P05480", "Q8BX10", "P48678", "Q9BCZ4", "Q99KC8", "Q920I9", "Q91Z67", "Q7TN98", "P70303", "Q5DTX6", "Q99LD9", "Q8CAS9", "Q9QYH6", "Q8K363", "P62073", "Q9Z0R4", "Q8CIB5", "Q9CSU0", "Q9D8B3", "Q9Z2A0", "Q9D711", "Q924K8", "Q923G2", "Q80V26", "A2BH40", "Q80TM9", "P62082", "Q01853", "P35980", "Q8CI11", "Q8BZM1", "Q69Z38", "O70139", "Q9EPA7", "Q9D8V7", "Q9DBC3", "Q9CQA9", "Q8VCA8", "Q7TMB8", "P12815", "Q8BRF7", "P23949", "Q6A0A9", "O08759", "Q64152", "P62748", "Q61792", "P0C0S6", "P27046", "P20152", "Q8BMA6", "Q8BMJ2", "P14869", "Q8BXZ1", "Q6NS82", "E9QAM5", "Q9CQX8", "Q6ZWV3", "Q8K327", "Q80ZQ9", "Q8BX90", "Q61624", "Q3U829", "Q8VEH3", "O35972", "Q8VDW0", "Q99KK2", "P35505", "Q8CIL4", "P58021", "Q8WTY4", "Q8BUY5", "P24472", "Q9QX47", "Q9WUQ2", "Q6ZQ38", "Q6KAR6", "Q9JMG7", "Q8VC30", "Q07076", "P70414", "P57716", "Q61144", "Q9CQ45", "Q9CR95", "Q8K183", "Q5SS00", "P07901", "P97351", "Q8BVQ5", "Q8R1V4", "Q921X9", "Q9DBA6", "Q6PDM2", "P50427", "Q9JJ80", "Q9ESW4", "Q9CQH3", "Q9DAW6", "Q9EPL9", "O54781", "O54916", "Q8CI61", "Q9Z1R2", "P62196", "Q8K072", "Q9JLR9", "O09106", "Q921G8", "Q640M1", "P62257", "P70340", "Q9CPN8", "P30412", "Q9JM96", "Q9EQ28", "Q3UJU9", "P58252", "Q9CWI3", "P63321", "Q8C0E2", "Q8BUR4", "Q9D136", "Q8K1H1", "Q921K9", "Q9D273", "Q8C407", "Q91WD5", "Q3TEA8", "Q6PFQ7", "Q3URE1", "Q9D0M3", "Q8VDT9", "Q8BG05", "P83940", "P31750", "P62774", "Q8BP40", "O35704", "Q9QUM9", "Q1HFZ0", "B2RY04", "P60904", "Q91VS8", "Q8C7Q4", "P63248", "Q02614", "Q80XI4", "P42128", "Q91VM3", "O35345", "Q99JB2", "Q8CCS6", "Q9JJX6", "P47713", "Q9CXI5", "Q9DC60", "Q62383", "Q8R4N0", "P14148", "Q8BZH4", "P28352", "P57759", "Q9ES28", "Q9EST5", "Q9CRA5", "Q8CIN4", "Q6PIP5", "Q9EQU5", "Q9QZD9", "Q6Y685", "Q9CWK3", "Q03145", "P97863", "Q78IK2", "Q9CQJ6", "P62311", "O08917", "Q9D824", "Q8BT07", "Q9CQZ5", "P06537", "Q9JIS8", "Q3TZM9", "Q99KR3", "Q9JIF0", "Q9D554", "Q8BHL5", "O54946", "Q9DBF1", "Q9Z1D1", "P29341", "Q9JLB0", "P11499", "Q8BVG4", "Q8CJ53", "P35979", "Q9D903", "P67984", "Q9QY76", "Q8R5A3", "O08915", "Q3TDN2", "Q9CWU9", "P62071", "Q61387", "Q6DIC0", "O08919", "Q6P9Q4", "Q9D7S9", "Q91YI1", "Q8QZY1", "O70493", "Q148V8", "O88952", "O35075", "Q9DCT5", "P47811", "Q8BW75", "Q5SWD9", "P47962", "Q8BU30", "Q9JI11", "Q9CQ13", "Q6P9J9", "P19253", "Q9CWS0", "Q9DBG9", "Q9DCL8", "Q9DC53", "P60335", "Q8VCF1", "Q6P4T2", "Q6A065", "Q8VI75", "P58871", "Q8K3K8", "P62627", "O08912", "Q9D0J4", "Q80Y98", "P47955", "Q5SSI6", "P55012", "Q921H9", "Q6IR34", "Q9JJW0", "O09012", "Q91WK2", "Q8BK64", "A2AVA0", "P80317", "Q3TWN3", "Q99LJ6", "Q3UGR5", "O70458", "Q9WUN2", "Q3UYV9", "B9EKI3", "Q8R326", "P84096", "P29533", "Q9QZC8", "Q62448", "P47199", "Q99LI7", "Q99KY4", "Q9CZG3", "Q63810", "Q9Z0N1", "Q9CYN2", "O88543", "Q91YR9", "Q8K4G5", "Q8CIF4", "P62849", "Q7TMM9", "O88291", "Q03147", "Q8R2W9", "Q91YS8", "Q99PU8", "P58059", "Q8VDP6", "P47941", "Q9JKB1", "Q8R3C6", "Q62084", "P62806", "P52479", "P70362", "Q9JHK4", "Q921K8", "A2A8Z1", "P60867", "P28481", "Q9CQE1", "Q921I9", "O88413", "Q61029", "Q7TSH3", "Q3U2P1", "Q9DC51", "Q9JI46", "P36552", "Q80XR2", "P10852", "Q99LJ1", "Q9JIX8", "P00375", "Q921C5", "Q62086", "Q3UPH7", "Q9CY18", "Q8R003", "Q8BGB5", "Q9JKP5", "Q9QXK3", "Q8CIG8", "Q8BG67", "P12023", "Q5SSW2", "Q921Z5", "P21107", "Q9D0B0", "Q9JJ89", "Q9CQ79", "Q9CYH6", "Q9CW46", "Q99JT1", "Q80XI3", "Q64010", "Q9DC61", "Q8BR63", "Q99LE1", "Q9JII5", "Q6NWV3", "Q9D8M4", "Q9Z2G9", "P30280", "P10630", "Q91W86", "Q9JK92", "Q9CZP5", "Q9CQ40", "Q9DC29", "P59808", "P35601", "Q8BYN3", "Q9EP69", "Q810L4", "P68372", "Q91WG8", "Q8BGC4", "P57080", "O89032", "Q99LF4", "Q8K2T8", "Q9JIG7", "Q99MI1", "Q6P069", "P39054", "Q61127", "Q91WC0", "D3YZP9", "P70399", "P70445", "Q9DBG7", "P26450", "Q9R1T4", "P27601", "Q6P5D8", "Q9D172", "Q5SRX1", "Q99JX3", "P02468", "P43883", "Q9R1P1", "Q62311", "P47791", "P50543", "Q61102", "P54729", "Q9CYD3", "O88207", "O54941", "P52332", "O55106", "Q60855", "Q5SSL4", "P54754", "Q9D1I6", "Q9D1I5", "Q8R3Q6", "Q9Z2R6", "Q9WV32", "P08556", "O08529", "Q9DCT2", "P70699", "Q3UPH1", "P22437", "Q8JZM7", "P32507", "Q91WM3", "Q8BZS9", "P59266", "Q8BZX4", "Q91YR7", "Q5Y5T1", "Q9EPU4", "Q80YE7", "Q8K021", "Q99J45", "Q9CQY5", "Q9D3D0", "P0CW02", "Q3TFQ1", "Q59J78", "A2AJI0", "Q99KU1", "P21619", "Q8R1J9", "Q4PJX1", "Q9WUR2", "P18826", "Q6PE15", "Q9ESU6", "Q9CR89", "Q8BP27", "P63024", "Q9CQA6", "Q8K009", "Q6VN19", "Q6NZF1", "P0C0A3", "Q8BGE6", "P17182", "Q61578", "Q9CXF4", "Q8BUH8", "Q9Z1G4", "Q8BIH0", "Q9DBZ1", "Q61097", "P59481", "Q9D706", "Q9ERU9", "Q99LM2", "Q5F2E7", "Q8VDS4", "Q9CZG9", "Q8R016", "Q922J3", "B1AQJ2", "Q3UGS4", "Q8BM85", "Q64310", "Q6P9Q6", "O54998", "Q9QZQ1", "Q9D6U8", "Q9CQU5", "Q6ZQI3", "Q91WT8", "Q9WUR9", "P52927", "Q8BY71", "Q8K4Q0", "Q3UQ84", "Q6ZQ29", "P09411", "P58801", "Q9QZ06", "Q69ZQ2", "Q9JHR7", "Q5SVQ0", "Q921I2", "O88689", "Q9QWV9", "Q9CQ06", "Q9R1T2", "Q9D1M4", "Q91YD9", "Q08501", "O35593", "Q62419", "Q7TQK5", "Q9CQF3", "Q91VI7", "Q9D4H1", "Q6PFR5", "Q00560", "Q91VM9", "P16675", "Q80UG5", "P47911", "Q8K1M6", "Q9DB96", "Q8BGR9", "Q8VEJ9", "Q64105", "Q9CQM9", "O55060", "P97770", "Q99JP4", "Q60931", "P23116", "Q9QYI3", "Q9Z0H8", "Q9D0R2", "Q6NV83", "Q791V5", "P42859", "Q9D1C8", "Q9JKX4", "Q8BKF1", "P62932", "P35821", "Q9DBR1", "E1U8D0", "Q8BVU0", "P10922", "Q8BTS4", "Q9JJT0", "Q9ERK4", "P97817", "O70496", "Q6P9S0", "Q8R0H9", "Q8CDG3", "Q05CL8", "Q8C7K6", "Q9WTR2", "Q99LS3", "Q9WTU6", "P62046", "Q8CIV8", "Q80W93", "Q8BJ03", "Q8BMP6", "Q9D8N2", "P10853", "Q61464", "O35250", "B1AZI6", "Q9CRB2", "O35691", "Q8VD62", "Q64331", "Q6P549", "P83917", "O70133", "P62717", "O88456", "Q9JHJ3", "P17742", "P62983", "Q9D8B4", "Q922H1", "O88696", "Q78YY6", "Q9D8T0", "Q8K1L5", "P47963", "P03888", "Q9EQG9", "P40338", "Q8R001", "P62334", "Q8VEE0", "Q810V0", "Q7TND5", "P53994", "Q9QZQ8", "O35239", "P63011", "Q8BGX2", "Q9D8E6", "P58064", "P70265", "Q8R574", "Q8R5H1", "Q6ZQ58", "Q8VCH5", "P63085", "Q5FWK3", "Q9DCZ4", "Q8JZR0", "P32261", "Q9DCZ1", "Q8C167", "Q922B9", "Q3U182", "Q9EP89", "Q9CSH3", "Q8K2T1", "Q99KP6", "Q61586", "Q64674", "Q8CH25", "Q9CYP7", "Q7TMQ7", "Q01341", "Q80TV8", "Q8K4L3", "Q9QYC7", "Q8BJ34", "Q8JZQ9", "O70194", "Q8BKC5", "Q6P9R2", "Q6P9R1", "Q8BU11", "Q6PA06", "P98083", "P26443", "Q3URQ0", "Q8BV57", "Q6ZWZ2", "Q60520", "B1AVY7", "Q8BH74", "Q9D1P4", "Q8VDQ9", "Q811D0", "Q9DAM7", "Q9R060", "Q99LN9", "O89086", "Q8BRN9", "Q8R3Y8", "Q8C7E9", "Q3U4G3", "P24527", "P18052", "Q7TQH0", "Q9Z2Q6", "Q8VDC0", "Q8BGT7", "Q9CPR4", "P62827", "P47857", "Q99N85", "Q8BVI5", "Q9QYE6", "Q9QUR6", "Q921M4", "O08808", "Q3UM18", "Q9QY40", "Q6P5F9", "O55029", "Q8BIP0", "Q9WTX5", "Q61584", "P32067", "Q9Z223", "Q6P1H6", "Q8K3J1", "P70698", "Q148V7", "Q9ESP1", "P28798", "Q8BMB3", "Q8K4R9", "O35448", "Q8BFU3", "Q8K0C1", "P54923", "Q7TPV4", "B2RRE7", "Q9WUP7", "O70591", "Q9JMG1", "Q8K2C6", "Q9CQJ4", "P42337", "P70295", "Q8BM72", "Q8CGC7", "Q99NB9", "Q80XQ2", "Q8R0S2", "Q9ERB0", "Q69ZR2", "Q9DBU3", "Q3UHX0", "P28656", "Q9D7A8", "Q925T6", "Q80Y81", "Q60865", "Q9D1K7", "Q8C7R4", "Q8CG76", "Q8BI84", "Q62376", "Q9D1J3", "Q5SYD0", "Q9Z315", "Q9WV03", "Q8BJW6", "Q99LR1", "Q8K2H2", "Q2L4X1", "Q9CZ28", "Q80YR4", "Q91VK4", "Q8R010", "P12367", "P13439", "Q8BIJ7", "Q8BL66", "O35286", "Q9JK48", "Q6NZN0", "Q8K212", "P47758", "Q8BYY4", "P58742", "Q9Z2Y8", "Q9JKX6", "Q9D4V7", "Q9QXB9", "Q9JJN0", "Q91VC3", "Q60749", "Q6PHZ5", "Q8R332", "Q14AX6", "P48036", "Q9R1R2", "Q5QD15", "P34884", "Q9D0Z3", "Q689Z5", "Q9ESJ4", "O08663", "Q3UVG3", "Q6EDY6", "Q8R1B4", "A2A6T1", "Q922Q4", "Q8BFR4", "Q9WUL7", "Q5FWH2", "Q05D44", "Q9CQ48", "Q08288", "Q99LB6", "P47739", "Q5U458", "Q9CR41", "Q80UY2", "Q6NVF9", "P51163", "Q80X50", "Q8BWF0", "Q8BTI8", "Q62203", "O54794", "O35609", "Q8BNV1", "O70251", "Q921J2", "Q9R0Q3", "Q3TCN2", "Q0VGB7", "O88522", "O35226", "Q8C0E3", "Q9R0M5", "Q69ZU8", "Q61183", "Q7TN29", "Q9WVE8", "Q5SUE8", "Q9Z130", "Q8VDG3", "Q8R2Y2", "Q9R0Q7", "Q80X41", "P63101", "P09602", "Q9QZB7", "Q8BH60", "P11881", "P46096", "Q60864", "Q9Z2Q5", "Q3TTA7", "Q3V009", "P61211", "Q9CPX6", "Q80SY5", "Q8CBY8", "P59326", "O88531", "Q3UMB9", "P62911", "O35598", "Q9JJK2", "Q61166", "P97287", "O35382", "P70335", "Q9CZU3", "Q9JK23", "Q9CS42", "P20352", "Q91VT4", "P97393", "Q9Z0L8", "Q8K301", "Q9R0X4", "P12970", "P61202", "Q3UHA3", "P61965", "P70398", "P63260", "Q8K2K6", "Q6PGG6", "Q80W00", "Q8R2Y8", "Q9D0L4", "P50247", "Q61136", "Q9EPU0", "P60898", "Q9WUA2", "Q4KML4", "Q9DB25", "Q8K1R3", "Q9JL16", "Q9CR27", "Q3UMU9", "Q8BPG6", "Q9D3E6", "Q8K354", "Q58A65", "Q91ZX7", "Q6ZWX6", "Q62348", "Q60973", "P04223", "Q8VCL2", "Q9WV86", "Q3UZ39", "P49312", "Q8K2A8", "P35279", "Q8C754", "Q6ZPS6", "Q80W68", "P23591", "Q8VC48", "P54726", "Q8BX94", "Q3UM45", "Q64213", "Q6P9Z1", "P61971", "Q8BJ05", "P0C8B4", "P54728", "Q8R0K4", "Q8VDZ4", "Q791T5", "Q61074", "A2ADA5", "Q9DC71", "Q8C078", "Q9DCB8", "Q9D4V0", "Q9CQI7", "P53811", "Q4KWH5", "Q497V5", "Q61140", "Q8K120", "Q9JM58", "P52623", "Q8BU03", "O54984", "Q8K274", "Q7TNV0", "P33215", "Q9CY34", "Q8C2K1", "Q05860", "P53995", "Q8VDG5", "O35166", "Q3U0V1", "O08739", "O70439", "P70245", "Q8VH51", "Q9Z1P4", "Q922R1", "Q5SV80", "B2RXR6", "P61327", "P18608", "Q9QYF9", "Q91W67", "Q8CFI7", "P29758", "Q61990", "P26883", "Q9Z108", "Q8VHK9", "P84091", "Q9CR02", "Q6NSR8", "O08997", "Q9DBS5", "Q8C181", "Q8BFP9", "Q8VCQ3", "Q9QYB5", "P49722", "P51125", "O88569", "P11031", "P45481", "Q78HU3", "Q60668", "Q8R143", "Q8K4X7", "Q6ZQB6", "P61600", "Q9Z0G0", "Q9D8W5", "Q9DCI3", "Q9CR20", "Q8CFJ9", "Q8CI08", "Q8K4B0", "E2JF22", "Q8VEH5", "P48453", "Q9JHE3", "O08664", "P61021", "Q8CHW4", "Q9CQ10", "Q8K0V4", "P28740", "Q3B7Z2", "Q6RHR9", "O09111", "Q9D7B6", "Q9DC48", "P85094", "P62204", "Q9EQP2", "Q60598", "Q9DB90", "Q9DB73", "O35316", "Q8BKE9", "Q8CF66", "Q9Z2D1", "F8VQB6", "Q9EQH2", "Q60597", "Q8CIM3", "P26043", "Q3V4B5", "Q8CG72", "P97363", "P61963", "Q9CQX5", "Q7TNE3", "Q9JHH9", "Q8BFU2", "A2RSJ4", "Q9CY58", "P47879", "Q8VEK3", "Q91Z49", "Q8BU85", "Q9Z120", "Q99LL5", "Q9R0A0", "Q8BUB4", "Q9DCB1", "Q9CR57", "Q91VY9", "P60605", "P61922", "P22315", "Q9CRT8", "Q8VHR5", "Q99PV0", "O54833", "Q9CR61", "Q9D920", "Q9EP97", "P62830", "Q9CR51", "Q922M7", "Q8WUR0", "Q9CX13", "Q9Z0U1", "P10833", "Q9WU60", "Q61037", "Q91ZN5", "P61804", "Q6PB60", "Q6P4S8", "Q6PAV2", "Q9D7M8", "Q6ZPJ3", "Q91VE0", "P63005", "P14115", "Q3UDW8", "Q91VN6", "A2AJ15", "Q8CJ96", "P43247", "Q8BHN1", "Q9ET30", "Q9WV02", "Q91WB7", "Q8CD10", "Q9Z2D6", "P62843", "Q8CBG9", "O54774", "Q9DB05", "Q922B1", "Q8BKS9", "P62751", "Q9CZX8", "P61226", "P59048", "O88477", "O08709", "Q8K2F8", "Q91ZU6", "Q80X71", "Q01405", "Q7TPH6", "Q8CGF6", "P07607", "Q5SFM8", "Q9Z1K6", "Q9QYJ3", "Q9CQI6", "P67871", "Q9R099", "Q8BW10", "Q8K2B0", "P11859", "Q9CZ42", "P50396", "P02088", "Q91XC9", "Q8VE73", "Q99MR1", "Q9D5V5", "Q3TMH2", "Q8BLF1", "Q8BNI4", "Q9JHE7", "Q9CR96", "P62270", "P51859", "O88958", "Q99JP0", "Q60716", "P52293", "Q6NZB0", "Q6GU68", "Q3UDE2", "Q8R2U0", "Q9WUM3", "Q8C3F2", "Q9WVL6", "Q924H7", "Q9CZM2", "Q3TWF6", "Q8VD65", "Q64G17", "Q9JJL8", "Q3UHN9", "Q9CQE8", "P68254", "Q921R8", "Q9Z1K5", "Q8VE18", "P42567", "Q7TPM1", "Q64520", "P02301", "Q9EQS3", "Q99PP6", "Q9CX30", "P42208", "P59222", "P60766", "Q61655", "P80313", "Q9JLJ5", "Q8CI59", "Q64737", "P49442", "Q9CR26", "Q8BG54", "Q8N9S3", "Q9CQR6", "O08599", "Q9D0F6", "Q3TXT3", "Q9DC70", "P97370", "O88676", "Q61191", "Q99KK9", "Q9CQ01", "P18760", "Q64327", "Q61420", "Q6PIU9", "Q61147", "Q8BQU3", "Q810U5", "Q3TIV5", "Q9D7E4", "O70421", "Q811U4", "O35864", "P97760", "Q7TSG2", "Q9CW79", "Q8VCN9", "Q8BU88", "Q9JME7", "P61358", "O88839", "Q9ES00", "P62488", "Q61189", "Q9JK53", "P97480", "Q9CQB2", "P97868", "Q9D0Q7", "Q8VDD9", "Q6P5G6", "P40237", "Q64430", "Q9WVR4", "P46061", "Q9Z0J0", "Q8VDK1", "Q5DTM8", "Q9JI75", "P47753", "Q9CYC5", "Q9D4H2", "P00329", "Q9Z2D0", "Q8K1N4", "O35643", "Q9EPL8", "P36536", "Q05793", "Q8VCY6", "Q6P2K6", "Q8C1D8", "P84099", "Q8BP67", "Q9CZ04", "Q6ZQ08", "Q8BSK8", "P63154", "Q91W90", "O88685", "Q8C547", "Q6GQT1", "P97825", "Q8BX17", "Q62048", "Q9D0W5", "Q08639", "O55201", "Q61937", "Q8VDV8", "Q9EPK7", "Q80YS6", "P19324", "Q9EQN3", "Q99MR6", "P62869", "Q8C854", "P63073", "Q04750", "Q3V1V3", "Q9CZ62", "Q91WG2", "Q8R1N0", "Q8BRH0", "Q8K221", "Q60649", "A2BE28", "P97494", "Q99JY4", "Q8R5L3", "Q8CFI2", "P26369", "Q8BGU5", "Q8BGQ7", "Q3TYS2", "Q91Z96", "Q9JI78", "Q9Z0F7", "Q63932", "Q9D773", "P70671", "Q60872", "Q9JLV6", "Q5U5Q9", "Q9CQ36", "Q61024", "Q8C0L8", "Q80YR5", "Q91XE4", "Q571I9", "Q3UL36", "Q80VJ2", "A6H5Y3", "Q9D924", "P20664", "Q62059", "Q03958", "Q8CGA0", "Q8R123", "A2BDX3", "Q2NL51", "Q9DCI9", "Q9EST4", "Q9D071", "P13864", "Q99J93", "Q80UJ7", "Q9WUA3", "P59235", "Q922Q2", "P99027", "Q0P678", "Q8BRM2", "Q3U487", "Q5PRF0", "Q8BYB9", "O88587", "Q8R0F6", "Q8BFZ9", "Q7TQK1", "Q8BWY7", "Q9ERA0", "Q5D1E7", "Q8BU14", "Q61572", "Q9QUR7", "Q9QUH0", "Q99KH8", "Q61171", "Q920Q4", "P55264", "Q9DBE8", "Q8BHL8", "Q9CY57", "Q6PD26", "P23708", "Q8VE19", "P55258", "Q99J09", "Q9CQL7", "Q9ESV0", "Q9CR21", "P53996", "Q9Z321", "Q80TL7", "P08074", "Q64112", "P61222", "Q8BMK4", "Q8R015", "Q9WTX2", "Q3THK7", "Q99L27", "Q7TPS5", "Q9WTX8", "Q02789", "Q9DA77", "P48377", "Q6PAM1", "Q99ME2", "Q9JKY5", "Q9DC28", "P48759", "P70302", "Q9D8X5", "Q80U87", "P27773", "Q70FJ1", "P63242", "P56695", "Q8R1K1", "P50431", "P97822", "Q8K4Z5", "Q9CR47", "Q80U95", "Q99L02", "Q9R0Q6", "Q9Z2A9", "Q8VCM4", "Q3UPF5", "P70459", "Q9QZM0", "Q8BXQ2", "Q9DBM2", "Q9D168", "Q9CX00", "Q9CQL4", "Q9D7J9", "Q9JHI7", "Q6PDL0", "Q3TCH7", "Q03963", "Q8BZ20", "Q08509", "Q925J9", "Q99K28", "Q9CZ44", "Q91W34", "P23506", "Q9CQM5", "Q61687", "P97797", "Q80U93", "Q6PFD9", "Q91V61", "Q5SF07", "Q8C6G8", "Q64337", "Q9JL35", "Q922D8", "Q8K2V6", "P70333", "O70480", "Q8BWU5", "Q8BJY1", "Q05186", "Q5SUF2", "Q9D114", "Q61554", "Q8CGF7", "P49443", "Q61595", "Q5NCQ5", "Q6PF93", "Q8N7N5", "Q8BJM3", "Q9CR09", "P41731", "Q9D7Z3", "Q60974", "P24668", "Q8C3X8", "Q61239", "Q9CWZ7", "Q3UKC1", "Q9D1F4", "Q9CZH7", "Q9CR00", "P63089", "Q9CYX7", "Q91WJ8", "Q922Y1", "P11103", "Q8C7V8", "Q8VDQ1", "Q78PG9", "O35326", "Q62465", "P11930", "Q9ESJ0", "Q3TGF2", "Q8CH18", "Q8VE88", "Q91VR7", "Q8R088", "P01029", "Q6A009", "Q8BPA8", "Q3U2A8", "P0C872", "Q9D0R8", "O88572", "P70318", "Q9CWP6", "Q8C761", "Q9CQH7", "Q9CQQ8", "Q8BFR6", "Q7TS68", "Q9D892", "Q9CQT1", "Q91Z38", "Q8BHC7", "P12265", "Q9Z2H5", "Q9D5J6", "Q8CB77", "Q80TN7", "Q9Z0M5", "Q61771", "A8C756", "Q9DBA9", "Q9ER81", "Q8BTZ7", "Q9CYU6", "Q9WTX6", "B1AY13", "Q8VC70", "A2A6Q5", "P59708", "P70452", "Q9DB34", "Q9CPZ8", "Q8BXL7", "Q61712", "O35387", "O89051", "Q924M7", "Q922U1", "Q6NXN1", "Q9EQI8", "Q8VD04", "Q6PDH0", "Q3THS6", "P51410", "Q60676", "Q9CXZ1", "Q78ZA7", "Q99N95", "Q810J8", "P62852", "P21300", "Q03719", "Q9CPT5", "P48722", "Q9CZA6", "Q8CJF7", "O88448", "P62900", "P46737", "O89090", "O35144", "P11157", "Q8BYL4", "P97432", "Q8VD12", "Q9ES46", "Q6P2B1", "Q9QXV1", "O35218", "Q9D1H7", "Q8CCP0", "Q9Z1T1", "Q9CXK9", "Q8BWH0", "Q9D2V7", "Q61508", "Q80XK6", "Q9CZH3", "Q9ER41", "Q9R0U0", "Q5XJY5", "P62855", "Q8VDQ8", "Q8VEJ4", "Q6PNC0", "Q4VBE8", "P51807", "P83877", "Q01721", "Q99KE1", "O88665", "Q9DBR0", "P97390", "P20065", "Q9DBT5", "Q8BH43", "P14824", "Q69Z89", "Q8BQZ4", "Q64191", "Q8R146", "Q8R4R6", "Q14AI6", "Q4PZA2", "Q7M6Y3", "Q9DB43", "Q9WUD1", "O09174", "Q8R2N2", "P70297", "Q9JIW9", "Q6NWW9", "Q4FZC9", "Q61334", "Q8K1J6", "Q9CZX5", "Q8BX80", "Q6DVA0", "O09117", "Q9WU42", "Q9CQN7", "Q9DC16", "P48410", "Q9QYF1", "Q8BG81", "O89017", "P40630", "Q9ERE7", "Q8K010", "P28658", "Q3UI43", "Q8BP47", "Q9QYL7", "Q91ZW2", "P61759", "Q8BSI6", "P24788", "O09172", "Q80US4", "A2AGH6", "Q5XF89", "Q80YD1", "Q923D2", "P49586", "G5E870", "Q61216", "P32233", "P35762", "Q99PW4", "Q9JK38", "O35900", "Q8CFE4", "P43275", "Q3U1V6", "Q3U7U3", "P46978", "Q9Z204", "Q06138", "Q80UW8", "Q3TIX9", "Q9CYZ2", "P45377", "P51175", "P68181", "Q91VJ4", "Q8VDG7", "Q9CYI4", "Q9D8P4", "Q08091", "Q9JLV1", "Q99N87", "Q6P8M1", "P45591", "Q99JI6", "P11370", "Q6NSU3", "Q9WUK4", "Q91YJ3", "P48758", "Q9JIQ3", "Q9DBR3", "Q3TPX4", "Q99J36", "P97855", "Q8BLH7", "Q91X96", "P62892", "Q9CPS7", "Q8VH37", "Q8VC19", "P62264", "O70309", "P82349", "Q9D020", "P08207", "Q9QWR8", "Q5SV66", "Q8VIJ6", "Q9WTL7", "P58854", "E9PZJ8", "Q61062", "P70268", "P57776", "P01942", "Q9CR59", "Q03350", "Q9DCT8", "Q8R0W6", "Q9DB26", "Q8BPM2", "P98203", "Q9JIB4", "Q9D9K3", "O54750", "Q61542", "Q8CHG3", "P46718", "Q9D8V0", "Q9D4H8", "Q5SVR0", "Q80WG5", "Q99104", "Q3UH06", "Q9DB10", "Q5H8C4", "Q9Z1J1", "Q3V0K9", "Q9D1G1", "Q8VIM9", "Q9JLI6", "Q8K0D7", "P97434", "Q3U1J4", "Q8CHT0", "P07724", "Q8K4Z3", "Q9D883", "Q99P31", "Q61160", "Q91W59", "Q8BUE4", "Q923D5", "P97471", "Q62245", "Q6P8I4", "Q8K2V1", "Q9EP53", "Q9EPK2", "Q9QXT0", "Q6PAQ4", "Q9CPT3", "Q8VCE2", "P56542", "Q9ET26", "Q9CXI0", "Q8BTY2", "Q8VCW8", "P28867", "Q9R0K7", "Q8VBT9", "P68037", "Q5U4C3", "Q6TXD4", "Q9JHW2", "Q9QZF2", "Q8BSQ9", "Q9R0Q9", "Q3TKT4", "P50637", "Q3UHX2", "Q61749", "Q8R422", "P70188", "P46467", "Q8BJZ3", "Q5I043", "Q7TSH2", "Q60961", "P99028", "P61082", "O88983", "Q8VCB2", "P62702", "Q921F4", "Q9CYI0", "Q9Z2C5", "Q9Z266", "O35969", "Q99LC2", "Q5ND52", "P51863", "Q9D198", "Q3TQI7", "Q9CPY1", "Q3UIA2", "Q9JJC9", "P68040", "Q8K2C9", "Q9CPU0", "A3KGS3", "P17879", "O88736", "Q9DBG3", "Q9CQJ2", "Q8CFE3", "P97352", "Q9CWZ3", "Q99N57", "Q9WTI7", "Q62313", "Q9QXN3", "Q80Y55", "Q99JX7", "Q9QXA5", "Q9R1S8", "Q9JJT9", "Q9EQ80", "Q7TPE5", "Q9QXD8", "Q9D710", "P30285", "Q9R210", "Q8BIQ5", "Q62446", "P28741", "O09159", "O35177", "Q8BWY3", "B2RXC1", "Q61165", "P62918", "Q9JMB0", "P32921", "Q80VI1", "Q99LI2", "P56212", "P31938", "Q80XU3", "Q9CWU4", "Q99JY8", "P97499", "Q05512"],
    ),
    "regulated_proteins": fields.List(
        fields.String,
        description = "All regulated protein accession in experiment",
        example = ["Q9QXS1", "Q8BTM8", "Q8VDD5", "Q80X90", "P63038", "Q05920", "E9Q555", "Q69ZN7", "P41216", "Q52KR3", "Q7TPR4", "Q5SWU9", "P11087", "Q61391", "Q01149", "Q04857", "Q5SW19", "O88492", "Q8BMS1", "O35206", "P50544", "Q02788", "O35945", "Q8BFW7", "Q9CZU6", "Q921H8", "P24270", "P28271", "O88844", "P04117", "P54310", "Q8CGN5", "Q8CC88", "P14873", "P32020", "Q8BWT1", "P26231", "Q3U7R1", "P13707", "Q9D819", "Q08857", "Q71LX4", "P25085", "Q99KQ4", "Q62523", "Q9WTP6", "P35564", "Q8C129", "Q9R0E1", "P06801", "P12382", "P02469", "P05622", "Q9WVK4", "Q99JY0", "P52825", "Q62417", "P30999", "Q9ERG0", "P25206", "Q99K51", "P59017", "Q9WV92", "P37040", "P51660", "O54724", "Q8BMF4", "O08749", "P08122", "Q91YH5", "Q9ET54", "Q8BH64", "P08121", "Q60605", "Q9D0F9", "Q05816", "Q9QXS6", "Q8QZT1", "P51174", "B2RXS4", "Q9JHS4", "Q8CI51", "Q91YP2", "Q9CQ65", "P50136", "Q9D7N9", "P35486", "Q60847", "Q8R3V5", "Q99NH0", "Q9R0H0", "Q8C0N2", "Q9CQ62", "Q02248", "P97310", "Q9EQ20", "Q9D6R2", "Q9D281", "Q5SSZ5", "P97315", "P49717", "P31324", "Q9DAW9", "Q63918", "Q9DBN5", "P11152", "P49718", "Q8K0D5", "Q9EP71", "Q9DBR7", "Q9QVP9", "Q8BH97", "Q99NB9", "Q99L13", "Q62188", "Q8VDN2", "P37804", "Q9D051", "Q9JLM8", "Q99L88", "Q9QXZ0", "Q8JZK9", "Q9Z0P4", "Q9Z1E4", "Q2TPA8", "P97346", "Q64339", "Q9D9V3", "P28301", "O35855", "P53395", "P97449", "Q64429", "Q60994", "Q921G7", "P35922", "Q920E5", "P45952", "Q93092", "P15626", "F8VPU2", "Q8R180", "Q9Z2V4", "Q64449", "Q8BJ56", "Q8CC35", "Q3U962", "P18155", "Q9R1J0", "Q8VCI5", "Q6GQT9", "Q9CQ60", "Q6IRU2", "Q80UU9", "Q62351", "Q9QZA0", "P56395", "P51881", "P61022", "Q99JR1", "E9Q634", "Q07417", "P20108", "Q9JLJ2", "P47934", "Q91W92", "O08715", "O08528", "Q922J9", "Q8QZS1", "Q61881", "Q91YT0", "Q9DCW4", "P00493", "Q3UGP9", "Q91V41", "Q99MN9", "P81122", "Q9D0K1", "P28653", "Q925I1", "P58771", "Q91V01", "O88967", "Q3UVK0", "Q80TH2", "Q00915", "Q69Z38", "P42125", "Q9DBL9", "Q8K1N2", "Q8BVA4", "Q9R1Z8", "D0QMC3", "Q9WVB0", "Q8BWW4", "Q9DB20", "Q8CAY6", "Q64433", "Q8R550", "Q99KI3", "Q6ZQ73", "Q2PZL6", "Q8K2Z4", "Q91YM4", "P67778", "Q99LP6", "Q9QUJ7", "Q9WVA4", "O35129", "P70336", "P18654", "P55096", "Q61425", "Q9CR62", "Q9D6Y7", "P16332", "Q7TNP2", "Q3U0J8", "Q99M87", "Q8BIW1", "Q08093", "Q9CPQ3", "P58774", "Q3TMW1", "Q3UMC0", "P10493", "P16406", "Q8CIB5", "Q9EQU5", "Q8R2Y2", "P62843", "Q05186", "Q8K1N1", "Q9D517", "P53690", "Q5XG73", "Q8K3K7", "Q8K0C4", "P02463", "O35405", "P54116", "Q8R059", "Q99J87", "Q8BH95", "P01027", "O70423", "P97465", "Q60634", "P98078", "O70583", "P50518", "P70404", "Q8CG48", "Q6ZPJ0", "Q8K124", "Q9ER88", "P14152", "Q8CFI0", "P98192", "Q9ERT9", "P97379", "Q60823", "Q6P3A8", "O70566", "Q9D1L0", "Q6R891", "O35459", "Q9QYC0", "Q8BLN5", "Q8VEE1", "P48024", "P59325", "Q60790", "Q8BH79", "Q62000", "P58044", "Q8C5Q4", "P12787", "Q8BHD7", "Q91Z67", "Q8R1X6", "Q8R050", "Q9DB15", "Q9CXW2", "Q99J47", "P51912", "Q8CGU1", "Q8K1Z0", "Q9CXE7", "Q8BML1", "Q8BZA9", "Q64516", "O70503", "Q8BKZ9", "P62075", "P50172", "Q99J56", "Q5SNZ0", "Q9DBE0", "Q5F2E8", "Q9Z191", "Q9DBS9", "Q9QXE0", "Q9QXG4", "P09671", "Q8BVI4", "Q9CYV5", "Q8K4Q8", "Q9QWY8", "Q6ZQL4", "Q8K2C7", "Q8C156", "Q8CCF0", "Q61107", "Q8R1G6", "Q9JLC4", "Q8BPB5", "Q80WQ2", "Q61245", "Q08619;Q8CGE8", "Q9WUM4", "Q9Z2U1", "Q8VCF0", "Q920A7", "Q61490", "Q9QXL2", "Q91YR1", "Q925B0", "Q9DB27", "Q8BYK6", "P70227", "Q8CDJ8", "Q8R1U1", "Q8VE22", "Q9ER38", "Q9D7V9", "Q99JF5", "Q8CCK0", "O09061", "Q91V76", "Q9QY24", "Q9D1J1", "Q06335", "Q9QYJ0", "Q8CG47", "Q9QXY9", "Q91YW3", "Q9DCG2", "Q8C4J7", "Q8BWM0", "Q99KW3", "Q9DAI2", "Q8K0C9", "Q5RKZ7", "P49070", "Q8VEH3;Q9CQW2", "Q9ESW4", "P47713", "P29533", "Q9QZC8", "Q5SSW2", "Q80XI3", "O35382", "Q8CFI7", "Q9DB73", "Q9R0A0", "P97822", "Q62219", "Q5U4D8", "Q60953", "Q8C1E7", "Q9CQ19", "P62737;P68033", "A2AAY5", "P13595", "Q9QZD8", "Q9QYR9", "P18872", "Q9D880", "Q61009", "Q91V64", "P17809", "P82347", "P63328", "Q3UJD6", "Q71FD7", "P26645", "Q9CQC9", "Q6PIE5", "Q8R2Q4", "Q9CQX2", "Q9QZH6", "Q8BMS9", "Q04899", "Q8C3W1", "P11440", "Q01279", "Q91YX5", "A2AF47", "Q03173", "Q8R5K4", "P54227", "P20444", "Q3TJD7", "Q9CQI3", "Q61703", "Q9CZ52", "Q9DCX2", "O35099", "O35343", "Q8BWS5", "Q9R062", "Q71RI9", "Q91VA6", "Q9CZR8", "O35215", "Q8VE62", "Q8BQ30", "Q8CFH6", "Q9Z0P5", "Q9QZK2", "Q9D1K2", "Q9D7S7", "Q9WVL0", "Q9JI39", "Q8BKI2", "Q9Z0H1", "Q99LJ5", "B1AZP2", "Q9DBP5", "Q9Z0R9", "Q3TH73", "Q68FH4", "Q61398", "Q6PHU5", "P30681", "Q8R2G6", "Q99NB8", "Q64345", "P97927", "Q8BSF4", "Q91YY4", "Q5F2F2", "Q61235", "Q9CX97", "Q8BNU0", "Q8R3D1", "P97493", "Q3TBW2", "Q9D8C4", "Q9CQT5", "P11862", "P16125", "Q5FWI3", "Q9R1P3", "Q501J7", "Q6P4S6", "Q9D0I8", "Q9QYK7", "Q8R3N1", "Q9DB60", "P15092", "Q8R151", "G5E8K5", "Q8VEB4", "Q6ZQK5", "Q8R409", "Q8BJZ4", "Q9CR68", "Q9D8B3", "Q8CI61", "Q640M1", "Q99JB2", "Q8BT07", "Q8VI75", "Q6IR34", "Q91YR9", "Q8R2W9", "P28481", "P70399", "Q9D172", "Q6NV83", "Q9JKX4", "Q05CL8", "O88456", "Q8R5H1", "Q922B9", "P62827", "P12367", "Q3UMB9", "Q8VCL2", "P63005", "Q9Z2D6", "Q8VE73", "Q05793", "Q69Z89", "Q9Z2A7", "O08756", "Q9DBW0", "Q923Z0", "Q9JJE7", "O55239", "Q9DB76", "Q8K5B2", "O08573", "Q61739", "Q68FL4", "Q58NB6", "Q8BQS5", "Q9CZU4", "P97314", "Q91VS7", "Q8C838", "Q91WK5", "P41778", "Q9CR98", "Q9QUN9", "P31786", "Q9JKB3", "Q8VCN5", "Q99MQ1", "Q9CZE3", "Q8C7V3", "Q8VEE4", "Q9D365", "Q8BGI4", "Q3TBT3", "Q8VE99", "Q9R008", "Q62356", "Q3UFY8", "Q6DID7", "Q99LD8", "P27048;P63163;Q925U4", "Q99LJ0", "Q9D1N9", "Q99KR7", "A6X8Z5", "Q9WVA2", "Q8BHB4", "Q8C6B9", "Q8R0N6", "Q4LDD4", "P52432", "Q80Y17", "Q9QX11", "Q8CCI5", "Q3TN34", "Q91VU0", "Q3TDD9", "P19426", "Q8K1I7", "Q8C0L9", "Q8C6I2", "P62320", "Q9CQQ7", "Q80V53", "Q9D023", "P63087", "Q8R138", "Q80Y14", "Q60714", "Q3TB82", "Q8BJH1", "Q9D7G0", "P25799", "O88441", "Q8C033", "Q70E20", "P39053", "Q80TA1", "O54890", "Q3UBZ5", "Q9DBX2", "Q2KN98", "Q80ZJ1", "Q3V1L4", "Q8CI95", "Q9EPB5", "P22682", "Q61810", "Q60838", "P60824", "P16330", "Q8R3F5", "Q61733", "O88286", "Q9Z0S1", "P70257", "Q9R1Q9", "Q9JI13", "Q8BPU7", "Q9JKV1", "Q64261", "Q9Z0F8", "P27641", "A6H5X4", "Q9WUZ9", "Q8K157", "Q8BTY1", "Q9JIM1", "P62715;P63330", "Q8C166", "P10518", "Q8BK03", "Q9D1G2", "Q69ZW3", "Q3TMX7", "Q8K298", "Q9CZL2", "Q3UTJ2", "Q9EPK6", "Q9WU56", "O70433", "Q921G6", "Q3UH60", "Q91WS0", "Q99LH1", "Q9CX60", "A2ASS6", "P10923", "Q9CYA0", "O88322", "Q9QZH3", "P97819", "P35951", "Q9D2C7", "Q68FF6", "P14142", "Q9WV84", "Q3UH93", "Q8BV66", "Q3TPE9", "Q499X9", "Q9JIK9", "Q6P6L0", "Q924L1", "Q812A2", "Q921U8", "P97376", "O35710", "P50285", "Q9D6M3", "P41233", "Q6A058", "Q80UK0", "Q99PP9", "Q9D0T1", "P20357", "Q9R0Z9", "P15105", "Q99N93", "Q61824", "Q05BC3", "Q8QZR5", "Q9JJG9", "P57746", "P61025", "Q9JIZ9", "Q7TN98", "P62082", "Q8BX90", "Q8CIL4", "Q9WUQ2", "P70414", "P57716", "Q9DBA6", "Q9EPL9", "O54916", "Q6PFQ7", "Q8BZH4", "Q9CQZ5", "Q9CWU9", "Q6DIC0", "Q148V8", "Q9JI11", "O08912", "Q921K8", "O88413", "P36552", "Q62086", "Q61102", "O55106", "Q9Z2R6", "Q9EPU4", "Q3TFQ1", "Q8R1J9", "Q9ESU6", "Q9CR89", "Q8BP27", "P16675", "Q791V5", "Q8BKF1", "Q8BVU0", "Q8R0H9", "Q99LS3", "Q9D8N2", "P63011", "Q01341", "Q9QYC7", "Q3URQ0", "Q8R3Y8", "Q3U4G3", "Q7TQH0", "Q8VDC0;Q9D9J8", "Q6P1H6", "P70295", "P28656", "Q8C7R4", "Q9Z315", "Q8BIJ7", "Q3TCN2", "P61211", "O88531", "Q9R0X4", "Q8R2Y8", "P50247", "P61971", "Q8VDG5", "P84091", "Q9D7B6", "P85094", "P61922", "Q9WU60", "P61804", "Q8CD10", "Q9CQI6", "Q6GU68", "Q8VE18", "P42567", "Q8BG54", "Q61147", "Q9D071", "P55264", "Q8CGF7", "P20065", "Q9QWR8", "Q8K4Z3", "Q9D883"],
    ),
    "background_metabolites": fields.List(
        fields.String,
        description = "All metabolite ChEBI Ids in experiment (background)",
        example = ["1,3BPG", "3-Phosphoglyceric acid", "Acetylcarnitine DL", "Acetylcholine", "Adenine", "Adenosine", "Adenosine  phosphosulfate", "ADP-Glucose", "ADP-Rib", "Agmatine", "Alanine", "AMP", "Arginine", "Asparagine", "Aspartic acid", "Beta-alanine", "Carbamoylphosphate", "Carnitine", "Choline", "CMP", "Creatine", "cyclic AMP", "Cytidine", "Cytosine", "dCMP", "Deoxyguanosine", "Dephospho-CoA", "dTDP", "dTTP", "Ethanolamine", "FAD", "Galactose", "Gamma-Aminobutyric acid", "Glucosamine 6-phosphate", "Glutamic acid", "Glutamine", "Glyceraldehyde 3-phosphate", "Glycerol", "Glycerol 3-phosphate", "Glycerophosphocholine", "Glycine", "GMP", "GSH", "GSSG", "Guanine", "Guanosine", "Hexose monophosphate", "Hexose pool", "Histidine", "Homocysteine", "Hypoxanthine", "IMP", "Indole", "Inosine", "Leucine", "Isoleucine", "Lactic acid", "Lysine", "Methionine", "methylnicotinamide", "Myoinositol", "N-Acetyl-glucosamine-1-phosphate", "N-Acetylornithine", "N-acetylserine", "NAD+", "NADP+", "Nicotinic Acid", "Ornithine", "Orotidylic Acid", "Pantothenic acid", "Phenylalanine", "Phosphorylcholine", "Phosphoserine", "Pipecolic acid (DL)", "Prephenate", "Proline", "Pyroglutamic acid", "Pyruvic acid", "Riboflavin", "Ribose-5-phosphate", "S-Methyl-5-thioadenosine", "S-Adenosylhomocysteine", "S-Adenosylmethionine", "Sarcosine", "Sedoheptulose monophosphate", "Serine", "Shikimic acid", "Succinic acid", "Succinyl-CoA", "Taurine", "Thiamine", "Threonine", "Thymine", "Tryptophan", "Tyrosine", "UDP", "UDP-D-glucose", "UDP-D-glucuronate", "UDP-N-acetyl-D-glucosamine", "UMP", "Uracil", "Uridine", "Uridine diphosphate glucuronic", "Valine", "Xanthine", "Xanthosine", "Homoserine", "Erythrose 4-phosphate", "Ribulose-1,5-bisphosphate", "Aminoimidazole carboxamide ribonucleotide"],
    ),
    "regulated_metabolites": fields.List(
        fields.String,
        description = "All regulated metabolite ChEBI Ids in experiment",
        example = ["Erythrose 4-phosphate", "Acetylcholine", "Adenosine", "Riboflavin", "Guanosine", "dTDP", "Thiamine", "FAD", "Xanthosine", "UMP", "Carnitine", "Inosine", "Uracil", "Uridine", "Adenine", "ADP-Glucose", "IMP", "Cytidine", "Cytosine", "Xanthine", "Hexose monophosphate", "GMP", "AMP", "Hypoxanthine", "Pantothenic acid", "dTTP", "Ethanolamine", "Succinyl-CoA", "Dephospho-CoA"],
    ),
    "background_transcripts": fields.List(
        fields.String,
        description = "All ensembl Ids in experiment (background)",
        example = [],
    ),
    "regulated_transcripts": fields.List(
        fields.String,
        description = "All regulated ensembl Ids in experiment",
        example = [],
    ),
    "organism_taxonomy": fields.String(
        description = "Select organism (Taxonomic number), default: 'NCBITaxon:9606' (Homo sapiens)",
        enum = [v for k, v in organisms.items()],
        example = "NCBITaxon:10090",
    ),
    "domains": fields.List(
        fields.String,
        description = "Select domain(s), default: 'biological_process'",
        enum = sorted(list(d.lower().replace(" ", "_") for d in enrichment_ontologies[INIT_ORGANISM].domains)),
        example = ["biological_process"],
    ),
    "pvalue_correction": fields.String(
        description = "Select method for p-value correction, default: 'fdr_bh' (Benjamini-Hochberg)",
        enum = [c["value"] for c in correction_list],
        example = "fdr_bh",
    ),
    "term_representation": fields.String(
        description = "Term representation, default: greater",
        enum = [t["value"] for t in term_representation],
        example = "greater",
    ),
    "unrecognizable_molecules": fields.String(
        description = f"Handling of unrecognizable molecules, default: {MOLECULE_HANDLING_REMOVE}",
        enum = [m["value"] for m in molecule_handling],
        example = MOLECULE_HANDLING_REMOVE,
    ),
    "non_background_molecules": fields.String(
        description = f"Handling of non-background regulated molecules, default: {MOLECULE_HANDLING_REMOVE}",
        enum = [r["value"] for r in regulated_molecule_handling],
        example = MOLECULE_HANDLING_REMOVE,
    ),
    "bounded_fatty_acyls": fields.Boolean(
        description = f"Use bounded fatty acyls for analysis too, default: false",
        example = False,
    ),
})


@ns.route("/enrichment")
class EnrichmentResource(Resource):
    @api.expect(enrichment_model)
    def post(self):

        logger.info(f"New API access: enrichment")

        try:
            data = api.payload  # JSON body
            background_lipids = data.get("background_lipids", [])
            regulated_lipids = data.get("regulated_lipids", [])
            background_proteins = data.get("background_proteins", [])
            regulated_proteins = data.get("regulated_proteins", [])
            background_metabolites = data.get("background_metabolites", [])
            regulated_metabolites = data.get("regulated_metabolites", [])
            background_transcripts = data.get("background_transcripts", [])
            regulated_transcripts = data.get("regulated_transcripts", [])

            if type(background_lipids) != list:
                return {"error_message": "'background_lipids' needs to be a list", "result": []}, 422
            if type(regulated_lipids) != list:
                return {"error_message": "'regulated_lipids' needs to be a list", "result": []}, 422
            if type(background_proteins) != list:
                return {"error_message": "'background_proteins' needs to be a list", "result": []}, 422
            if type(regulated_proteins) != list:
                return {"error_message": "'regulated_proteins' needs to be a list", "result": []}, 422
            if type(background_metabolites) != list:
                return {"error_message": "'background_metabolites' needs to be a list", "result": []}, 422
            if type(regulated_metabolites) != list:
                return {"error_message": "'regulated_metabolites' needs to be a list", "result": []}, 422
            if type(background_transcripts) != list:
                return {"error_message": "'background_transcripts' needs to be a list", "result": []}, 422
            if type(regulated_transcripts) != list:
                return {"error_message": "'regulated_transcripts' needs to be a list", "result": []}, 422

            organism_api = data.get("organism_taxonomy", "NCBITaxon:9606")
            domains_api = data.get("domains", ["biological_process"])
            accepted_domains = {d.lower().replace(" ", "_"): d for d in enrichment_ontologies[INIT_ORGANISM].domains}
            pvalue_correction_api = data.get("pvalue_correction", "fdr_bh")
            term_representation_api = data.get("term_representation", "greater")
            unrecognizable_molecules_api = data.get("unrecognizable_molecules", MOLECULE_HANDLING_IGNORE)
            non_background_molecules_api = data.get("non_background_molecules", MOLECULE_HANDLING_REMOVE)
            bounded_fatty_acyls_api = data.get("bounded_fatty_acyls", False)

            if str(organism_api) not in organisms.values():
                return {"error_message": f"'organism_taxonomy' must be one of {set(organisms.values())}", "result": []}, 422
            organism_api = str(organism_api)

            if type(domains_api) != list:
                return {"error_message": "'domains' needs to be a list", "result": []}, 422

            if len(set(domains_api) - accepted_domains.keys()) > 0:
                return {"error_message": f"accepted 'domains' entries are {accepted_domains.keys()}", "result": []}, 422

            if len(domains_api) == 0:
                return {"error_message": "No domain(s) selected.", "result": []}, 422
            domains_api = [accepted_domains[d] for d in domains_api]

            if type(pvalue_correction_api) != str:
                return {"error_message": "'pvalue_correction' must be a string", "result": []}, 422

            if pvalue_correction_api not in set(c["value"] for c in correction_list):
                return {"error_message": f"accepted 'pvalue_correction' entries are {set(c["value"] for c in correction_list)}", "result": []}, 422

            if type(term_representation_api) != str:
                return {"error_message": "'term_representation' must be a string", "result": []}, 422

            if term_representation_api not in set(t["value"] for t in term_representation):
                return {"error_message": f"accepted 'term_representation' entries are {set(t["value"] for t in term_representation)}", "result": []}, 422

            if type(unrecognizable_molecules_api) != str:
                return {"error_message": "'unrecognizable_molecules' must be a string", "result": []}, 422

            if unrecognizable_molecules_api not in set(m["value"] for m in molecule_handling):
                return {"error_message": f"accepted 'unrecognizable_molecules' entries are {set(m["value"] for m in molecule_handling)}", "result": []}, 422

            if type(non_background_molecules_api) != str:
                return {"error_message": "'non_background_molecules' must be a string", "result": []}, 422

            if non_background_molecules_api not in set(r["value"] for r in regulated_molecule_handling):
                return {"error_message": f"accepted 'non_background_molecules' entries are {set(r["value"] for r in regulated_molecule_handling)}", "result": []}, 422

            if type(bounded_fatty_acyls_api) not in {int, bool}:
                return {"error_message": "'bounded_fatty_acyls' must be a boolean (0 | 1)", "result": []}, 422

            if type(bounded_fatty_acyls_api) == int:
                bounded_fatty_acyls_api = bounded_fatty_acyls_api > 0

            with_lipids = len(background_lipids) > 0 or len(regulated_lipids)
            with_proteins = len(background_proteins) > 0 or len(regulated_proteins)
            with_metabolites = len(background_metabolites) > 0 or len(regulated_metabolites)
            with_transcripts = len(background_transcripts) > 0 or len(regulated_transcripts)

            omics_included = [with_lipids, with_proteins, with_metabolites, with_transcripts]

            omics_lists = [
                background_lipids,
                regulated_lipids,
                background_proteins,
                regulated_proteins,
                background_metabolites,
                regulated_metabolites,
                background_transcripts,
                regulated_transcripts,
            ]

            ontology = enrichment_ontologies[organism_api]

            try:
                (
                    target_set,
                    lipidome,
                    regulated_lipids,
                    proteome,
                    regulated_proteins,
                    metabolome,
                    regulated_metabolites,
                    transcriptome,
                    regulated_transcripts,
                    background_list,
                ) = check_user_input(
                    omics_included,
                    omics_lists,
                    ontology,
                    unrecognizable_molecules_api,
                    non_background_molecules_api,
                )
            except Exception as error_message:
                return {"error_message": str(error_message), "result": []}, 422

            session = SessionEntry()
            session.time = time.time()
            session.use_bounded_fatty_acyls = bounded_fatty_acyls_api
            (
                session.search_terms,
                session.all_parent_nodes,
            ) = ontology.set_background(
                lipid_dict = lipidome,
                protein_set = proteome,
                metabolite_set = metabolome,
                transcript_set = transcriptome,
                use_bounded_fatty_acyls = session.use_bounded_fatty_acyls,
            )
            session.num_background = len(lipidome) + len(proteome) + len(metabolome) + len(transcriptome)
            session.ontology = ontology
            session.data_loaded = True
            analytics("api_enrichment_analysis")
            session.background_lipids = lipidome if with_lipids else None
            session.regulated_lipids = regulated_lipids if with_lipids else None
            session.background_proteins = proteome if with_proteins else None
            session.regulated_proteins = regulated_proteins if with_proteins else None
            session.background_metabolites = metabolome if with_metabolites else None
            session.regulated_metabolites = regulated_metabolites if with_metabolites else None
            session.background_transcripts = transcriptome if with_transcripts else None
            session.regulated_transcripts = regulated_transcripts if with_transcripts else None
            background_list.sort(key = lambda row: row["value"])
            session.background_list = background_list

            session.domains = set(domains_api)
            results = ontology.enrichment_analysis(session.search_terms, session.num_background, target_set, domains_api, term_representation_api, pvalue_correction_api)
            session.result = results

            result_data = []
            for result in results:
                expected = round((result.fisher_data[0] + result.fisher_data[1]) * (result.fisher_data[0] + result.fisher_data[2]) / sum(result.fisher_data))
                row = {
                    "domain": " | ".join(result.term.domain),
                    "term": result.term.name,
                    "termid": result.term.term_id_str,
                    "count": f"{result.fisher_data[0]} ({expected}) / {len(result.source_terms)}",
                    "pvalue": result.pvalue_corrected,
                    "log_odds_ratio": result.lor,
                }
                result_data.append(row)


            return {"error_message": "", "result": result_data}, 200

        except Exception as e:
            logger.error("".join(traceback.format_tb(e.__traceback__)))
            return {"error_message": f"{e}", "result": []}, 500



@callback(
    Output("field_regulated_lipids", "style", allow_duplicate = True),
    Output("field_up_down_regulated_lipids", "style", allow_duplicate = True),
    Output("field_regulated_proteins", "style", allow_duplicate = True),
    Output("field_up_down_regulated_proteins", "style", allow_duplicate = True),
    Output("field_regulated_metabolites", "style", allow_duplicate = True),
    Output("field_up_down_regulated_metabolites", "style", allow_duplicate = True),
    Output("field_regulated_transcripts", "style", allow_duplicate = True),
    Output("field_up_down_regulated_transcripts", "style", allow_duplicate = True),
    Output("sankey_upregulated", "style", allow_duplicate = True),
    Output("sankey_downregulated", "style", allow_duplicate = True),
    Input("separate_updown_switch", "checked"),
    State("field_regulated_lipids", "style"),
    State("field_up_down_regulated_lipids", "style"),
    State("field_regulated_proteins", "style"),
    State("field_up_down_regulated_proteins", "style"),
    State("field_regulated_metabolites", "style"),
    State("field_up_down_regulated_metabolites", "style"),
    State("field_regulated_transcripts", "style"),
    State("field_up_down_regulated_transcripts", "style"),
    State("sankey_upregulated", "style"),
    State("sankey_downregulated", "style"),
    prevent_initial_call = True,
)
def separate_updown_switch_changed(
    separate_updown_switch,
    field_regulated_lipids_style,
    field_up_down_regulated_lipids_style,
    field_regulated_proteins_style,
    field_up_down_regulated_proteins_style,
    field_regulated_metabolites_style,
    field_up_down_regulated_metabolites_style,
    field_regulated_transcripts_style,
    field_up_down_regulated_transcripts_style,
    sankey_upregulated_style,
    sankey_downregulated_style,
):
    field_regulated_lipids_style["display"] = "none" if separate_updown_switch else "flex"
    field_up_down_regulated_lipids_style["display"] = "flex" if separate_updown_switch else "none"
    field_regulated_proteins_style["display"] = "none" if separate_updown_switch else "flex"
    field_up_down_regulated_proteins_style["display"] = "flex" if separate_updown_switch else "none"
    field_regulated_metabolites_style["display"] = "none" if separate_updown_switch else "flex"
    field_up_down_regulated_metabolites_style["display"] = "flex" if separate_updown_switch else "none"
    field_regulated_transcripts_style["display"] = "none" if separate_updown_switch else "flex"
    field_up_down_regulated_transcripts_style["display"] = "flex" if separate_updown_switch else "none"
    sankey_upregulated_style["display"] = "flex" if separate_updown_switch else "none"
    sankey_downregulated_style["display"] = "flex" if separate_updown_switch else "none"

    return (
        field_regulated_lipids_style,
        field_up_down_regulated_lipids_style,
        field_regulated_proteins_style,
        field_up_down_regulated_proteins_style,
        field_regulated_metabolites_style,
        field_up_down_regulated_metabolites_style,
        field_regulated_transcripts_style,
        field_up_down_regulated_transcripts_style,
        sankey_upregulated_style,
        sankey_downregulated_style,
    )



