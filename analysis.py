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
from flask import Flask, request, jsonify, redirect, session as uuid_session
import dash_mantine_components as dmc
import plotly.graph_objs as go
import dash_ag_grid as dag
from dash_iconify import DashIconify
import json
import numpy as np
import pandas as pd
import io, os
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
import numpy as np
from ApiRequest import *
import base64
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.backends import default_backend
from cryptography.fernet import Fernet


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
        'Homo sapiens': 'NCBITaxon:9606',
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
INNER_CIRCLE = 350
CIRCLE_WIDTH = 350
BAR_SORTING_PVALUE = "pvalue"
BAR_SORTING_SIMILARITY = "sim"

LIPIDS_COLOR = "#A3DC9A"
PROTEINS_COLOR = "#DEE791"
METABOLITES_COLOR = "#FFF9BD"
TRANSCRIPTS_COLOR = "#FFD6BA"
NoneType = type(None)

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

server = Flask(__name__)

@server.route("/submit", methods=["POST"], strict_slashes=False)
def submit():
    data = request.get_json()
    if not data:
        return "No JSON received", 400

    if "uid" not in uuid_session:
        return jsonify({"status": "ok"})
    session_id = uuid_session.get("uid")
    if not session_id:
        return jsonify({"status": "error", "reason": "no uid"}), 403

    if session_id not in sessions:
        sessions[session_id] = SessionEntry()
        logger.info(f"New session via submit: {session_id}")

    session = sessions[session_id]
    ui = session.ui
    organism = None
    domain_values = []

    def is_true(value):
        return (type(value) == bool and param_value) or (type(value) == str and value.lower() in {"true", "1", "yes"})

    try:
        for param_key, param_value in data.items():
            match param_key:
                case "separate_updown_switch":
                    ui["separate_updown_switch"] = is_true(param_value)
                case "all_lipids":
                    ui["all_lipids_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "regulated_lipids":
                    ui["regulated_lipids_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "upregulated_lipids":
                    ui["upregulated_lipids_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "downregulated_lipids":
                    ui["downregulated_lipids_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "all_proteins":
                    ui["all_proteins_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "regulated_proteins":
                    ui["regulated_proteins_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "upregulated_proteins":
                    ui["upregulated_proteins_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "downregulated_proteins":
                    ui["downregulated_proteins_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "all_metabolites":
                    ui["all_metabolites_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "regulated_metabolites":
                    ui["regulated_metabolites_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "upregulated_metabolites":
                    ui["upregulated_metabolites_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "downregulated_metabolites":
                    ui["downregulated_metabolites_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "all_transcripts":
                    ui["all_transcripts_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "regulated_transcripts":
                    ui["regulated_transcripts_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "upregulated_transcripts":
                    ui["upregulated_transcripts_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "downregulated_transcripts":
                    ui["downregulated_transcripts_list"] = "\n".join(param_value) if type(param_value) == list else param_value
                case "use_bounded_fatty_acyls":
                    session.use_bounded_fatty_acyls = is_true(param_value)
                case "use_lipids":
                    ui["with_lipids"] = is_true(param_value)
                case "use_proteins":
                    ui["with_proteins"] = is_true(param_value)
                case "use_metabolites":
                    ui["with_metabolites"] = is_true(param_value)
                case "use_transcripts":
                    ui["with_transcripts"] = is_true(param_value)
                case "organism":
                    if param_value in enrichment_ontologies:
                        ui["select_organism"] = param_value
                        organism = param_value
                case "molecule_handling":
                    if param_value in {MOLECULE_HANDLING_ERROR, MOLECULE_HANDLING_REMOVE, MOLECULE_HANDLING_IGNORE}:
                        ui["select_molecule_handling"] = param_value
                case "regulated_molecule_handling":
                    if param_value in {MOLECULE_HANDLING_ERROR, MOLECULE_HANDLING_REMOVE}:
                        ui["select_regulated_molecule_handling"] = param_value
                case "term_representation":
                    if param_value in {t["value"] for t in term_representation}:
                        ui["select_term_representation"] = param_value
                case "test_method":
                    if param_value in {c["value"] for c in correction_list}:
                        ui["select_test_method"] = param_value
                case "domains":
                    domain_values = param_value if type(param_value) == list else [param_value]
                case "use_upregulated_lipids":
                    ui["check_textarea_upregulated_lipids"] = is_true(param_value)
                case "use_downregulated_lipids":
                    ui["check_textarea_downregulated_lipids"] = is_true(param_value)
                case "use_upregulated_proteins":
                    ui["check_textarea_upregulated_proteins"] = is_true(param_value)
                case "use_downregulated_proteins":
                    ui["check_textarea_downregulated_proteins"] = is_true(param_value)
                case "use_upregulated_metabolites":
                     ui["check_textarea_upregulated_metabolites"] = is_true(param_value)
                case "use_downregulated_metabolites":
                    ui["check_textarea_downregulated_metabolites"] = is_true(param_value)
                case "use_upregulated_transcripts":
                    ui["check_textarea_upregulated_transcripts"] = is_true(param_value)
                case "use_downregulated_transcripts":
                    ui["check_textarea_downregulated_transcripts"] = is_true(param_value)

    except:
        pass

    if organism and domain_values:
        ontology = enrichment_ontologies[organism]
        ontology_domains = sorted(ontology.domains)
        ui["select_domains"] = [v for v in domain_values if v in ontology_domains]

    return jsonify({"status": "ok", "uid": session_id})



# Create the Dash app
app = Dash("app", update_title = None, server = server)
app.secret_key = "".join(chr(c) for c in np.random.randint(33, 127, 100))
server.secret_key = app.secret_key
app.title = APPLICATION_SHORT_TITLE


# --- Key derivation ---
def derive_key(password: str, salt: bytes) -> bytes:
    kdf = PBKDF2HMAC(
        algorithm = hashes.SHA256(),
        length = 32,
        salt = salt,
        iterations = 100_000,
        backend = default_backend()
    )
    return base64.urlsafe_b64encode(kdf.derive(password.encode()))

# --- Encrypt ---
def encrypt_uuid(uuid_str: str, password: str) -> str:
    salt = os.urandom(16)  # random salt
    key = derive_key(password, salt)
    f = Fernet(key)
    token = f.encrypt(uuid_str.encode())
    return base64.urlsafe_b64encode(salt + token).decode()


# --- Decrypt ---
def decrypt_uuid(encrypted_data: str, password: str) -> str:
    data = base64.urlsafe_b64decode(encrypted_data.encode())
    salt = data[:16]
    token = data[16:]
    key = derive_key(password, salt)
    f = Fernet(key)
    return f.decrypt(token).decode()





@server.before_request
def set_udata():
    uid = request.args.get("uid")
    if uid:
        uuid_session["uid"] = uid
    else:
        uid = request.cookies.get("uid")
        if uid:
            try:
                uuid_session["uid"] = decrypt_uuid(uid, app.secret_key)
            except Exception as e:
                uuid_session["uid"] = str(uuid.uuid4())  # brand new user
        elif "uid" not in uuid_session:
            uuid_session["uid"] = str(uuid.uuid4())  # brand new user


@server.after_request
def save_uid_cookie(response):
    if "uid" in uuid_session:
        response.set_cookie(
            "uid",
            encrypt_uuid(uuid_session["uid"], app.secret_key),
            max_age=60*60*24*365, # 1 year
            httponly=True,
            samesite="Lax"
        )
    return response


api = Api(
    server,
    version = get_latest_tag(),
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
                html.Div([
                    html.Div(
                        dcc.Graph(
                            id = "barplot_terms",
                            config = {**plotly_config, **{"responsive": True}},
                            style = {
                                "height": "100%",
                                "width": "100%",
                            },
                        ),
                        style = {
                            "resize": "both",
                            "overflow": "hidden",
                            "width": "100%",
                            "height": f"{CIRCLE_WIDTH * 2}px",
                        }
                    ),
                    html.Div(
                        [
                            dmc.NumberInput(
                                id = "barplot_numberinput_font_size",
                                label = "Font size:",
                                value = 14,
                                min = 1,
                                max = 20,
                            ),
                            dmc.NumberInput(
                                id = "barplot_numberinput_connect_ths",
                                label = "Connectivity Threshold:",
                                value = 80,
                                min = 0,
                                max = 100,
                                style = {"marginTop": "5px"},
                            ),
                            dmc.Select(
                                id = "barplot_select_name",
                                label = "Select bar label:",
                                data = [{"value": "id", "label": "Term ID"}, {"value": "name", "label": "Term name"}],
                                value = "id",
                                style = {"marginTop": "5px"},
                            ),
                            dmc.Select(
                                id = "barplot_select_sorting",
                                label = "Select sorting:",
                                data = [{"value": BAR_SORTING_PVALUE, "label": "p-value"}, {"value": BAR_SORTING_SIMILARITY, "label": "Molecule similarity in term"}],
                                value = BAR_SORTING_PVALUE,
                                style = {"marginTop": "5px"},
                            ),
                            dmc.Text(
                                "Figure legend:",
                                fw = 500,
                                style = {"width": "350px", "marginTop": "25px", "fontSize": 14},
                            ),
                            dmc.Image(
                                src = dash.get_app().get_asset_url(f"barplot-legend.png"),
                                style = {"width": "350px", "marginTop": "10px"},
                            ),
                        ],
                        id = "barplot_controls",
                        style = {"marginTop": "10px"},
                    ),
                    ], style={"display": "flex"}, id = "barplot_terms_wrapper",
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
                                "flex": "5",
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
                    "Providing a web service for multiomics enrichment analyses. This website uses cookies that are essential for its operation. These cookies enable core functionalities for, i.e., session storing. By using this website, you consent to the use of these essential cookies. No personal information are stored in these cookies.",
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
                                        html.Div(
                                            [
                                                dmc.Title(
                                                    "All up-regulated lipid names in experiment",
                                                    order = 5,
                                                ),
                                                dmc.Checkbox(
                                                    id = "check_textarea_upregulated_lipids",
                                                    checked = True,
                                                ),
                                            ],
                                            style = {
                                                "width": "100%",
                                                "marginTop": "10px",
                                                "display": "flex",
                                                "alignItems": "center",
                                                "gap": "10px",
                                            },
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

                                        html.Div(
                                            [
                                                dmc.Title(
                                                    "All down-regulated lipid names in experiment",
                                                    order = 5,
                                                ),
                                                dmc.Checkbox(
                                                    id = "check_textarea_downregulated_lipids",
                                                    checked = True,
                                                ),
                                            ],
                                            style = {
                                                "width": "100%",
                                                "marginTop": "10px",
                                                "display": "flex",
                                                "alignItems": "center",
                                                "gap": "10px",
                                            },
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
                                            "All regulated protein accessions in experiment",
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
                                        html.Div(
                                            [
                                                dmc.Title(
                                                    "All up-regulated protein accessions in experiment",
                                                    order = 5,
                                                ),
                                                dmc.Checkbox(
                                                    id = "check_textarea_upregulated_proteins",
                                                    checked = True,
                                                ),
                                            ],
                                            style = {
                                                "width": "100%",
                                                "marginTop": "10px",
                                                "display": "flex",
                                                "alignItems": "center",
                                                "gap": "10px",
                                            },
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

                                        html.Div(
                                            [
                                                dmc.Title(
                                                    "All down-regulated protein accession in experiment",
                                                    order = 5,
                                                ),
                                                dmc.Checkbox(
                                                    id = "check_textarea_downregulated_proteins",
                                                    checked = True,
                                                ),
                                            ],
                                            style = {
                                                "width": "100%",
                                                "marginTop": "10px",
                                                "display": "flex",
                                                "alignItems": "center",
                                                "gap": "10px",
                                            },
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
                                        html.Div(
                                            [
                                                dmc.Title(
                                                    "All up-regulated metabolite ChEBI Ids in experiment",
                                                    order = 5,
                                                ),
                                                dmc.Checkbox(
                                                    id = "check_textarea_upregulated_metabolites",
                                                    checked = True,
                                                ),
                                            ],
                                            style = {
                                                "width": "100%",
                                                "marginTop": "10px",
                                                "display": "flex",
                                                "alignItems": "center",
                                                "gap": "10px",
                                            },
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

                                        html.Div(
                                            [
                                                dmc.Title(
                                                    "All down-regulated metabolite ChEBI Ids in experiment",
                                                    order = 5,
                                                ),
                                                dmc.Checkbox(
                                                    id = "check_textarea_downregulated_metabolites",
                                                    checked = True,
                                                ),
                                            ],
                                            style = {
                                                "width": "100%",
                                                "marginTop": "10px",
                                                "display": "flex",
                                                "alignItems": "center",
                                                "gap": "10px",
                                            },
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
                                        html.Div(
                                            [
                                                dmc.Title(
                                                    "All up-regulated ensembl Ids in experiment",
                                                    order = 5,
                                                ),
                                                dmc.Checkbox(
                                                    id = "check_textarea_upregulated_transcripts",
                                                    checked = True,
                                                ),
                                            ],
                                            style = {
                                                "width": "100%",
                                                "marginTop": "10px",
                                                "display": "flex",
                                                "alignItems": "center",
                                                "gap": "10px",
                                            },
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

                                        html.Div(
                                            [
                                                dmc.Title(
                                                    "All down-regulated ensembl Ids in experiment",
                                                    order = 5,
                                                ),
                                                dmc.Checkbox(
                                                    id = "check_textarea_downregulated_transcripts",
                                                    checked = True,
                                                ),
                                            ],
                                            style = {
                                                "width": "100%",
                                                "marginTop": "10px",
                                                "display": "flex",
                                                "alignItems": "center",
                                                "gap": "10px",
                                            },
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
                                    data = sorted(enrichment_ontologies[INIT_ORGANISM].domains),
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
                            html.Div([
                                dmc.Button(
                                    "Reset",
                                    id = "button_reset",
                                    variant = "outline",
                                    style = {"marginRight": "10px"},
                                ),
                                dmc.Button(
                                    "Run analysis",
                                    id = "button_run_enrichment",

                                ),
                            ], style = {"marginLeft": "auto"}),
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
    data = sorted(ontology.domains)
    values = [v for v in domain_values if v in data]

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
    Output("check_textarea_upregulated_lipids", "checked", allow_duplicate = True),
    Output("check_textarea_downregulated_lipids", "checked", allow_duplicate = True),
    Output("check_textarea_upregulated_proteins", "checked", allow_duplicate = True),
    Output("check_textarea_downregulated_proteins", "checked", allow_duplicate = True),
    Output("check_textarea_upregulated_metabolites", "checked", allow_duplicate = True),
    Output("check_textarea_downregulated_metabolites", "checked", allow_duplicate = True),
    Output("check_textarea_upregulated_transcripts", "checked", allow_duplicate = True),
    Output("check_textarea_downregulated_transcripts", "checked", allow_duplicate = True),
    Input("url", "pathname"),
    prevent_initial_call = True,
)
def load_uid(_):
    session_id = uuid_session.get("uid")
    if session_id not in sessions:
        sessions[session_id] = SessionEntry()
        logger.info(f"New session: {session_id}")
        return (session_id, *([no_update] * 55))

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

    check_textarea_upregulated_lipids = ui["check_textarea_upregulated_lipids"] if "check_textarea_upregulated_lipids" in ui else True
    check_textarea_downregulated_lipids = ui["check_textarea_downregulated_lipids"] if "check_textarea_downregulated_lipids" in ui else True
    check_textarea_upregulated_proteins = ui["check_textarea_upregulated_proteins"] if "check_textarea_upregulated_proteins" in ui else True
    check_textarea_downregulated_proteins = ui["check_textarea_downregulated_proteins"] if "check_textarea_downregulated_proteins" in ui else True
    check_textarea_upregulated_metabolites = ui["check_textarea_upregulated_metabolites"] if "check_textarea_upregulated_metabolites" in ui else True
    check_textarea_downregulated_metabolites = ui["check_textarea_downregulated_metabolites"] if "check_textarea_downregulated_metabolites" in ui else True
    check_textarea_upregulated_transcripts = ui["check_textarea_upregulated_transcripts"] if "check_textarea_upregulated_transcripts" in ui else True
    check_textarea_downregulated_transcripts = ui["check_textarea_downregulated_transcripts"] if "check_textarea_downregulated_transcripts" in ui else True

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
        check_textarea_upregulated_lipids,
        check_textarea_downregulated_lipids,
        check_textarea_upregulated_proteins,
        check_textarea_downregulated_proteins,
        check_textarea_upregulated_metabolites,
        check_textarea_downregulated_metabolites,
        check_textarea_upregulated_transcripts,
        check_textarea_downregulated_transcripts,
    )


@callback(
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("separate_updown_switch", "checked", allow_duplicate = True),
    Output("num_enrichment_terms", "children", allow_duplicate = True),
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
    Output("multiselect_filter_molecules", "data", allow_duplicate = True),
    Output("select_term_representation", "value", allow_duplicate = True),
    Output("select_test_method", "value", allow_duplicate = True),
    Output("select_domains", "value", allow_duplicate = True),
    Output("check_textarea_upregulated_lipids", "checked", allow_duplicate = True),
    Output("check_textarea_downregulated_lipids", "checked", allow_duplicate = True),
    Output("check_textarea_upregulated_proteins", "checked", allow_duplicate = True),
    Output("check_textarea_downregulated_proteins", "checked", allow_duplicate = True),
    Output("check_textarea_upregulated_metabolites", "checked", allow_duplicate = True),
    Output("check_textarea_downregulated_metabolites", "checked", allow_duplicate = True),
    Output("check_textarea_upregulated_transcripts", "checked", allow_duplicate = True),
    Output("check_textarea_downregulated_transcripts", "checked", allow_duplicate = True),
    Output("graph_enrichment_results", "rowData", allow_duplicate = True),
    Output("graph_enrichment_results", "selectedRows", allow_duplicate = True),
    Output("graph_enrichment_results", "filterModel", allow_duplicate = True),
    Input("button_reset", "n_clicks"),
    prevent_initial_call = True,
)
def reset_maptogo(n_clicks):
    if n_clicks == None:
        raise dash.exceptions.PreventUpdate

    return (
        False,
        "",
        False,
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "Entries: 0",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        False,
        False,
        False,
        False,
        False,
        INIT_ORGANISM,
        MOLECULE_HANDLING_REMOVE,
        MOLECULE_HANDLING_REMOVE,
        [],
        [],
        "greater",
        "fdr_bh",
        ["Biological process"],
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        True,
        [],
        [],
        {},
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
    State("check_textarea_upregulated_lipids", "checked"),
    State("check_textarea_downregulated_lipids", "checked"),
    State("check_textarea_upregulated_proteins", "checked"),
    State("check_textarea_downregulated_proteins", "checked"),
    State("check_textarea_upregulated_metabolites", "checked"),
    State("check_textarea_downregulated_metabolites", "checked"),
    State("check_textarea_upregulated_transcripts", "checked"),
    State("check_textarea_downregulated_transcripts", "checked"),
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
    check_textarea_upregulated_lipids,
    check_textarea_downregulated_lipids,
    check_textarea_upregulated_proteins,
    check_textarea_downregulated_proteins,
    check_textarea_upregulated_metabolites,
    check_textarea_downregulated_metabolites,
    check_textarea_upregulated_transcripts,
    check_textarea_downregulated_transcripts,
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
            upregulated_lipids_list if check_textarea_upregulated_lipids else [],
            downregulated_lipids_list if check_textarea_downregulated_lipids else [],
            all_proteins_list,
            regulated_proteins_list,
            upregulated_proteins_list if check_textarea_upregulated_proteins else [],
            downregulated_proteins_list if check_textarea_downregulated_proteins else [],
            all_metabolites_list,
            regulated_metabolites_list,
            upregulated_metabolites_list if check_textarea_upregulated_metabolites else [],
            downregulated_metabolites_list if check_textarea_downregulated_metabolites else [],
            all_transcripts_list,
            regulated_transcripts_list,
            upregulated_transcripts_list if check_textarea_upregulated_transcripts else [],
            downregulated_transcripts_list if check_textarea_downregulated_transcripts else [],
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
                "domain": " | ".join(ontology.get_domains(result.term.domains)),
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
                "domain": " | ".join(ontology.get_domains(result.term.domains)),
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
    Input("check_textarea_upregulated_lipids", "checked"),
    Input("check_textarea_downregulated_lipids", "checked"),
    Input("check_textarea_upregulated_proteins", "checked"),
    Input("check_textarea_downregulated_proteins", "checked"),
    Input("check_textarea_upregulated_metabolites", "checked"),
    Input("check_textarea_downregulated_metabolites", "checked"),
    Input("check_textarea_upregulated_transcripts", "checked"),
    Input("check_textarea_downregulated_transcripts", "checked"),
    Input("sessionid", "children"),
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
    check_textarea_upregulated_lipids,
    check_textarea_downregulated_lipids,
    check_textarea_upregulated_proteins,
    check_textarea_downregulated_proteins,
    check_textarea_upregulated_metabolites,
    check_textarea_downregulated_metabolites,
    check_textarea_upregulated_transcripts,
    check_textarea_downregulated_transcripts,
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
    session.ui["check_textarea_upregulated_lipids"] = check_textarea_upregulated_lipids
    session.ui["check_textarea_downregulated_lipids"] = check_textarea_downregulated_lipids
    session.ui["check_textarea_upregulated_proteins"] = check_textarea_upregulated_proteins
    session.ui["check_textarea_downregulated_proteins"] = check_textarea_downregulated_proteins
    session.ui["check_textarea_upregulated_metabolites"] = check_textarea_upregulated_metabolites
    session.ui["check_textarea_downregulated_metabolites"] = check_textarea_downregulated_metabolites
    session.ui["check_textarea_upregulated_transcripts"] = check_textarea_upregulated_transcripts
    session.ui["check_textarea_downregulated_transcripts"] = check_textarea_downregulated_transcripts

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
    State("sunburst_controls", "style"),
    prevent_initial_call = True,
)
def open_sankeyplot(
    n_clicks,
    radiogroup_sankey,
    row_data,
    selected_rows,
    session_id,
    sunburst_controls_style,
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
                            if term.categories > 0:
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
                if type(layer_categories) == int:
                    layer_categories = [i for i in range(layer_categories.bit_length()) if (layer_categories >> i) & 1]
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
                    layer_categories_prev = path_layers[i - 1][1]
                    if type(layer_categories_prev) == int:
                        layer_categories_prev = [i for i in range(layer_categories_prev.bit_length()) if (layer_categories_prev >> i) & 1]
                    for layer_category_prev in layer_categories_prev:
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

    return (
        True,
        fig,
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
    Output("barplot_numberinput_connect_ths", "disabled", allow_duplicate = True),
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
            no_update,
        )

    session = sessions[session_id]

    ontology = session.ontology
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

    # distance matrix
    jaccard_values = np.zeros((n , n))
    for i in range(0, n - 1):
        i_terms = source_terms[i]
        for j in range(i + 1, n):
            j_terms = source_terms[j]
            jaccard_values[j, i] = jaccard_values[i, j] = 1 - len(i_terms & j_terms) / len(i_terms | j_terms)


    Z, sort_order = None, None
    if bar_sorting == BAR_SORTING_SIMILARITY and n > 1:
        if jaccard_values.var() < 1e-20:
            jaccard_values += (r := np.random.uniform(0, 1e-5, (n, n))) * r.T * (1 - np.diag([1] * n))

        try:
            Z = linkage(squareform(jaccard_values), method = 'average')
            sort_order = leaves_list(Z)
            selected_term_ids = selected_term_ids[sort_order]
            jaccard_values = jaccard_values[sort_order]
            jaccard_values = jaccard_values[:, sort_order]
        except Exception as e:
            pass

    if session.separate_updown_switch:
        pvalues = -np.log10([session_data[term_id].pvalue_corrected[0] for term_id in selected_term_ids])
        number_regulated_entities = np.array([session_data[term_id].fisher_data[0][0] + session_data[term_id].fisher_data[1][0] for term_id in selected_term_ids])
    else:
        pvalues = -np.log10([session_data[term_id].pvalue_corrected for term_id in selected_term_ids])
        number_regulated_entities = np.array([session_data[term_id].fisher_data[0] for term_id in selected_term_ids])
    number_entities = np.array([len(session_data[term_id].source_terms) for term_id in selected_term_ids])

    domains = ["|".join(ontology.get_domains(session_data[term_id].term.domains)) for term_id in selected_term_ids]
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
    sorted_domains = sorted(list(domain_colors.keys()))

    max_domain_num = 0
    max_categories_num = 0
    remaining_domains_categories = [[0, {dom: 0 for dom in domain_colors.keys()}, {}] for i in range(n)]
    domain_bar = bar_width * 0.95 / len(sorted_domains)
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
            domain_colors[ontology.get_domains(result.term.domains)[0]][0],
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
                domain_colors[ontology.get_domains(result.term.domains)[0]][1],
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

            for domain in ontology.get_domains(foreign_term.domains):
                remaining_domains_categories[i][1][domain] += 1
                max_domain_num = max(max_domain_num, remaining_domains_categories[i][1][domain])
        remaining_domains_categories[i][0] = description_start_angle - bar_width * 0.025

        for molecule in molecules:
            for term in get_path(session.all_parent_nodes[molecule], result.term):
                if type(term) == OntologyTerm: break

            for category in [ontology.categories[i] for i in range(term.categories.bit_length()) if (term.categories >> i) & 1]:
                if category not in remaining_domains_categories[i][2]:
                    remaining_domains_categories[i][2][category] = 0
                remaining_domains_categories[i][2][category] += 1
                max_categories_num = max(max_categories_num, remaining_domains_categories[i][2][category])

    for i in range(n):
        domains_position = remaining_domains_categories[i][0]
        for domain_name in sorted_domains:
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

    if bar_sorting == BAR_SORTING_SIMILARITY and n > 1 and type(Z) != NoneType:
        add_arc(fig, 0, 360, 0, outer_inner_radius * 0.05, "#000000")

        dendrogram_points = np.zeros(((n - 1) * 4, 2))
        ddata = dendrogram(Z, no_plot = True)
        icoord = ddata['icoord']  # x-coordinates for each link
        dcoord = ddata['dcoord']  # y-coordinates (heights) for each link

        i = 0
        for xs, ys in zip(icoord, dcoord):
            for x, y in zip(xs, ys):
                dendrogram_points[i, 0] = x
                dendrogram_points[i, 1] = outer_inner_radius * (1 - y)
                i += 1

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

        draw_line(fig, (dendrogram_points[-2, 0] + dendrogram_points[-3, 0]) / 2., dendrogram_points[-2, 1], 0)

    elif bar_sorting == BAR_SORTING_PVALUE:
        jaccard_max_saturation = 90
        for i in range(0, n - 1):
            for j in range(i + 1, n):
                jaccard = 1 - jaccard_values[i, j]
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
    barplot_terms_wrapper_style["height"] = "70%"

    return (
        True,
        fig,
        CIRCLE_WIDTH * 2 + 450,
        barplot_terms_wrapper_style,
        no_update,
        no_update,
        barplot_controls_style,
        sunburst_controls_style,
        bar_sorting == BAR_SORTING_SIMILARITY,
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
        example = api_background_lipids_example,

    ),
    "regulated_lipids": fields.List(
        fields.String,
        description = "All regulated lipid names in experiment",
        example = api_regulated_lipids_example,
    ),
    "upregulated_lipids": fields.List(
        fields.String,
        description = "All up-regulated lipid names in experiment",
        example = [],
    ),
    "downregulated_lipids": fields.List(
        fields.String,
        description = "All down-regulated lipid names in experiment",
        example = [],
    ),
    "background_proteins": fields.List(
        fields.String,
        description = "All protein accessions in experiment (background)",
        example = api_background_proteins_example,
    ),
    "regulated_proteins": fields.List(
        fields.String,
        description = "All regulated protein accession in experiment",
        example = api_regulated_proteins_example,
    ),
    "upregulated_proteins": fields.List(
        fields.String,
        description = "All up-regulated protein accession in experiment",
        example = [],
    ),
    "downregulated_proteins": fields.List(
        fields.String,
        description = "All down-regulated protein accession in experiment",
        example = [],
    ),
    "background_metabolites": fields.List(
        fields.String,
        description = "All metabolite ChEBI Ids in experiment (background)",
        example = api_background_metabolites_example,
    ),
    "regulated_metabolites": fields.List(
        fields.String,
        description = "All regulated metabolite ChEBI Ids in experiment",
        example = api_regulated_metabolites_example,
    ),
    "upregulated_metabolites": fields.List(
        fields.String,
        description = "All up-regulated metabolite ChEBI Ids in experiment",
        example = [],
    ),
    "downregulated_metabolites": fields.List(
        fields.String,
        description = "All down-regulated metabolite ChEBI Ids in experiment",
        example = [],
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
    "upregulated_transcripts": fields.List(
        fields.String,
        description = "All up-regulated ensembl Ids in experiment",
        example = [],
    ),
    "downregulated_transcripts": fields.List(
        fields.String,
        description = "All down-regulated ensembl Ids in experiment",
        example = [],
    ),
    "separate_updown": fields.Boolean(
        description = f"Or regulated molecules are available (False) or up / down-regulated (True), default: false",
        example = False,
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
            upregulated_lipids = data.get("upregulated_lipids", [])
            downregulated_lipids = data.get("downregulated_lipids", [])
            background_proteins = data.get("background_proteins", [])
            regulated_proteins = data.get("regulated_proteins", [])
            upregulated_proteins = data.get("upregulated_proteins", [])
            downregulated_proteins = data.get("downregulated_proteins", [])
            background_metabolites = data.get("background_metabolites", [])
            regulated_metabolites = data.get("regulated_metabolites", [])
            upregulated_metabolites = data.get("upregulated_metabolites", [])
            downregulated_metabolites = data.get("downregulated_metabolites", [])
            background_transcripts = data.get("background_transcripts", [])
            regulated_transcripts = data.get("regulated_transcripts", [])
            upregulated_transcripts = data.get("upregulated_transcripts", [])
            downregulated_transcripts = data.get("downregulated_transcripts", [])

            if type(background_lipids) != list:
                return {"error_message": "'background_lipids' needs to be a list", "result": []}, 422
            if type(regulated_lipids) != list:
                return {"error_message": "'regulated_lipids' needs to be a list", "result": []}, 422
            if type(upregulated_lipids) != list:
                return {"error_message": "'upregulated_lipids' needs to be a list", "result": []}, 422
            if type(downregulated_lipids) != list:
                return {"error_message": "'downregulated_lipids' needs to be a list", "result": []}, 422
            if type(background_proteins) != list:
                return {"error_message": "'background_proteins' needs to be a list", "result": []}, 422
            if type(regulated_proteins) != list:
                return {"error_message": "'regulated_proteins' needs to be a list", "result": []}, 422
            if type(upregulated_proteins) != list:
                return {"error_message": "'upregulated_proteins' needs to be a list", "result": []}, 422
            if type(downregulated_proteins) != list:
                return {"error_message": "'downregulated_proteins' needs to be a list", "result": []}, 422
            if type(background_metabolites) != list:
                return {"error_message": "'background_metabolites' needs to be a list", "result": []}, 422
            if type(regulated_metabolites) != list:
                return {"error_message": "'regulated_metabolites' needs to be a list", "result": []}, 422
            if type(upregulated_metabolites) != list:
                return {"error_message": "'upregulated_metabolites' needs to be a list", "result": []}, 422
            if type(downregulated_metabolites) != list:
                return {"error_message": "'downregulated_metabolites' needs to be a list", "result": []}, 422
            if type(background_transcripts) != list:
                return {"error_message": "'background_transcripts' needs to be a list", "result": []}, 422
            if type(regulated_transcripts) != list:
                return {"error_message": "'regulated_transcripts' needs to be a list", "result": []}, 422
            if type(upregulated_transcripts) != list:
                return {"error_message": "'upregulated_transcripts' needs to be a list", "result": []}, 422
            if type(downregulated_transcripts) != list:
                return {"error_message": "'downregulated_transcripts' needs to be a list", "result": []}, 422

            organism_api = data.get("organism_taxonomy", "NCBITaxon:9606")
            domains_api = data.get("domains", ["biological_process"])
            accepted_domains = {d.lower().replace(" ", "_"): d for d in enrichment_ontologies[INIT_ORGANISM].domains}
            pvalue_correction_api = data.get("pvalue_correction", "fdr_bh")
            term_representation_api = data.get("term_representation", "greater")
            unrecognizable_molecules_api = data.get("unrecognizable_molecules", MOLECULE_HANDLING_IGNORE)
            non_background_molecules_api = data.get("non_background_molecules", MOLECULE_HANDLING_REMOVE)
            bounded_fatty_acyls_api = data.get("bounded_fatty_acyls", False)
            separate_updown_api = data.get("separate_updown", False)

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

            if type(separate_updown_api) not in {int, bool}:
                return {"error_message": "'separate_updown' must be a boolean (0 | 1)", "result": []}, 422

            if type(separate_updown_api) == int:
                v_api = separate_updown_api > 0

            with_lipids = len(background_lipids) > 0
            with_proteins = len(background_proteins) > 0
            with_metabolites = len(background_metabolites) > 0
            with_transcripts = len(background_transcripts) > 0

            omics_included = [with_lipids, with_proteins, with_metabolites, with_transcripts]

            omics_lists = [
                background_lipids,
                regulated_lipids,
                upregulated_lipids,
                downregulated_lipids,
                background_proteins,
                regulated_proteins,
                upregulated_proteins,
                downregulated_proteins,
                background_metabolites,
                regulated_metabolites,
                upregulated_metabolites,
                downregulated_metabolites,
                background_transcripts,
                regulated_transcripts,
                upregulated_transcripts,
                downregulated_transcripts,
            ]

            ontology = enrichment_ontologies[organism_api]

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
                    separate_updown_api,
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
            session.separate_updown_switch = separate_updown_api
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
            results = ontology.enrichment_analysis(session.separate_updown_switch, session.search_terms, session.num_background, target_set, domains_api, term_representation_api, pvalue_correction_api)
            session.result = results

            result_data = []
            for result in results:
                expected = round((result.fisher_data[0] + result.fisher_data[1]) * (result.fisher_data[0] + result.fisher_data[2]) / sum(result.fisher_data))
                row = {
                    "domain": " | ".join(ontology.get_domains(result.term.domains)),
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



