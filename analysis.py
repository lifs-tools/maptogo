import logging
import os

logger = logging.getLogger(__name__)
logging.basicConfig(
    level = os.environ.get("LOGLEVEL", "INFO"),
    format = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
)
logger.info("Started enrichment server")


from dash import Dash, dcc, html, Input, Output, State, callback, exceptions, no_update, MATCH, ALL, callback_context, clientside_callback, dash_table, ctx
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
from statsmodels.stats.multitest import multipletests
from EnrichmentDataStructure import EnrichmentOntology, current_path, SessionEntry, OntologyTerm, TermType
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidAdduct import LipidAdduct
import pathlib
import time
import hashlib
import threading
import freetype
from math import ceil
from flask_restx import Api, Resource, fields
import traceback
import requests
import threading


INIT_ORGANISM = "10090"

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
    from organisms import organisms
except Exception as e:
    organisms = {
        #'Homo sapiens': '9606',
        'Mus musculus': '10090',
        # 'Bacillus cereus': "405534",
        # 'Saccharomyces cerevisiae': '4932',
        # 'Escherichia coli': '562',
        # 'Drosophila melanogaster': '7227',
        'Rattus norvegicus': '10116',
        # 'Bos taurus': '9913',
        # 'Caenorhabditis elegans': '6239',
        # 'Pseudomonas aeruginosa': '287',
        # 'Arabidopsis thaliana': '3702',
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



hash_function = hashlib.new('sha256')
LINK_COLOR = "#2980B9"
SESSION_DURATION_TIME = 60 * 60 * 2 # two hours
domain_colors = {
    "Biological process": ["#F09EA7", "#e5a9b0"],
    "Cellular component": ["#F6CA94", "#eac9a0"],
    "Disease": ["#FAFABE", "#f3f3c5"],
    "Metabolic and signalling pathway": ["#C1EBC0", "#c9e3c8"],
    "Molecular function": ["#C7CAFF", "#cdcff9"],
    "Phenotype": ["#CDABEB", "#ccb5e1"],
    "Physical or chemical properties": ["#F6C2FA", "#f0c9f3"],
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

MOLECULE_HANDLING_ERROR = "molecule_handling_error"
MOLECULE_HANDLING_REMOVE = "molecule_handling_remove"
MOLECULE_HANDLING_IGNORE = "molecule_handling_ignore"



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



def analytics(action):
    return

    def send_action(recorded_action):
        try:
            url = "https://lifs-tools.org/matomo/matomo.php?idsite=17&rec=1&e_c=MOEA-1.0.0&e_a=" + recorded_action
            response = requests.get(url, timeout = 5)
        except Exception as e:
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
    if worksheet_name[0] == "_": continue
    df = xl.parse(worksheet_name)
    worksheet = {
        "bgl": [], "regl": [],
        "bgp": [], "regp": [],
        "bgm": [], "regm": [],
        "bgt": [], "regt": [],
        "title": "", "desc": "", "comp": "", "doi": "", "org": INIT_ORGANISM
    }
    if "BackgroundLipids" in df: worksheet["bgl"] = list(df["BackgroundLipids"].dropna())
    if "RegulatedLipids" in df: worksheet["regl"] = list(df["RegulatedLipids"].dropna())
    if "BackgroundProteins" in df: worksheet["bgp"] = list(df["BackgroundProteins"].dropna())
    if "RegulatedProteins" in df: worksheet["regp"] = list(df["RegulatedProteins"].dropna())
    if "BackgroundMetabolites" in df: worksheet["bgm"] = list(df["BackgroundMetabolites"].dropna())
    if "RegulatedMetabolites" in df: worksheet["regm"] = list(df["RegulatedMetabolites"].dropna())
    if "BackgroundTranscripts" in df: worksheet["bgt"] = list(df["BackgroundTranscripts"].dropna())
    if "RegulatedTranscripts" in df: worksheet["regt"] = list(df["RegulatedTranscripts"].dropna())
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
app.title = "GO multiomics"

api = Api(
    app.server,
    version = "1.0",
    title = "Gene Ontology (GO) multiomics enrichment analysis - REST API",
    description = "Application programming interface for the GO multiomics enrichment analysis",
    doc = "/api/docs"
)

lipid_parser = LipidParser()
enrichment_ontologies = {}
for tax_name, tax_id in organisms.items():
    logger.info(f"Loading {tax_name}")
    enrichment_ontologies[tax_id] = EnrichmentOntology(f"{current_path}/Data/ontology_{tax_id}.gz", tax_name, lipid_parser = lipid_parser)


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
    current_time = time.time()
    hash_function.update(f"{current_time}t".encode())
    session_id = hash_function.hexdigest()
    sessions[session_id] = SessionEntry()
    logger.info(f"New session: {session_id}")

    ontology = enrichment_ontologies[INIT_ORGANISM]
    predefined_proteins = [
        len(ontology.reviewed_proteins),
        len(set(ontology.proteins.keys()) - ontology.reviewed_proteins),
        len(ontology.proteins),
    ]

    return html.Div([
        dcc.Download(id = "download_data"),
        dmc.TextInput(id = "html_download_trigger", style = {"display": "none"}),
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
                    dmc.Title("Gene Ontology (GO) multiomics enrichment analysis"),
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
                ),
                dmc.Text("|", style = {"marginLeft": "20px", "marginRight": "20px"}),
                html.A(
                    dmc.Text(
                        "REST API",
                    ),
                    id = "api_link",
                    href = "/api/docs",
                    target = "_blank",
                    style = {"color": LINK_COLOR, "cursor": "pointer"},
                ),
                dmc.Text(" ", style = {"marginLeft": "5px"}),
                ], style = {"height": "100%", "display": "flex", "alignItems": "flex-start", "justifyContent": "right"},
            ),
        ], cols = 2),
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
            #title = "Load",
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
                            dmc.Text("Maximum p-value", size = "xs", w = "500"),
                            dmc.Text("Minimum p-value", size = "xs", w = "500"),
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
                html.Div(
                    dmc.Switch(
                        id = "switch_sankey_regulated_only",
                        checked = False,
                        label = "Regulated molecules only",
                        style = {"paddingBottom": "8px", "paddingRight": "2px"},
                    ),
                    id = "sankey_controls",
                    style = {"marginTop": "10px"},
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
                dmc.Title("About GO multiomics", order = 4),
                dmc.Text(
                    "GO multiomics is the first implementation that enables GO term enrichment analysis on lipids. Starting with a list of identified lipids (at least on the species level) as background and a list with only the (differentially) regulated lipids, each lipid is mapped to proteins involved in metabolic reactions. On this base, GO term analysis will then be performed with a statistical test. Depending on the selected domain (biological process, molecular function, cellular compartment, physical or chemical properties, or metabolic and signaling pathways), the analysis provides a sorted list of GO terms with a p-value for each term.",
                    style = {"textAlign": "justify"},
                ),
            ]),
            dmc.Text("Responsible people for this database are:"),
            dmc.Text("Dominik Kopczynski: Implementation of both front- and backend"),
            dmc.Text("Cristina Coman: Method validation"),
            dmc.Text("Robert Ahrends: Project leader"),

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
                            ], cols = 2),
                            dmc.SimpleGrid([
                                dmc.Group(
                                    dmc.Text(
                                        "Entries: 0",
                                        id = "num_all_lipids",
                                        style = {"color": "#808080"},
                                        size = "12px",
                                    ),
                                    position = "right",
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
                            ], cols = 2),
                            dmc.SimpleGrid([
                                dmc.Group(
                                    dmc.Text(
                                        "Entries: 0",
                                        id = "num_all_proteins",
                                        style = {"color": "#808080"},
                                        size = "12px",
                                    ),
                                    position = "right",
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
                            ], cols = 2),
                            dmc.SimpleGrid([
                                dmc.Group(
                                    dmc.Text(
                                        "Entries: 0",
                                        id = "num_all_metabolites",
                                        style = {"color": "#808080"},
                                        size = "12px",
                                    ),
                                    position = "right",
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
                            ], cols = 2)],
                            value = "metabolites_tab",
                        ),
                        dmc.TabsPanel([
                            dmc.SimpleGrid([
                                dmc.Title(
                                    "All ensembl Ids in experiment (background)",
                                    order = 5,
                                    style = {"marginTop": "10px"},
                                ),
                                dmc.Title(
                                    "All regulated ensembl Ids in experiment",
                                    order = 5,
                                    style = {"marginTop": "10px"},
                                ),
                            ], cols = 2),
                            dmc.SimpleGrid([
                                dmc.Textarea(
                                    id = "textarea_all_transcripts",
                                    style = {"height": "100%", "display": "inline"},
                                    minRows = 15,
                                ),
                                dmc.Textarea(
                                    id = "textarea_regulated_transcripts",
                                    minRows = 15,
                                ),
                            ], cols = 2),
                            dmc.SimpleGrid([
                                dmc.Group(
                                    dmc.Text(
                                        "Entries: 0",
                                        id = "num_all_transcripts",
                                        style = {"color": "#808080"},
                                        size = "12px",
                                    ),
                                    position = "right",
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
                                dmc.Select(
                                    id = "select_term_representation",
                                    data = term_representation,
                                    value = "greater",
                                    label = "Term representation:",
                                ),
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
                                                    dmc.Text("When activated, fatty chains in lipids such as the 20:4 (AA) in 'PI 18:0/20:4' will be considered for the analysis."),
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
                html.Div(
                    dag.AgGrid(
                        id = "graph_enrichment_results",
                        columnDefs = [
                            {
                                'field': "domain",
                                "headerName": "Domain",
                                "maxWidth": 200,
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
                                'field': "count",
                                "headerName": "Count",
                                "width": 80,
                                "headerTooltip": "Number of regulated associated molecules (Expected number of regulated associated molecules) / Number of all associated molecules",
                            },
                            {
                                'field': "pvalue",
                                "headerName": "p-value",
                                "width": 100,
                                "valueFormatter": {"function": "params.value != null ? Number(params.value).toPrecision(6) : ''"},
                                "headerTooltip": "The p-value is the statistical significance, the q-value is the adjusted p-value after multiple testing correction",
                            },
                            {
                                'field': "log_odds_ratio",
                                "headerName": "Log Odds ratio",
                                "width": 100,
                                "valueFormatter": {"function": "params.value != null ? Number(params.value).toPrecision(6) : ''"},
                                "headerTooltip": "The log odds ratio determines how strong the association is (effect size). Positive: enriched; Zero: no enrichment; Negative: depleted.",
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
                        className = "ag-theme-alpine terms-theme",
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
                    className = "terms-overlay",
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



def check_user_input(
    omics_included,
    omics_lists,
    ontology,
    ignore_unrecognizable_molecules,
    ignore_unknown,
):
    with_lipids, with_proteins, with_metabolites, with_transcripts = omics_included
    target_set = set()
    lipidome, regulated_lipids = {}, set()
    proteome, regulated_proteins = set(), set()
    metabolome, regulated_metabolites = set(), set()
    transcriptome, regulated_transcripts = set(), set()
    background_list = []

    (
        all_lipids_list,
        regulated_lipids_list,
        all_proteins_list,
        regulated_proteins_list,
        all_metabolites_list,
        regulated_metabolites_list,
        all_transcripts_list,
        regulated_transcripts_list,
    ) = omics_lists

    if with_lipids:
        if type(all_lipids_list) == str: all_lipids_list = all_lipids_list.split("\n")
        elif len(all_lipids_list) == 0:
            return True, "Please paste lipid names into the first text area.", []
        if type(regulated_lipids_list) == str: regulated_lipids_list = regulated_lipids_list.split("\n")
        elif len(regulated_lipids_list) == 0:
            return True, "Please paste lipid names into the second text area.", []

        for lipid_name in all_lipids_list:
            if len(lipid_name) == 0: continue
            try:
                lipid = lipid_parser.parse(lipid_name)
            except Exception as e:
                lower_lipid_name = lipid_name.lower()
                if lower_lipid_name in ontology.unspecific_lipids:
                    lipid = ontology.unspecific_lipids[lower_lipid_name]
                else:
                    if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                        lipid = None
                    elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                        continue
                    else:
                        return True, f"Lipid name '{lipid_name}' unrecognizable! Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.", []
            lipidome[lipid_name] = lipid

        for lipid_name in regulated_lipids_list:
            if len(lipid_name) == 0: continue
            if lipid_name not in lipidome:
                if ignore_unknown == MOLECULE_HANDLING_REMOVE: continue
                return True, f"The regulated lipid '{lipid_name}' does not occur in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of non-background regulated molecules' setting.", []
            regulated_lipids.add(lipid_name)

        if len(lipidome) == 0:
            return True, "No background lipid left after lipid recognition.", []

        if len(regulated_lipids) == 0:
            return True, "No regulated lipid left after lipid recognition.", []

        if len(regulated_lipids) > len(lipidome):
            return True, "Length of regulated lipid list must be smaller than background list.", []

        left_lipids = regulated_lipids - lipidome.keys()
        if len(left_lipids) > 0:
            if ignore_unknown == MOLECULE_HANDLING_REMOVE:
                for lipid_name in left_lipids:
                    del lipidome[lipid_name]
            else:
                return True, "The regulated lipid" + (' ' if len(left_lipids) == 1 else 's ') + "'" + "', '".join(left_lipids) + ("' does" if len(left_lipids) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of non-background regulated molecules' setting.", []

        target_set |= regulated_lipids
        background_list += [{"value": k, "label": k} for k in lipidome.keys()]

    if with_proteins:
        if type(all_proteins_list) == str: all_proteins_list = all_proteins_list.split("\n")
        elif len(all_proteins_list) == 0:
            return True, "Please paste protein accession into the first text area.", []

        if type(regulated_proteins_list) == str: regulated_proteins_list = regulated_proteins_list.split("\n")
        elif len(regulated_proteins_list) == 0:
            return True, "Please paste protein accessions into the second text area.", []

        proteome = set(protein for protein in all_proteins_list if len(protein) > 0)
        regulated_proteins = set(protein for protein in regulated_proteins_list if len(protein) > 0)

        background_list += [{"value": pp, "label": p + (' (' + ontology.proteins[pp].name + ')' if (pp in ontology.proteins) else '')} for p in proteome if (pp := "UNIPROT:" + p)]
        proteome = set(p.split("-")[0] for p in proteome)
        regulated_proteins = set(rp.split("-")[0] for rp in regulated_proteins)
        left_proteins = proteome - ontology.clean_protein_ids
        if len(left_proteins) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                proteome -= left_proteins
            else:
                return True, "The protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' is" if len(left_proteins) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.", []

        left_proteins = regulated_proteins - ontology.clean_protein_ids
        if len(left_proteins) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                regulated_proteins -= left_proteins
            else:
                return True, "The protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' is" if len(left_proteins) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.", []

        left_proteins = regulated_proteins - proteome
        if len(left_proteins) > 0:
            if ignore_unknown == MOLECULE_HANDLING_REMOVE:
                regulated_proteins -= left_proteins
            else:
                return True, "The regulated protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' does" if len(left_proteins) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of non-background regulated molecules' setting.", []

        if len(proteome) == 0:
            return True, "No background protein left after protein recognition.", []

        if len(regulated_proteins) == 0:
            return True, "No regulated protein left after protein recognition.", []

        if len(regulated_proteins) > len(proteome):
            return True, "Length of regulated protein list must be smaller than background list.", []

        proteome = set([f"UNIPROT:{protein}" for protein in proteome])
        regulated_proteins = set([f"UNIPROT:{protein}" for protein in regulated_proteins])
        target_set |= regulated_proteins

    if with_metabolites:
        if type(all_metabolites_list) == str: all_metabolites_list = all_metabolites_list.split("\n")
        elif len(all_metabolites_list) == 0:
            return True, "Please paste metabolite ChEBI Ids into the first text area.", []

        if type(regulated_metabolites_list) == str: regulated_metabolites_list = regulated_metabolites_list.split("\n")
        elif len(regulated_metabolites_list) == 0:
            return True, "Please paste metabolite ChEBI Ids into the second text area.", []

        metabolome = set(metabolite for metabolite in all_metabolites_list if len(metabolite) > 0)
        regulated_metabolites = set(metabolite for metabolite in regulated_metabolites_list if len(metabolite) > 0)

        background_list += [{"value": m, "label": m} for m in metabolome]
        left_metabolites = metabolome - ontology.clean_metabolite_ids - ontology.metabolites.keys()
        left_metabolites -= set([m for m in left_metabolites if m.lower() in ontology.metabolite_names.keys()])
        if len(left_metabolites) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                metabolome -= left_metabolites
            else:
                return True, "The metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' is" if len(left_metabolites) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.", []

        left_metabolites = regulated_metabolites - ontology.clean_metabolite_ids - ontology.metabolites.keys()
        left_metabolites -= set([m for m in left_metabolites if m.lower() in ontology.metabolite_names.keys()])
        if len(left_metabolites) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                regulated_metabolites -= left_metabolites
            else:
                return True, "The metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' is" if len(left_metabolites) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.", []

        left_metabolites = regulated_metabolites - metabolome
        if len(left_metabolites) > 0:
            if ignore_unknown == MOLECULE_HANDLING_REMOVE:
                regulated_metabolites -= left_metabolites
            else:
                return True, "The regulated metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' does" if len(left_metabolites) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", []

        if len(metabolome) == 0:
            return True, "No background metabolite left after metabolite recognition.", []

        if len(regulated_metabolites) == 0:
            return True, "No regulated metabolite left after metabolite recognition.", []

        if len(regulated_metabolites) > len(metabolome):
            return True, "Length of regulated metabolite list must be smaller than background list.", []

        metabolome = set([f"CHEBI:{metabolite}" if (type(metabolite) == int or not metabolite.startswith("CHEBI:")) and metabolite.lower() not in ontology.metabolite_names else metabolite for metabolite in metabolome])
        regulated_metabolites = set([f"CHEBI:{metabolite}" if (type(metabolite) == int or not metabolite.startswith("CHEBI:")) and metabolite.lower() not in ontology.metabolite_names else metabolite for metabolite in regulated_metabolites])
        target_set |= regulated_metabolites

    if with_transcripts:
        if type(all_transcripts_list) == str: all_transcripts_list = all_transcripts_list.split("\n")
        elif len(all_transcripts_list) == 0:
            return True, "Please paste transcript ensembl Ids into the first text area.", []

        if type(regulated_transcripts_list) == str: regulated_transcripts_list = regulated_transcripts_list.split("\n")
        elif len(regulated_transcripts_list) == 0:
            return True, "Please paste transcript ensembl Ids into the second text area.", []

        transcriptome = set(transcript for transcript in all_transcripts_list if len(transcript) > 0)
        regulated_transcripts = set(transcript for transcript in regulated_transcripts_list if len(transcript) > 0)

        background_list += [{"value": t, "label": t + (' (' + ontology.transcripts[tt].name + ')' if (tt in ontology.transcripts) else '')} for t in transcriptome if (tt := t.split(".")[0])]
        transcript_keys = set(ontology.transcripts.keys())
        left_transcripts = set(t for t in transcriptome if t.split(".")[0] not in transcript_keys)
        if len(left_transcripts) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                transcriptome -= left_transcripts
            else:
                return True, "The transcript" + (' ' if len(left_transcripts) == 1 else 's ') + "'" + "', '".join(left_transcripts) + ("' is" if len(left_transcripts) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.", []

        left_transcripts = set(t for t in regulated_transcripts if t.split(".")[0] not in transcript_keys)
        if len(left_transcripts) > 0:
            if ignore_unrecognizable_molecules == MOLECULE_HANDLING_IGNORE:
                pass
            elif ignore_unrecognizable_molecules == MOLECULE_HANDLING_REMOVE:
                regulated_transcripts -= left_transcripts
            else:
                return True, "The transcript" + (' ' if len(left_transcripts) == 1 else 's ') + "'" + "', '".join(left_transcripts) + ("' is" if len(left_transcripts) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Remove for analysis' option in the 'Handling of unrecognizable molecules' setting.", []

        left_transcripts = regulated_transcripts - transcriptome
        if len(left_transcripts) > 0:
            if ignore_unknown == MOLECULE_HANDLING_REMOVE:
                regulated_transcripts -= left_transcripts
            else:
                return True, "The regulated transcript" + (' ' if len(left_transcripts) == 1 else 's ') + "'" + "', '".join(left_transcripts) + ("' does" if len(left_transcripts) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", []
        if len(transcriptome) == 0:
            return True, "No background transcript left after transcript recognition.", []

        if len(regulated_transcripts) == 0:
            return True, "No regulated transcript left after transcript recognition.", []

        if len(regulated_transcripts) > len(transcriptome):
            return True, "Length of regulated transcript list must be smaller than background list.", []

        target_set |= regulated_transcripts

    molecule_tables = [
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
    ]

    return False, "", molecule_tables






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
    Input("button_run_enrichment", "n_clicks"),
    State("textarea_all_lipids", "value"),
    State("textarea_regulated_lipids", "value"),
    State("textarea_all_proteins", "value"),
    State("textarea_regulated_proteins", "value"),
    State("textarea_all_metabolites", "value"),
    State("textarea_regulated_metabolites", "value"),
    State("textarea_all_transcripts", "value"),
    State("textarea_regulated_transcripts", "value"),
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
    all_transcripts_list,
    regulated_transcripts_list,
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
):
    histogram_disabled = True

    if session_id not in sessions:
        return "", [], [], {}, True, "Your session has expired. Please refresh the website.",  histogram_disabled, [], []

    logger.info(f"Enrichment session: {session_id}")
    session = sessions[session_id]
    session.time = time.time()

    if not with_lipids and not with_proteins and not with_metabolites and not with_transcripts:
        return "", [], [], {}, True, "No omics data is selected.", histogram_disabled, [], []

    if len(domains) == 0:
        return "", [], [], {}, True, "No domain(s) selected.", histogram_disabled, [], []

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
            all_proteins_list,
            regulated_proteins_list,
            all_metabolites_list,
            regulated_metabolites_list,
            all_transcripts_list,
            regulated_transcripts_list,
        ]

        do_activate_alert, error_message, molecule_tables = check_user_input(
            omics_included,
            omics_lists,
            ontology,
            ignore_unrecognizable_molecules,
            ignore_unknown,
        )

        if do_activate_alert:
            return "", [], [], {}, do_activate_alert, error_message, histogram_disabled, [], []

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
        ) = molecule_tables

        ontology.set_background(session, lipid_dict = lipidome, protein_set = proteome, metabolite_set = metabolome, transcript_set = transcriptome)
        session.ontology = ontology
        session.data_loaded = True

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

    else:
        target_set = set()
        if with_lipids: target_set |= session.regulated_lipids
        if with_proteins: target_set |= session.regulated_proteins
        if with_metabolites: target_set |= session.regulated_metabolites
        if with_transcripts: target_set |= session.regulated_transcripts
        background_list = no_update

    analytics("enrichment_analysis")
    session.domains = set(domains)
    results = ontology.enrichment_analysis(session, target_set, domains, term_regulation)
    session.result = results

    if correction_method != "no" and len(results) > 1:
        pvalues = [r.pvalue for r in results]
        pvalues = multipletests(pvalues, method = correction_method)[1]
        for pvalue, r in zip(pvalues, results): r.pvalue_corrected = pvalue

    results.sort(key = lambda row: (row.pvalue_corrected, row.term.name))

    data = []
    session.data = {}
    session.results = []
    result_map = {}
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
        data.append(row)
        session.results.append((result, row))

    histogram_disabled = False
    for result in results:
        for term_id in result.term.term_id:
            session.data[term_id] = result


    return (
        "",
        data,
        [],
        {},
        False,
        "",
        histogram_disabled,
        background_list,
        [],
    )



@callback(
    Output("graph_enrichment_results", "rowData", allow_duplicate = True),
    Output("graph_enrichment_results", "selectedRows", allow_duplicate = True),
    Output("graph_enrichment_results", "filterModel", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Input("multiselect_filter_molecules", "value"),
    State("multiselect_filter_molecules", "data"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def filter_result_table(multiselect_values, multiselect_data, session_id):
    if multiselect_data == None or len(multiselect_data) == 0 or multiselect_values == None or session_id == None:
        raise exceptions.PreventUpdate

    if session_id not in sessions:
        return (
            no_update,
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
        )

    session = sessions[session_id]
    session.time = time.time()

    if session.results is None or len(session.results) == 0:
        raise exceptions.PreventUpdate

    if len(multiselect_values) > 0:
        multiselect_values = set(multiselect_values)
        for r, d in session.results:
            print(r.source_terms)
        data = [row for result, row in session.results if not set(result.source_terms).isdisjoint(multiselect_values)]
    else:
        data = [row for _, row in session.results]
    return data, [], {}, False, ""



@callback(
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("num_all_lipids", "children", allow_duplicate = True),
    Output("num_regulated_lipids", "children", allow_duplicate = True),
    Output("num_all_proteins", "children", allow_duplicate = True),
    Output("num_regulated_proteins", "children", allow_duplicate = True),
    Output("num_all_metabolites", "children", allow_duplicate = True),
    Output("num_regulated_metabolites", "children", allow_duplicate = True),
    Output("num_all_transcripts", "children", allow_duplicate = True),
    Output("num_regulated_transcripts", "children", allow_duplicate = True),
    Input("textarea_all_lipids", "value"),
    Input("textarea_regulated_lipids", "value"),
    Input("textarea_all_proteins", "value"),
    Input("textarea_regulated_proteins", "value"),
    Input("textarea_all_metabolites", "value"),
    Input("textarea_regulated_metabolites", "value"),
    Input("textarea_all_transcripts", "value"),
    Input("textarea_regulated_transcripts", "value"),
    Input("use_bounded_fatty_acyls", "checked"),
    Input("checkbox_use_lipids", "checked"),
    Input("checkbox_use_proteins", "checked"),
    Input("checkbox_use_metabolites", "checked"),
    Input("checkbox_use_transcripts", "checked"),
    Input("select_organism", "value"),
    Input("select_molecule_handling", "value"),
    Input("select_regulated_molecule_handling", "value"),
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
    all_transcripts_list,
    regulated_transcripts_list,
    use_bounded_fatty_acyls,
    with_lipids,
    with_proteins,
    with_metabolites,
    with_transcripts,
    select_organism,
    select_molecule_handling,
    select_regulated_molecule_handling,
    session_id,
):
    num_all_lipids = sum(len(line) > 0 for line in all_lipids_list.split("\n"))
    num_regulated_lipids = sum(len(line) > 0 for line in regulated_lipids_list.split("\n"))
    num_all_proteins = sum(len(line) > 0 for line in all_proteins_list.split("\n"))
    num_regulated_proteins = sum(len(line) > 0 for line in regulated_proteins_list.split("\n"))
    num_all_metabolites = sum(len(line) > 0 for line in all_metabolites_list.split("\n"))
    num_regulated_metabolites = sum(len(line) > 0 for line in regulated_metabolites_list.split("\n"))
    num_all_transcripts = sum(len(line) > 0 for line in all_transcripts_list.split("\n"))
    num_regulated_transcripts = sum(len(line) > 0 for line in regulated_transcripts_list.split("\n"))

    if session_id not in sessions:
        return (
            True,
            "Your session has expired. Please refresh the website.",
            f"Entries: {num_all_lipids}",
            f"Entries: {num_regulated_lipids}",
            f"Entries: {num_all_proteins}",
            f"Entries: {num_regulated_proteins}",
            f"Entries: {num_all_metabolites}",
            f"Entries: {num_regulated_metabolites}",
            f"Entries: {num_all_transcripts}",
            f"Entries: {num_regulated_transcripts}",
        )
    session = sessions[session_id]
    session.time = time.time()
    session.data_loaded = False
    session.use_bounded_fatty_acyls = use_bounded_fatty_acyls

    return (
        False,
        "",
        f"Entries: {num_all_lipids}",
        f"Entries: {num_regulated_lipids}",
        f"Entries: {num_all_proteins}",
        f"Entries: {num_regulated_proteins}",
        f"Entries: {num_all_metabolites}",
        f"Entries: {num_regulated_metabolites}",
        f"Entries: {num_all_transcripts}",
        f"Entries: {num_regulated_transcripts}",
    )



@callback(
    Output("loading_output", "children", allow_duplicate = True),
    Output("download_data", "data", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Input("icon_download_results", "n_clicks"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def download_table(
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

    return "", dcc.send_bytes(output.getvalue(), "GO_multiomics_results.xlsx"), False, ""



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
    Output("textarea_all_proteins", "value", allow_duplicate = True),
    Output("textarea_regulated_proteins", "value", allow_duplicate = True),
    Output("textarea_all_metabolites", "value", allow_duplicate = True),
    Output("textarea_regulated_metabolites", "value", allow_duplicate = True),
    Output("textarea_all_transcripts", "value", allow_duplicate = True),
    Output("textarea_regulated_transcripts", "value", allow_duplicate = True),
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
        "\n".join(examples[index]["bgp"]),
        "\n".join(examples[index]["regp"]),
        "\n".join(examples[index]["bgm"]),
        "\n".join(examples[index]["regm"]),
        "\n".join(examples[index]["bgt"]),
        "\n".join(examples[index]["regt"]),
        examples[index]["org"],
        len(examples[index]["bgl"]) > 0 or len(examples[index]["regl"]) > 0,
        len(examples[index]["bgp"]) > 0 or len(examples[index]["regp"]) > 0,
        len(examples[index]["bgm"]) > 0 or len(examples[index]["regm"]) > 0,
        len(examples[index]["bgt"]) > 0 or len(examples[index]["regt"]) > 0,
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
    background_proteins = session.background_proteins
    regulated_proteins = session.regulated_proteins
    background_metabolites = session.background_metabolites
    regulated_metabolites = session.regulated_metabolites
    background_transcripts = session.background_transcripts
    regulated_transcripts = session.regulated_transcripts

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

    if with_lipids:
        background_lipids = set(background_lipids.keys())
        lipid_table = [{"molecule_id": molecule, "molecule": molecule, "regulated": ("X" if molecule in regulated_lipids else "")} for molecule in molecules if molecule in background_lipids]

    if with_proteins:
        protein_table = [
            {
                "molecule_id": molecule,
                "molecule": molecule.replace("UNIPROT:", "") + f" ({ontology.proteins[molecule].name})" if molecule in ontology.proteins else "",
                "regulated": ("X" if molecule in regulated_proteins else "")
            } for molecule in molecules if molecule in background_proteins
        ]

    if with_metabolites:
        metabolite_table = [{"molecule_id": molecule, "molecule": molecule, "regulated": ("X" if molecule in regulated_metabolites else "")} for molecule in molecules if molecule in background_metabolites]

    if with_transcripts:
        transcript_table = [
            {
                "molecule_id": molecule,
                "molecule": molecule + f" ({ontology.transcripts[m].name})" if m in ontology.transcripts else "",
                "regulated": ("X" if molecule in regulated_transcripts else "")
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
            {'if': {'row_index': 1}, 'borderBottom': '1px solid black'},
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
        False,
        "",
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
    Output("sankey_controls", "style", allow_duplicate = True),
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
    State("sankey_controls", "style"),
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
    sankey_controls_style,
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
    sankey_controls_style["display"] = "none"
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

            fisher_data = session_data[list(term.term_id)[0]].fisher_data
            n, K, N = fisher_data[0] + fisher_data[2], fisher_data[0] + fisher_data[1], sum(fisher_data)
            overrepresented = fisher_data[0] >= ceil((n + 1) * (K + 1) / (N + 2))
            # blue hex: 3a4cc0
            # red hex: b30d03
            if overrepresented:
                h, s, l = 3, int(ratio * 100), 75 - int(ratio * 39) # 39 = 75 - 36
            else:
                h, s, l = 207, int(ratio * 100), 75 - int(ratio * 39) # 39 = 75 - 36
            color = f"hsl({h}, {s}%, {l}%)"

            term_session_data = session_data[list(term.term_id)[0]]
            regulated_metabolites = term_session_data.fisher_data[0]
            all_metabolites = regulated_metabolites + term_session_data.fisher_data[1]
            custom_text += f"<br>p-value: {sunburst_term.pvalue}<br>Regulated: {regulated_metabolites} / {all_metabolites}"

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
        False,
        "",
        barplot_controls_style,
        sunburst_controls_style,
        sankey_controls_style,
        no_update,
        no_update,
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
    Output("sankey_controls", "style", allow_duplicate = True),
    Output("barplot_numberinput_max_pvalue", "value", allow_duplicate = True),
    Output("barplot_numberinput_min_pvalue", "value", allow_duplicate = True),
    Input("sankey_results", "n_clicks"),
    Input("switch_sankey_regulated_only", "checked"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    State("sessionid", "children"),
    State("barplot_controls", "style"),
    State("sunburst_controls", "style"),
    State("sankey_controls", "style"),
    State("barplot_terms_wrapper", "style"),
    prevent_initial_call = True,
)
def open_sankeyplot(
    n_clicks,
    switch_sankey_regulated_only,
    row_data,
    selected_rows,
    session_id,
    barplot_controls_style,
    sunburst_controls_style,
    sankey_controls_style,
    barplot_terms_wrapper_style,
):
    if session_id == None or n_clicks == None or switch_sankey_regulated_only == None:
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
    selected_term_ids = [row["termid"] for row in selected_rows]
    ontology = session.ontology
    domains = session.domains
    terms = ontology.ontology_terms

    barplot_controls_style["display"] = "none"
    sunburst_controls_style["display"] = "none"
    sankey_controls_style["display"] = "block"

    background_lipids = session.background_lipids
    regulated_lipids = session.regulated_lipids
    background_proteins = session.background_proteins
    regulated_proteins = session.regulated_proteins
    background_metabolites = session.background_metabolites
    regulated_metabolites = session.regulated_metabolites
    background_transcripts = session.background_transcripts
    regulated_transcripts = session.regulated_transcripts

    #if switch_sankey_regulated_only:
    regulated_molecules = (
        (regulated_lipids if regulated_lipids else set()) |
        (regulated_proteins if regulated_proteins else set()) |
        (regulated_metabolites if regulated_metabolites else set()) |
        (regulated_transcripts if regulated_transcripts else set())
    )




    flow_data = {}
    flow_data_id = {}
    max_path_length = 0
    flow_layers = [TermType.INPUT_TERM]
    for target_term_id in selected_term_ids:
        target_term_id_single = target_term_id.split("|")[0]
        target_term = terms[target_term_id_single]

        for molecule in session.data[target_term_id_single].source_terms:
            if switch_sankey_regulated_only and molecule not in regulated_molecules: continue

            if regulated_lipids and molecule in background_lipids.keys():
                input_molecule = "Input lipid"
            elif regulated_proteins and molecule in background_proteins:
                input_molecule = "Input protein"
            elif regulated_metabolites and molecule in background_metabolites:
                input_molecule = "Input metabolite"
            elif regulated_transcripts and molecule in background_transcripts:
                input_molecule = "Input transcript"
            path_layers = [[TermType.INPUT_TERM, [{input_molecule}]]]

            # determine category path
            for i, term in enumerate(get_path(sessions[session_id].all_parent_nodes[molecule], target_term)):
                if type(term) == str: term_id = term
                else: term_id = list(term.term_id)[0]

                if not term_id in ontology.ontology_terms:
                    if term_id in background_lipids and type(background_lipids[term_id]) == LipidAdduct:
                        lipid_category = background_lipids[term_id].get_lipid_string(LipidLevel.CATEGORY)
                        if not path_layers or path_layers[-1][0] != TermType.LIPID_SPECIES:
                            path_layers.append((TermType.LIPID_SPECIES, [{lipid_category}]))
                        else:
                            path_layers[-1][1].append({lipid_category})

                if term_id in ontology.ontology_terms:
                    term_type = ontology.ontology_terms[term_id].term_type
                    if term_type == TermType.LIPID_CLASS: term_type = TermType.LIPID_SPECIES
                    elif term_type in {TermType.UNREVIEWED_PROTEIN, TermType.ENSEMBLE_PROTEIN}: term_type = TermType.REVIEWED_PROTEIN
                    elif term_type == TermType.ENSEMBLE_GENE: term_type = TermType.GENE
                    elif term_type == TermType.METABOLITE: term_type = TermType.LIPID_SPECIES

                    if term_type != TermType.UNCLASSIFIED_TERM:
                        if not path_layers or path_layers[-1][0] != term_type:
                            path_layers.append((term_type, [ontology.ontology_terms[term_id].categories]))
                        else:
                            path_layers[-1][1].append(ontology.ontology_terms[term_id].categories)

                    else:
                        if not path_layers or path_layers[-1][0] != ontology.ontology_terms[term_id].domain:
                            path_layers.append((ontology.ontology_terms[term_id].domain, [{ontology.ontology_terms[term_id].name}]))
                        else:
                            path_layers[-1][1].append({ontology.ontology_terms[term_id].name})

            if len(path_layers) < 3: continue
            max_path_length = max(max_path_length, len(path_layers))
            flow_i = 0
            for pos in range(len(path_layers) - 1):
                path_key_prev, layers_list_prev = path_layers[pos]
                categories_prev = layers_list_prev[0] if type(path_key_prev) == TermType else layers_list_prev[-1]
                path_key_next, layers_list_next = path_layers[pos + 1]
                categories_next = layers_list_next[0] if type(path_key_next) == TermType else layers_list_next[-1]

                while flow_i + 1 < len(flow_layers) and flow_layers[flow_i] != path_key_prev:
                    flow_layers.insert(flow_i, path_key_prev)
                    flow_i += 1

                for category_next in categories_next:
                    if category_next not in flow_data_id:
                        flow_data_id[category_next] = [len(flow_data_id), path_key_next]
                        flow_data[flow_data_id[category_next][0]] = {}

                for category_prev in categories_prev:
                    if category_prev not in flow_data_id:
                        flow_data_id[category_prev] = [len(flow_data_id), path_key_prev]
                        flow_data[flow_data_id[category_prev][0]] = {}
                    last_flow_data = flow_data[flow_data_id[category_prev][0]]
                    for category_next in categories_next:
                        term_node_id = flow_data_id[category_next][0]
                        if term_node_id not in last_flow_data:
                            last_flow_data[term_node_id] = 0
                        last_flow_data[term_node_id] += 1

    pastel_colors = ["#FFD1DC", "#AAF0D1", "#FFB347", "#B5EAAA", "#CBAACB", "#FDFD96", "#CFCFC4", "#E3E4FA", "#AEC6CF", "#FF6961", "#FFE5B4", "#77DD77"]
    fig_source, fig_target, fig_value = [], [], []
    node_colors, edge_colors = [], []
    fig_label = sorted([f for f in flow_data_id], key = lambda f: flow_data_id[f][0])
    for i, source_node in enumerate(fig_label):
        source_node_id, path_key = flow_data_id[source_node]
        print(source_node_id, path_key)
        node_colors.append(pastel_colors[i % len(pastel_colors)])
        edge_colors += [pastel_colors[i % len(pastel_colors)]] * len(flow_data[source_node_id])
        for target_node_id, value in flow_data[source_node_id].items():
            fig_source.append(source_node_id)
            fig_target.append(target_node_id)
            fig_value.append(value)

    fig = go.Figure(
        go.Sankey(
            arrangement = 'snap',
            node = dict(
                label = fig_label,
                align = "right",
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
        "90%",
        barplot_terms_wrapper_style,
        False,
        "",
        barplot_controls_style,
        sunburst_controls_style,
        sankey_controls_style,
        no_update,
        no_update,
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
    Output("sankey_controls", "style", allow_duplicate = True),
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
    State("sankey_controls", "style"),
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
    sankey_controls_style,
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
    sankey_controls_style["display"] = "none"

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

    pvalues = -np.log10([session_data[term_id].pvalue_corrected for term_id in selected_term_ids])
    number_entities = np.array([len(session_data[term_id].source_terms) for term_id in selected_term_ids])
    number_regulated_entities = np.array([session_data[term_id].fisher_data[0] for term_id in selected_term_ids])

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
            add_arc(
                fig,
                pvalue_start_angle,
                pvalue_end_angle,
                pvalue_arc_inner_radius,
                pvalue_arc_outer_radius,
                domain_colors[list(result.term.domain)[0]][1],
                hoverinfo = "p-value: {:.6g}".format(result.pvalue_corrected)
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
        False,
        "",
        barplot_controls_style,
        sunburst_controls_style,
        sankey_controls_style,
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
    Output("sankey_controls", "style", allow_duplicate = True),
    Input("histogram_results", "n_clicks"),
    State("sessionid", "children"),
    State("barplot_controls", "style"),
    State("sunburst_controls", "style"),
    State("sankey_controls", "style"),
    State("barplot_terms_wrapper", "style"),
    prevent_initial_call = True,
)
def open_histogram(
    n_clicks,
    session_id,
    barplot_controls_style,
    sunburst_controls_style,
    sankey_controls_style,
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
    sankey_controls_style["display"] = "none"

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
        False,
        "",
        barplot_controls_style,
        sunburst_controls_style,
        sankey_controls_style,
    )



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
        href, term_name = ".", ""

        if term_id in ontology.ontology_terms:
            term_name = term.name if type(term) == OntologyTerm else ontology.ontology_terms[term_id].name
            if term_id.startswith("GO:"):
                href = "https://www.ebi.ac.uk/QuickGO/term/" + term_id

            elif term_id.startswith("SMP"):
                href = "https://pathbank.org/view/" + term_id

            elif term_id.startswith("LION:") or term_id.startswith("CAT:"):
                href = "https://bioportal.bioontology.org/ontologies/LION?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F" + term_id.replace(":", "_")

            elif term_id.startswith("RHEA:"):
                href = f"https://www.rhea-db.org/rhea/{term_id[5:]}"

            elif term_id.startswith("UNIPROT:"):
                href = f"https://www.uniprot.org/uniprotkb/{term_id[8:]}/entry"

            elif term_id.startswith("CHEBI:"):
                href = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={term_id}"

            elif term_id.startswith("DOID:"):
                href = f"https://disease-ontology.org/?id={term_id}"

            elif term_id.startswith("MONDO:"):
                href = f"https://monarchinitiative.org/{term_id}"

            elif term_id.startswith("HGNC:"):
                href = f"https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/{term_id}"

            elif term_id.startswith("HP:"):
                href = f"https://hpo.jax.org/browse/term/{term_id}"

            elif term_id.startswith("NCBI:"):
                href = f"https://www.ncbi.nlm.nih.gov/gene/{term_id.split(':')[1]}"

            elif term_id.startswith("ENS") or term_id.startswith("WBGene") or term_id.startswith("FBgn"):
                href = f"https://www.ensembl.org/id/{term_id}"

            elif term_id.startswith("R-"):
                href = f"https://reactome.org/PathwayBrowser/#/{term_id}"

            elif term_id.startswith("LM"):
                href = f"https://www.lipidmaps.org/databases/lmsd/{term_id}"

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
    return False, "", term_path



ns = api.namespace("api", description = "Operations")

# Model for GO analysis
enrichment_model = api.model("Enrichment", {
    "background_lipids": fields.List(
        fields.String,
        description = "All lipid names in experiment (background, required)",
        example = ["ST 27:1;O", "TG 16:0/18:1/20;0", "10-HDoHE", "12-HEPE", "12-HETE", "12-HHTrE", "12-OxoETE", "13-HODE", "15-HETE", "5,6-diHETE", "8,9-EET", "8-HETE", "9-HODE"],

    ),
    "regulated_lipids": fields.List(
        fields.String,
        description = "All regulated lipid names in experiment",
        example = ["ST 27:1;O", "TG 16:0/18:1/20;0", "10-HDoHE", "12-HEPE", "12-HETE", "12-HHTrE"],
    ),
    "background_proteins": fields.List(
        fields.String,
        description = "All protein accessions in experiment (background)",
        example = ["Q8TF30", "Q15465", "P37231", "P12958"],
    ),
    "regulated_proteins": fields.List(
        fields.String,
        description = "All regulated protein accession in experiment",
        example = ["P37231", "P12958"],
    ),
    "background_metabolites": fields.List(
        fields.String,
        description = "All metabolite ChEBI Ids in experiment (background)",
        example = ["CHEBI:46053", "CHEBI:27732", "methanol","CHEBI:27958", "CHEBI:17234"],
    ),
    "regulated_metabolites": fields.List(
        fields.String,
        description = "All regulated metabolite ChEBI Ids in experiment",
        example = ["CHEBI:27958", "CHEBI:17234"],
    ),
    "background_transcripts": fields.List(
        fields.String,
        description = "All ensembl Ids in experiment (background)",
        example = ["ENST00000515318.6", "ENST00000435425.1", "ENST00000438682.6", "ENSG00000156006.6"],
    ),
    "regulated_transcripts": fields.List(
        fields.String,
        description = "All regulated ensembl Ids in experiment",
        example = ["ENST00000438682.6", "ENSG00000156006.6"],
    ),
    "organism_taxonomy": fields.String(
        description = "Select organism (Taxonomic number), default: '9606' (Homo sapiens)",
        enum = [v for k, v in organisms.items()],
        example = "9606",
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
        description = f"Handling of unrecognizable molecules, default: {MOLECULE_HANDLING_IGNORE}",
        enum = [m["value"] for m in molecule_handling],
        example = MOLECULE_HANDLING_IGNORE,
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

            organism_api = data.get("organism_taxonomy", "9606")
            domains_api = data.get("domains", ["biological_process"])
            accepted_domains = {d.lower().replace(" ", "_"): d for d in enrichment_ontologies[INIT_ORGANISM].domains}
            pvalue_correction_api = data.get("pvalue_correction", "fdr_bh")
            term_representation_api = data.get("term_representation", "greater")
            unrecognizable_molecules_api = data.get("unrecognizable_molecules", MOLECULE_HANDLING_IGNORE)
            non_background_molecules_api = data.get("non_background_molecules", MOLECULE_HANDLING_REMOVE)
            bounded_fatty_acyls_api = data.get("bounded_fatty_acyls", False)

            if type(organism_api) not in {str, int} or (type(organism_api) == str and not organism_api.isnumeric()):
                return {"error_message": "'organism_taxonomy' must be a (string) number", "result": []}, 422

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
            check_failed, error_message, molecule_tables = check_user_input(
                omics_included,
                omics_lists,
                ontology,
                unrecognizable_molecules_api,
                non_background_molecules_api,
            )

            if check_failed:
                return {"error_message": error_message, "result": []}, 422

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
            ) = molecule_tables

            session = SessionEntry()
            session.time = time.time()
            session.use_bounded_fatty_acyls = bounded_fatty_acyls_api
            ontology.set_background(session, lipid_dict = lipidome, protein_set = proteome, metabolite_set = metabolome, transcript_set = transcriptome)
            session.ontology = ontology
            session.data_loaded = True

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
            results = ontology.enrichment_analysis(session, target_set, domains_api, term_representation_api)
            session.result = results

            if pvalue_correction_api != "no" and len(results) > 1:
                pvalues = [r.pvalue for r in results]
                pvalues = multipletests(pvalues, method = pvalue_correction_api)[1]
                for pvalue, r in zip(pvalues, results): r.pvalue_corrected = pvalue

            results.sort(key = lambda row: (row.pvalue_corrected, row.term.name))

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
            return {"error_message": f"{e}", "result": []}, 500






