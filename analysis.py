from dash import Dash, dcc, html, Input, Output, State, callback, exceptions, no_update, MATCH, ALL, callback_context, clientside_callback, dash_table
import dash_mantine_components as dmc
import plotly.graph_objs as go
import dash_ag_grid as dag
from dash_iconify import DashIconify
import json
import numpy as np
import pandas as pd
import io
import os
from statsmodels.stats.multitest import multipletests
from EnrichmentDataStructure import EnrichmentOntology, current_path, SessionEntry
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import LipidParser
import logging
import pathlib
import time
import hashlib
import threading
import freetype


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
SESSION_DURATION_TIME = 60 * 60 # one hour
domain_colors = {
    "Biological process": ["#fc7255", "#fc947e"],
    "Cellular component": ["#c86932", "#c87c50"],
    "Metabolic and signalling pathway": ["#4daf4a", "#68af66"],
    "Molecular function": ["#ffeda0", "#fff4c7"],
    "Physical or chemical properties": ["#377eb8", "#538ab8"],
    "Disease": ["#377eb8", "#538ab8"],
}
upregulated_color = "#F49E4C"
downregulated_color = "#87BCDE"
entities_number_color = "#F49E4C"
INNER_CIRCLE = 300
CIRCLE_WIDTH = 300



def annotate_polar(
    figure,
    annotations,
    text_to_add,
    mid_angle,
    angle_width,
    distance,
    font_size = 14,
    max_radius = 400,
    color = "black"
):
    text_width = sum(char_sizes[font_size][ord(char)] for char in text_to_add)
    if 360 / (2 * max_radius * np.pi) * text_width > angle_width:
        text_width = angle_width / 360 * (2 * max_radius * np.pi)

    # Calculate positions of each character along the arc
    start_pos = (2 * max_radius * np.pi) / 360 * mid_angle + text_width / 2
    # Compute the positions (x, y) for each character
    x, y, rotation_angles, text = [], [], [], []
    for i, char in enumerate(text_to_add):

        char_size = char_sizes[font_size][ord(char)]
        current_pos = start_pos - char_size / 2
        start_pos -= char_size

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



def add_arc(figure, start_angle, end_angle, arc_inner_radius, arc_outer_radius, color, arc_resolution = 50):
    theta_outer = np.linspace(start_angle, end_angle, arc_resolution)
    theta_inner = theta_outer[::-1]

    r_outer = [arc_outer_radius] * arc_resolution
    r_inner = [arc_inner_radius] * arc_resolution

    r_arc = np.concatenate([r_outer, r_inner])
    theta_arc = np.concatenate([theta_outer, theta_inner])

    figure.add_trace(go.Scatterpolar(
        r = r_arc,
        theta = theta_arc,
        mode = 'lines',
        fill = 'toself',
        fillcolor = color,
        line = dict(color = color),
        opacity = 1.0,
        hoverinfo = 'skip',
        showlegend = False,
    ))



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



organisms = {
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

try:
    from organisms import organisms
except Exception as e:
    pass

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
app.title = "GO multiomics"

lipid_parser = LipidParser()
enrichment_ontologies = {}
for tax_name, tax_id in organisms.items():
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
                dcc.Graph(
                    id = "barplot_terms",
                    config = plotly_config,
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
                        ],
                        cols = 2,
                    ),
                    id = "barplot_controls",
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
            children = dmc.ScrollArea([
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
                                        dmc.ActionIcon(
                                            DashIconify(icon = "simple-icons:helix", width = 16),
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
                                        data=[
                                            {"value": "no", "label": "No correction"},
                                            {"value": "bonferroni", "label": "Bonferroni"},
                                            {"value": "fdr_bh", "label": "Benjamini Hochberg"},
                                        ],
                                        value = "fdr_bh",
                                    ),
                                    dmc.Select(
                                        id = "select_term_representation",
                                        data = [
                                            {"value": "two-sided", "label": "Term over/under-represented"},
                                            {"value": "less", "label": "Term under-represented"},
                                            {"value": "greater", "label": "Term over-represented"},
                                        ],
                                        value = "two-sided",
                                        label = "Term representation:",
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
                                    html.Div(
                                        dmc.Switch(
                                            id = "use_bounded_fatty_acyls",
                                            checked = True,
                                            label = "Use bounded fatty acyls for analysis, too",
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
                                        DashIconify(icon="ion:bar-chart-sharp", width = 20),
                                        id = "histogram_results",
                                        title = "Show p-value histogram",
                                        disabled = True,
                                        style = {"marginRight": "5px"},
                                    ),
                                    dmc.ActionIcon(
                                        DashIconify(icon="tdesign:chart-column-filled", width = 20),
                                        id = "chart_results",
                                        title = "Show p-value chart",
                                        disabled = True,
                                        style = {"marginRight": "5px"},
                                    ),
                                    " ",
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
                                'field': "count",
                                "headerName": "Count",
                                "width": 80,
                                "headerTooltip": "Number of regulated associated molecules / Number of all associated molecules",
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
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("background_lipids", "children", allow_duplicate = True),
    Output("regulated_lipids", "children", allow_duplicate = True),
    Output("background_proteins", "children", allow_duplicate = True),
    Output("regulated_proteins", "children", allow_duplicate = True),
    Output("background_metabolites", "children", allow_duplicate = True),
    Output("regulated_metabolites", "children", allow_duplicate = True),
    Output("histogram_results", "disabled", allow_duplicate = True),
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
    State("select_term_representation", "value"),
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
    do_activate_alert = True
    histogram_disabled = True

    if session_id not in sessions:
        return "", [], do_activate_alert, "Your session has expired. Please refresh the website.", "", "", "", "", "", "", histogram_disabled

    logger.info(f"Enrichment session: {session_id}")
    sessions[session_id].time = time.time()

    if not with_lipids and not with_proteins and not with_metabolites:
        return "", [], do_activate_alert, "No omics data is selected.", "", "", "", "", "", "", histogram_disabled

    ontology = enrichment_ontologies[organism]

    target_set = set()
    lipidome, regulated_lipids = {}, set()
    proteome, regulated_proteins = set(), set()
    metabolome, regulated_metabolites = set(), set()
    if with_lipids:
        if len(all_lipids_list) == 0:
            return "", [], do_activate_alert, "Please paste lipid names into the first text area.", "", "", "", "", "", "", histogram_disabled

        if len(regulated_lipids_list) == 0:
            return "", [], do_activate_alert, "Please paste lipid names into the second text area.", "", "", "", "", "", "", histogram_disabled

        for lipid_name in all_lipids_list.split("\n"):
            if len(lipid_name) == 0: continue
            try:
                lipid = lipid_parser.parse(lipid_name)
                lipidome[lipid_name] = lipid

            except Exception as e:
                if not ignore_unrecognizable_molecules:
                    return "", [], do_activate_alert, f"Lipid name '{lipid_name}' unrecognizable in background list! Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", "", histogram_disabled

        for lipid_name in regulated_lipids_list.split("\n"):
            if len(lipid_name) == 0: continue
            if lipid_name not in lipidome:
                if ignore_unknown: continue
                return "", [], do_activate_alert, f"The regulated lipid '{lipid_name}' does not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", "", "", "", "", "", "", histogram_disabled
            regulated_lipids.add(lipid_name)

        if len(lipidome) == 0:
            return "", [], do_activate_alert, "No background lipid left after lipid recognition.", "", "", "", "", "", "", histogram_disabled

        if len(regulated_lipids) == 0:
            return "", [], do_activate_alert, "No regulated lipid left after lipid recognition.", "", "", "", "", "", "", histogram_disabled

        if len(regulated_lipids) > len(lipidome):
            return "", [], do_activate_alert, "Length of regulated lipid list must be smaller than background list.", "", "", "", "", "", "", histogram_disabled

        left_lipids = regulated_lipids - lipidome.keys()
        if len(left_lipids) > 0:
            if ignore_unknown:
                for lipid_name in left_lipids:
                    del lipidome[lipid_name]
            else:
                return "", [], do_activate_alert, "The regulated lipid" + (' ' if len(left_lipids) == 1 else 's ') + "'" + "', '".join(left_lipids) + ("' does" if len(left_lipids) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", "", "", "", "", "", "", histogram_disabled

        target_set |= regulated_lipids

    if with_proteins:
        if len(all_proteins_list) == 0:
            return "", [], do_activate_alert, "Please paste protein accession into the first text area.", "", "", "", "", "", "", histogram_disabled

        if len(regulated_proteins_list) == 0:
            return "", [], do_activate_alert, "Please paste protein accessions into the second text area.", "", "", "", "", "", "", histogram_disabled

        proteome = set(protein for protein in all_proteins_list.split("\n") if len(protein) > 0)
        regulated_proteins = set(protein for protein in regulated_proteins_list.split("\n") if len(protein) > 0)

        left_proteins = proteome - ontology.clean_protein_ids
        if len(left_proteins) > 0:
            if ignore_unrecognizable_molecules:
                proteome -= left_proteins
            else:
                return "", [], do_activate_alert, "The protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' is" if len(left_proteins) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", "", histogram_disabled

        left_proteins = regulated_proteins - ontology.clean_protein_ids
        if len(left_proteins) > 0:
            if ignore_unrecognizable_molecules:
                regulated_proteins -= left_proteins
            else:
                return "", [], do_activate_alert, "The protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' is" if len(left_proteins) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", "", histogram_disabled

        left_proteins = regulated_proteins - proteome
        if len(left_proteins) > 0:
            if ignore_unknown:
                proteome -= left_proteins
            else:
                return "", [], do_activate_alert, "The regulated protein" + (' ' if len(left_proteins) == 1 else 's ') + "'" + "', '".join(left_proteins) + ("' does" if len(left_proteins) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", "", "", "", "", "", "", histogram_disabled

        if len(proteome) == 0:
            return "", [], do_activate_alert, "No background protein left after protein recognition.", "", "", "", "", "", "", histogram_disabled

        if len(regulated_proteins) == 0:
            return "", [], do_activate_alert, "No regulated protein left after protein recognition.", "", "", "", "", "", "", histogram_disabled

        if len(regulated_proteins) > len(proteome):
            return "", [], do_activate_alert, "Length of regulated protein list must be smaller than background list.", "", "", "", "", "", "", histogram_disabled

        proteome = set([f"UNIPROT:{protein}" for protein in proteome])
        regulated_proteins = set([f"UNIPROT:{protein}" for protein in regulated_proteins])
        target_set |= regulated_proteins

    if with_metabolites:
        if len(all_metabolites_list) == 0:
            return "", [], do_activate_alert, "Please paste metabolite ChEBI Ids into the first text area.", "", "", "", "", "", "", histogram_disabled

        if len(regulated_metabolites_list) == 0:
            return "", [], do_activate_alert, "Please paste metabolite ChEBI Ids into the second text area.", "", "", "", "", "", "", histogram_disabled

        metabolome = set(metabolite for metabolite in all_metabolites_list.split("\n") if len(metabolite) > 0)
        regulated_metabolites = set(metabolite for metabolite in regulated_metabolites_list.split("\n") if len(metabolite) > 0)

        left_metabolites = metabolome - ontology.clean_metabolite_ids - ontology.metabolites.keys()
        left_metabolites -= set([m for m in left_metabolites if m.lower() in ontology.metabolite_names.keys()])
        if len(left_metabolites) > 0:
            if ignore_unrecognizable_molecules:
                metabolome -= left_metabolites
            else:
                return "", [], do_activate_alert, "The metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' is" if len(left_metabolites) == 1 else "' are") + " unrecognizable in the background list. Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", "", histogram_disabled

        left_metabolites = regulated_metabolites - ontology.clean_metabolite_ids - ontology.metabolites.keys()
        left_metabolites -= set([m for m in left_metabolites if m.lower() in ontology.metabolite_names.keys()])
        if len(left_metabolites) > 0:
            if ignore_unrecognizable_molecules:
                regulated_metabolites -= left_metabolites
            else:
                return "", [], do_activate_alert, "The metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' is" if len(left_metabolites) == 1 else "' are") + " unrecognizable in the regulated. Maybe enable the 'Ignore unrecognizable molecules' option.", "", "", "", "", "", "", histogram_disabled

        left_metabolites = regulated_metabolites - metabolome
        if len(left_metabolites) > 0:
            if ignore_unknown:
                metabolome -= left_metabolites
            else:
                return "", [], do_activate_alert, "The regulated metabolite" + (' ' if len(left_metabolites) == 1 else 's ') + "'" + "', '".join(left_metabolites) + ("' does" if len(left_metabolites) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated molecules that aren't in background' option.", "", "", "", "", "", "", histogram_disabled

        if len(metabolome) == 0:
            return "", [], do_activate_alert, "No background metabolite left after metabolite recognition.", "", "", "", "", "", "", histogram_disabled

        if len(regulated_metabolites) == 0:
            return "", [], do_activate_alert, "No regulated metabolite left after metabolite recognition.", "", "", "", "", "", "", histogram_disabled

        if len(regulated_metabolites) > len(metabolome):
            return "", [], do_activate_alert, "Length of regulated metabolite list must be smaller than background list.", "", "", "", "", "", "", histogram_disabled

        metabolome = set([f"CHEBI:{metabolite}" if (type(metabolite) == int or metabolite[:6] != "CHEBI:") and metabolite.lower() not in ontology.metabolite_names else metabolite for metabolite in metabolome])
        regulated_metabolites = set([f"CHEBI:{metabolite}" if (type(metabolite) == int or metabolite[:6] != "CHEBI:") and metabolite.lower() not in ontology.metabolite_names else metabolite for metabolite in regulated_metabolites])
        target_set |= regulated_metabolites

    if len(domains) == 0:
        return "", [], do_activate_alert, "No domain(s) selected.", "", "", "", "", "", "", histogram_disabled

    if sessions[session_id].data_loaded == False:
        ontology.set_background(sessions[session_id], lipid_dict = lipidome, protein_set = proteome, metabolite_set = metabolome)
        sessions[session_id].data_loaded = True

    ontology.enrichment_analysis(sessions[session_id], target_set, domains, term_regulation)
    results = sessions[session_id].result

    if correction_method != "no" and len(results) > 0:
        pvalues = [r.pvalue for r in results]
        pvalues = multipletests(pvalues, method = correction_method)[1]
        for pvalue, r in zip(pvalues, results): r.pvalue_corrected = pvalue

    results.sort(key = lambda row: row.pvalue_corrected)


    data = []
    for result in results:
        data.append(
            {
                "domain": result.term.domain,
                "term": result.term.name,
                "termid": result.term.term_id,
                "count": f"{result.fisher_data[0]} / {len(set(result.source_terms))}",
                "pvalue": "{:.6g}".format(result.pvalue_corrected)
            }
        )

    histogram_disabled = False
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
        histogram_disabled,
    )



@callback(
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("num_all_lipids", "children", allow_duplicate = True),
    Output("num_regulated_lipids", "children", allow_duplicate = True),
    Output("num_all_proteins", "children", allow_duplicate = True),
    Output("num_regulated_proteins", "children", allow_duplicate = True),
    Output("num_all_metabolites", "children", allow_duplicate = True),
    Output("num_regulated_metabolites", "children", allow_duplicate = True),
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
    num_all_lipids = sum(len(line) > 0 for line in all_lipids_list.split("\n"))
    num_regulated_lipids = sum(len(line) > 0 for line in regulated_lipids_list.split("\n"))
    num_all_proteins = sum(len(line) > 0 for line in all_proteins_list.split("\n"))
    num_regulated_proteins = sum(len(line) > 0 for line in regulated_proteins_list.split("\n"))
    num_all_metabolites = sum(len(line) > 0 for line in all_metabolites_list.split("\n"))
    num_regulated_metabolites = sum(len(line) > 0 for line in regulated_metabolites_list.split("\n"))

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
        )

    sessions[session_id].time = time.time()
    sessions[session_id].data_loaded = False

    return (
        False,
        "",
        f"Entries: {num_all_lipids}",
        f"Entries: {num_regulated_lipids}",
        f"Entries: {num_all_proteins}",
        f"Entries: {num_regulated_proteins}",
        f"Entries: {num_all_metabolites}",
        f"Entries: {num_regulated_metabolites}",
    )



@callback(
    Output("loading_output", "children", allow_duplicate = True),
    Output("download_data", "data", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
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
        return (
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
        )

    sessions[session_id].time = time.time()
    domains = []
    term_ids = []
    terms = []
    counts = []
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
        counts.append(row["count"])
        pvalues.append(row["pvalue"])

    df = pd.DataFrame({"Domain": domains, "Term ID": term_ids, "Term": terms, "Count": counts, "pValue": pvalues})
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
        molecules = data[term_id].source_terms.keys()

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

    return "", dcc.send_bytes(output.getvalue(), "GO_multiomics_results.xlsx"), False, ""



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
    Output("term_molecules_modal_id", "children", allow_duplicate = True),
    Output("lipid_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("protein_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("metabolite_tab_modal_tab", "disabled", allow_duplicate = True),
    Output("term_molecules_modal_tab", "value", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("term_path_area", "children", allow_duplicate = True),
    Output("contingency_table", "children", allow_duplicate = True),
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
        return (
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

    sessions[session_id].time = time.time()
    if "rowId" not in row_data:
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
    molecules = sorted(list(result.source_terms))

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

    tab_value = "lipid_tab_modal"
    if not with_lipids:
        tab_value = "protein_tab_modal"
        if not with_proteins:
            tab_value = "metabolite_tab_modal"

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
        result.term.term_id,
        not with_lipids,
        not with_proteins,
        not with_metabolites,
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
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("barplot_controls", "style", allow_duplicate = True),
    Input("chart_results", "n_clicks"),
    Input("barplot_numberinput_connect_ths", "value"),
    Input("barplot_numberinput_font_size", "value"),
    Input("barplot_select_name", "value"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    State("sessionid", "children"),
    State("barplot_controls", "style"),
    prevent_initial_call = True,
)
def open_barplot(n_clicks, jaccard_ths, font_size, bar_label, row_data, selected_rows, session_id, controls_style):
    if session_id == None or jaccard_ths == None or n_clicks == None or font_size == None or bar_label not in {"id", "name"}:
        raise exceptions.PreventUpdate

    if session_id not in sessions:
        return (
            no_update,
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            no_update,
        )

    fig = go.Figure()
    id_position = {row["termid"]: i for i, row in enumerate(row_data)}
    session_data = sessions[session_id].data
    selected_term_ids = sorted([row["termid"] for row in selected_rows], key = lambda x: id_position[x])
    pvalues = -np.log10([session_data[term_id].pvalue_corrected for term_id in selected_term_ids])
    number_entities = np.array([len(session_data[term_id].source_terms) for term_id in selected_term_ids])
    number_regulated_entities = np.array([session_data[term_id].fisher_data[0] for term_id in selected_term_ids])
    domains = [session_data[term_id].term.domain for term_id in selected_term_ids]
    term_names = [session_data[term_id].term.name for term_id in selected_term_ids]
    print(term_names)
    term_domain_colors_def = [domain_colors[domain][0] for domain in domains]
    term_domain_colors_bar = [domain_colors[domain][1] for domain in domains]
    max_axis = int(max(pvalues) + 1)
    controls_style["display"] = "block"
    jaccard_ths /= 100

    custom_data = [[s, row["term"], row["pvalue"], r, n, d] for s, row, r, n, d in zip(selected_term_ids, selected_rows, number_regulated_entities, number_entities, domains)]
    inner_margin = INNER_CIRCLE / CIRCLE_WIDTH * max(pvalues)
    pvalues = pvalues + inner_margin
    n = len(pvalues)
    max_height = max(pvalues) * 1.02
    angle_gap = 360 / n
    angles = [i * angle_gap for i in range(n)]
    bar_width = angle_gap - 1
    arc_outer_radius = max_height * 1.3
    axis_step_size = int(max(1, 0.7 * np.sqrt(max_axis)))

    # Arc settings
    description_arc_inner_radius = arc_outer_radius * 0.9
    entities_arc_outer_radius = arc_outer_radius * 0.86
    entities_arc_inner_radius = arc_outer_radius * 0.80
    #regulated_arc_inner_radius = arc_outer_radius * 0.76
    #regulated_arc_outer_radius = arc_outer_radius * 0.70

    #number_entities_sizes = np.log(number_entities + 2)
    #number_regulated_entities_sizes = np.log(number_regulated_entities + 2)
    number_entities_sizes = number_entities + 5
    number_regulated_entities_sizes = number_regulated_entities + 5

    annotations = []
    mnes = np.max(number_entities_sizes)
    number_entities_sizes = number_entities_sizes / mnes
    number_regulated_entities_sizes = number_regulated_entities_sizes / mnes

    fig = go.Figure(data = go.Bar())

    # Add radial bars
    fig.add_trace(go.Barpolar(
        r = [max_height] * n,
        theta = angles,
        width = [bar_width] * n,
        marker_color = "#eeeeee",
        hoverinfo = 'skip',
        showlegend = False,
    ))

    fig.add_trace(go.Barpolar(
        r = pvalues,
        theta = angles,
        width = [bar_width] * n,
        marker_color = term_domain_colors_bar,
        customdata = custom_data,
        hovertemplate = "Term ID: %{customdata[0]}<br />Term: %{customdata[1]}<br />p-value: %{customdata[2]}<br />Regulated molecules: %{customdata[3]}<br />Associated molecules: %{customdata[4]}<br />Domain: %{customdata[5]}<extra></extra>",
    ))

    # Add thick arcs above bars
    annotation_label = selected_term_ids if bar_label == "id" else term_names
    for i in range(n):
        center_angle = angles[i]
        description_start_angle = center_angle + bar_width / 2
        description_end_angle = center_angle - bar_width / 2
        add_arc(
            fig,
            description_start_angle,
            description_end_angle,
            description_arc_inner_radius,
            arc_outer_radius,
            term_domain_colors_def[i],
        )
        arc_mid_radius = (description_arc_inner_radius + arc_outer_radius) / 2
        annotate_polar(fig, annotations, annotation_label[i], angles[i], bar_width, arc_mid_radius, font_size = font_size)

        entities_start_angle = center_angle + bar_width / 2
        entities_mid_angle = entities_start_angle
        if number_regulated_entities[i] > 0:
            entities_end_angle = entities_start_angle - number_regulated_entities_sizes[i] * bar_width
            add_arc(
                fig,
                entities_start_angle,
                entities_end_angle,
                entities_arc_inner_radius,
                entities_arc_outer_radius,
                upregulated_color,
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
                downregulated_color,
            )

        arc_mid_radius = (entities_arc_inner_radius + entities_arc_outer_radius) / 2
        arc_angle = (entities_start_angle + entities_end_angle) / 2
        entities_number_text = f"{str(int(number_regulated_entities[i]))} / {str(int(number_entities[i]))}"
        annotate_polar(fig, annotations, entities_number_text, arc_angle, abs(entities_start_angle - entities_end_angle), arc_mid_radius, font_size = font_size)

    # Add white donut hole
    theta_circle = np.linspace(0, 360, 100)
    r_circle = [inner_margin] * 100

    fig.add_trace(go.Scatterpolar(
        r = r_circle,
        theta = theta_circle,
        mode = 'lines',
        fill = 'toself',
        fillcolor = 'white',
        line = dict(color = 'white'),
        hoverinfo = 'skip',
        showlegend = False
    ))

    source_terms = [
        set(session_data[term_id].source_terms.keys()) for term_id in selected_term_ids
    ]

    for i in range(0, n - 1):
        for j in range(i + 1, n):
            i_terms = source_terms[i]
            j_terms = source_terms[j]

            jaccard = len(i_terms & j_terms) / len(i_terms | j_terms)
            if jaccard < jaccard_ths: continue
            if jaccard_ths < 1:
                jaccard = int(min(90, 100 - 100 * (jaccard - jaccard_ths) / (1 - jaccard_ths)))
            else:
                jaccard = 90


            x0 = inner_margin * np.cos(np.radians(angles[i]))
            y0 = inner_margin * np.sin(np.radians(angles[i]))
            x2 = inner_margin * np.cos(np.radians(angles[j]))
            y2 = inner_margin * np.sin(np.radians(angles[j]))
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

            xx, yy = bezier_points[:, 0], bezier_points[:, 1]
            rr = np.sqrt(xx**2 + yy**2)
            theta = np.mod(np.arctan2(yy, xx) * (180 / np.pi) + 360, 360)

            fig.add_trace(go.Scatterpolar(
                r = rr,
                theta = theta,
                mode = 'lines',
                line_shape = 'spline',
                line = dict(color = f'hsv(0,0,{jaccard})', width = 2),
                name = 'Curved Line',
                hoverinfo = 'skip',
                showlegend = False
            ))

    # Layout settings with visible radial axis (height of bars)
    fig.update_layout(
        annotations = annotations,
        polar = dict(
            domain = dict(x = [0, 1], y = [0, 1]),  # same domain as Cartesian
            barmode = 'overlay',
            radialaxis = dict(
                visible = True,
                tickvals = [i + inner_margin for i in range(0, max_axis, axis_step_size)],
                ticktext = list(range(0, max_axis, axis_step_size)),
                showticklabels = True,
                ticks = 'inside',
                showline = False,
                showgrid = True,
                range = [0, arc_outer_radius],
            ),
            angularaxis = dict(
                visible = False,       # Hide angular axis ticks, labels, and grid
                showticklabels = False,
                ticks = '',
                showline = False,
                showgrid = False,
                categoryorder = 'array',
                categoryarray = selected_term_ids,
            ),
            bgcolor = 'white'
        ),

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
        paper_bgcolor='white',
        plot_bgcolor='white',
        width = CIRCLE_WIDTH * 2,
        height = CIRCLE_WIDTH * 2,
        margin = dict(
            l = 0,
            r = 0,
            t = 0,
            b = 0
        ),
    )

    return True, fig, CIRCLE_WIDTH * 2 + 50, False, "", controls_style



@callback(
    Output("barplot_terms_modal", "opened", allow_duplicate = True),
    Output("barplot_terms", "figure", allow_duplicate = True),
    Output("barplot_terms_modal", "size", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("barplot_controls", "style", allow_duplicate = True),
    Input("histogram_results", "n_clicks"),
    State("sessionid", "children"),
    State("barplot_controls", "style"),
    prevent_initial_call = True,
)
def open_histogram(n_clicks, session_id, controls_style):

    if session_id not in sessions:
        return (
            False,
            no_update,
            no_update,
            True,
            "Your session has expired. Please refresh the website.",
            no_update,
        )

    sessions[session_id].time = time.time()
    results = sessions[session_id].result
    controls_style["display"] = "none"

    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x = [r.pvalue for r in results],
            nbinsx = 20,
            marker = dict(color = "#1f77b4"),
        )
    )
    fig.update_layout(
        title = 'P-value Histogram',
        xaxis_title = 'Uncorrected p-values',
        yaxis_title = 'Count',
    )
    return True, fig, "90%", False, "", controls_style



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






@callback(
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Output("term_path_area", "children", allow_duplicate = True),
    Input("term_lipids_modal_grid", "cellRendererData"),
    Input("term_proteins_modal_grid", "cellRendererData"),
    Input("term_metabolites_modal_grid", "cellRendererData"),
    State("sessionid", "children"),
    State("term_molecules_modal_id", "children"),
    State("term_lipids_modal_grid", "rowData"),
    State("term_proteins_modal_grid", "rowData"),
    State("term_metabolites_modal_grid", "rowData"),
    State("select_organism", "value"),
    prevent_initial_call = True,
)
def show_molecule_term_path(
    renderer_data_lipids,
    renderer_data_proteins,
    renderer_data_metabolites,
    session_id,
    target_term_id,
    row_data_lipids,
    row_data_proteins,
    row_data_metabolites,
    organism,
):
    if session_id not in sessions:
        return True, "Your session has expired. Please refresh the website."

    if target_term_id not in sessions[session_id].search_terms:
        return True, "Your session has expired. Please refresh the website."

    trigger = callback_context.triggered[0]["prop_id"].split(".")[0]
    if trigger == "term_lipids_modal_grid":
        row_index = renderer_data_lipids["rowIndex"]
        if len(row_data_lipids) <= row_index:
            return True, "Your session has expired. Please refresh the website."
        molecule = row_data_lipids[row_index]["molecule"]

    elif trigger == "term_proteins_modal_grid":
        row_index = renderer_data_proteins["rowIndex"]
        if len(row_data_proteins) <= row_index:
            return True, "Your session has expired. Please refresh the website."
        molecule = "UNIPROT:" + row_data_proteins[row_index]["molecule"]

    elif trigger == "term_metabolites_modal_grid":
        row_index = renderer_data_metabolites["rowIndex"]
        if len(row_data_metabolites) <= row_index:
            return True, "Your session has expired. Please refresh the website."
        molecule = row_data_metabolites[row_index]["molecule"]

    ontology = enrichment_ontologies[organism]
    term_path = []
    for i, term_id in enumerate(sessions[session_id].search_terms[target_term_id][molecule]):
        if i > 0: term_path.append(dmc.Text("▼", style = {"textAlign": "center"}))
        href, term_name = ".", ""
        if term_id in ontology.ontology_terms:
            term_name = ontology.ontology_terms[term_id].name
            if term_id[:3] == "GO:":
                href = "https://amigo.geneontology.org/amigo/term/" + term_id

            elif term_id[:3] == "SMP":
                href = "https://pathbank.org/view/" + term_id

            elif term_id[:5] == "LION:" or term_id[:4] == "CAT:":
                href = "https://bioportal.bioontology.org/ontologies/LION?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F" + term_id.replace(":", "_")

            elif term_id[:5] == "RHEA:":
                href = f"https://www.rhea-db.org/rhea/{term_id[5:]}"

            elif term_id[:8] == "UNIPROT:":
                href = f"https://www.uniprot.org/uniprotkb/{term_id[8:]}/entry"

            elif term_id[:6] == "CHEBI:":
                href = f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={term_id}"

            elif term_id[:5] == "DOID:":
                href = f"https://disease-ontology.org/?id={term_id}"

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


@callback(
    Output("use_bounded_fatty_acyls", "checked", allow_duplicate = True),
    Output("info_modal", "opened", allow_duplicate = True),
    Output("info_modal_message", "children", allow_duplicate = True),
    Input("use_bounded_fatty_acyls", "checked"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def switch_bounded_fatty_acyls(checked, session_id):
    if checked == None:
        raise exceptions.PreventUpdate

    if session_id not in sessions:
        return dash.no_update, do_activate_alert, "Your session has expired. Please refresh the website."

    session = sessions[session_id]
    session.use_bounded_fatty_acyls = checked
    session.data_loaded = False

    raise exceptions.PreventUpdate



