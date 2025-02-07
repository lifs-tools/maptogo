from dash import Dash, dcc, html, Input, Output, State, callback, exceptions, no_update, MATCH, ALL, callback_context
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
from EnrichmentDataStructure import EnrichmentOntology
from pygoslin.domain.LipidFaBondType import LipidFaBondType
from pygoslin.domain.LipidLevel import LipidLevel
from pygoslin.parser.Parser import LipidParser
import logging
import pathlib
import time
import hashlib
import threading


hash_function = hashlib.new('sha256')
current_path = pathlib.Path(__file__).parent.resolve()
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

sessions, examples = {}, {}
xl = pd.ExcelFile(f"{current_path}/Data/examples.xlsx")
first_example = ""
for i, worksheet_name in enumerate(xl.sheet_names):
    if i == 0: first_example = worksheet_name
    df = xl.parse(worksheet_name)
    worksheet = {"bg": [], "reg": [], "title": "", "desc": "", "comp": "", "doi": ""}
    if "Background" in df: worksheet["bg"] = list(df["Background"].dropna())
    if "Regulated" in df: worksheet["reg"] = list(df["Regulated"].dropna())
    if "Title" in df: worksheet["title"] = df.loc[0, "Title"]
    if "Description" in df: worksheet["desc"] = df.loc[0, "Description"]
    if "Comparison" in df: worksheet["comp"] = df.loc[0, "Comparison"]
    if "DOI" in df: worksheet["doi"] = df.loc[0, "DOI"]
    if type(worksheet["title"]) in {float, np.float64} and np.isnan(worksheet["title"]): worksheet["title"] = ""
    if type(worksheet["desc"]) in {float, np.float64} and np.isnan(worksheet["desc"]): worksheet["desc"] = ""
    if type(worksheet["comp"]) in {float, np.float64} and np.isnan(worksheet["comp"]): worksheet["comp"] = ""
    if type(worksheet["doi"]) in {float, np.float64} and np.isnan(worksheet["doi"]): worksheet["doi"] = ""
    examples[worksheet_name] = worksheet



example_options = html.Div(
    children = [
        dmc.Group(
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
        ) for i, (example_name, example) in enumerate(examples.items())
    ]
)

# Create the Dash app
app = Dash("app", update_title = None)
app.title = "GO lipids"

lipid_parser = LipidParser()
enrichment_ontologies = {}
for tax_name, tax_id in species.items():
    logger.info(f"loading {tax_name}")
    enrichment_ontologies[tax_id] = EnrichmentOntology(f"{current_path}/Data/ontology_{tax_id}.gz", lipid_parser = lipid_parser)


#enrichment_ontologies[INIT_ORGANISM].set_background_lipids(["LPA 18:0", "LPA 20:0"])
#enrichment_ontologies[INIT_ORGANISM].enrichment_analysis(["LPA 18:0"], ["Biological process"])




def session_timer_trigger(sessions):
    while True:
        time.sleep(10)
        current_time = time.time()
        sessions_to_delete = []
        for session_id, session_data in sessions.items():
            if "time" not in session_data or current_time - session_data["time"] > 60 * 60: # one hour
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
    sessions[session_id] = {"time": current_time, "data": None}

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
            html.Div(
                html.A(
                    dmc.Text(
                        "Description & disclaimer",
                    ),
                    id = "disclaimer",
                    style = {"color": "#2980B9", "cursor": "pointer"},
                ),
                style = {"height": "100%", "display": "flex", "alignItems": "flex-start", "justifyContent": "right"},
            ),
        ], cols = 2),
        html.Div(
            id = "background_lipids",
            style = {"display": "none"},
        ),
        html.Div(
            id = "regulated_lipids",
            style = {"display": "none"},
        ),
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
                            "Close",
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
            id = "term_lipids_modal",
            zIndex = 10000,
            children = [
                dag.AgGrid(
                    id = "term_lipids_modal_grid",
                    columnDefs = [
                        {
                            'field': "lipid",
                            "headerName": "Lipid",
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
                    },
                ),
            ],
            size = "40%",
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
            children = [
                dmc.SimpleGrid(
                    cols = 2,
                    children = [
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
                    ],
                ),
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
                                ),
                                dmc.ActionIcon(
                                    DashIconify(icon="material-symbols:download-rounded", width = 20),
                                    id = "icon_download_results",
                                    title = "Download table",
                                ),
                            ], spacing = 0),
                        ),
                        style = {"height": "100%", "display": "flex", "alignItems": "flex-end", "justifyContent": "right"},
                    ),
                ], cols = 2),
            ],
        ),
        dmc.SimpleGrid(
            cols = 2,
            children = [
                html.Div([
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
                        dmc.Button(
                            "Load example datasets",
                            id = "load_examples_button",
                        ),
                    ], cols = 2),
                    dmc.Title(
                        "Analysis parameters",
                        order = 5,
                        style = {"marginTop": "20px"},
                    ),
                    html.Div(
                        [
                            dmc.Alert(
                                "Something happened! You made a mistake and there is no going back, your data was lost forever!",
                                id = "alert_enrichment",
                                title = "Enrichment analysis error!",
                                color = "red",
                                withCloseButton = True,
                                hide = True,
                                style = {"marginTop": "10px", "marginBottom": "10px"},
                            ),
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
                                    html.Div(
                                        dmc.Switch(
                                            id = "switch_ignore_unparsable_lipids",
                                            checked = False,
                                            label = "Ignore unrecognizable lipid names",
                                            style = {"paddingBottom": "8px"},
                                        ),
                                        style = {"height": "100%", "display": "flex", "alignItems": "flex-end", "paddingLeft": "10px"},
                                    ),
                                    dmc.MultiSelect(
                                        id = "select_domains",
                                        data = sorted(list(enrichment_ontologies[INIT_ORGANISM].domains)),
                                        value = ["Biological process"],
                                        label = "Select domain(s):",
                                    ),
                                    html.Div(
                                        dmc.Switch(
                                            id = "switch_ignore_unknown_regulated_lipids",
                                            checked = False,
                                            label = "Ignore regulated lipids that aren't in background",
                                            style = {"paddingBottom": "8px"},
                                        ),
                                        style = {"height": "100%", "display": "flex", "alignItems": "flex-end", "paddingLeft": "10px"},
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
                                        dmc.Button(
                                            "Run enrichment analysis",
                                            id = "button_run_enrichment",
                                        ),
                                        style = {"height": "100%", "display": "flex", "alignItems": "flex-end", "paddingLeft": "10px"},
                                    ),
                                ],
                            ),
                        ],
                        style = {"border": "1px solid #dddddd", "radius": "10px", "padding": "10px"},
                    ),
                ]),
                html.Div([
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
    Output("graph_enrichment_results", "selectedRows", allow_duplicate = True),
    Output("alert_enrichment", "hide", allow_duplicate = True),
    Output("alert_enrichment", "children", allow_duplicate = True),
    Output("background_lipids", "children", allow_duplicate = True),
    Output("regulated_lipids", "children", allow_duplicate = True),
    Input("button_run_enrichment", "n_clicks"),
    State("textarea_all_lipids", "value"),
    State("textarea_regulated_lipids", "value"),
    State("select_organism", "value"),
    State("select_domains", "value"),
    State("switch_ignore_unparsable_lipids", "checked"),
    State("select_test_method", "value"),
    State("switch_ignore_unknown_regulated_lipids", "checked"),
    State("select_term_regulation", "value"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def run_enrichment(
    n_clicks,
    all_lipids_list,
    regulated_lipids_list,
    organism,
    domains,
    ignore_unrecognizable_lipids,
    correction_method,
    ignore_unknown,
    term_regulation,
    session_id,
):
    if len(all_lipids_list) == 0:
        return "", [], False, "Please paste lipid names into the first text area.", "", ""

    if len(regulated_lipids_list) == 0:
        return "", [], False, "Please paste lipid names into the second text area."

    if len(domains) == 0:
        return "", [], False, "No domain(s) selected.", "", ""

    lipidome, regulated_lipids = set(), set()
    for lipid_name in all_lipids_list.split("\n"):
        if len(lipid_name) == 0: continue
        try:
            lipidome.add(lipid_parser.parse(lipid_name).get_lipid_string())
        except Exception as e:
            if not ignore_unrecognizable_lipids:
                return "", [], False, f"Lipid name '{lipid_name}' unrecognizable! Maybe enable the 'Ignore unrecognizable lipid names' option.", "", ""

    for lipid_name in regulated_lipids_list.split("\n"):
        if len(lipid_name) == 0: continue
        try:
            regulated_lipids.add(lipid_parser.parse(lipid_name).get_lipid_string())
        except Exception as e:
            if not ignore_unrecognizable_lipids:
                return "", [], False, f"Lipid name '{lipid_name}' unrecognizable! Maybe enable the 'Ignore unrecognizable lipid names' option.", "", ""

    if len(lipidome) == 0:
        return "", [], False, "No background lipid left after lipid recognition.", "", ""

    if len(regulated_lipids) == 0:
        return "", [], False, "No regulated lipid left after lipid recognition.", "", ""

    if len(regulated_lipids) > len(lipidome):
        return "", [], False, "Length of regulated lipid list must be smaller than background list.", "", ""

    left_lipids = regulated_lipids - lipidome
    if len(left_lipids) > 0:
        if ignore_unknown:
            lipidome -= left_lipids
        else:
            return "", [], False, "The lipid" + (' ' if len(left_lipids) == 1 else 's ') + "'" + "', '".join(left_lipids) + ("' does" if len(left_lipids) == 1 else "' do") + " not occur in the background list. Maybe enable the 'Ignore regulated lipids that aren't in background' option.", "", ""

    ontology = enrichment_ontologies[organism]
    ontology.set_background_lipids(lipidome)

    results = ontology.enrichment_analysis(regulated_lipids, domains, term_regulation)
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

    sessions[session_id] = {"time": time.time(), "data": {result.term.term_id: result for result in results}}

    return (
        "",
        data,
        data,
        True,
        "",
        "|".join(lipidome),
        "|".join(regulated_lipids),
    )



@callback(
    Output("loading_output", "children", allow_duplicate = True),
    Output("download_data", "data", allow_duplicate = True),
    Input("icon_download_results", "n_clicks"),
    State("graph_enrichment_results", "virtualRowData"),
    State("graph_enrichment_results", "selectedRows"),
    State("background_lipids", "children"),
    State("regulated_lipids", "children"),
    State("sessionid", "children"),
    prevent_initial_call = True,
)
def download_table(
    _,
    graph_enrichment_results,
    selected_rows,
    background_lipids,
    regulated_lipids,
    session_id
):
    if session_id not in sessions or "data" not in sessions[session_id]:
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

    background_lipids = sorted(background_lipids.split("|"))
    regulated_lipids = sorted(regulated_lipids.split("|"))
    regulated_lipid_set = set(regulated_lipids)
    associated_lipids = {}
    data = sessions[session_id]["data"]

    for term_id, term in zip(term_ids, terms):
        if term_id not in data: continue
        lipids = sorted(list(data[term_id].term.term_paths.keys()))
        associated_lipids[term] = lipids
        associated_lipids[f"{term}, regulated"] = ["X" if lipid in regulated_lipid_set else "" for lipid in lipids]

    output = io.BytesIO()
    writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
    df.to_excel(writer, sheet_name = "GO lipidomics results", index = False)
    pd.DataFrame.from_dict(associated_lipids, orient = "index").transpose().to_excel(writer, sheet_name = "Associated lipids", index = False)
    pd.DataFrame({"LipidName": background_lipids}).to_excel(writer, sheet_name = "Background lipids", index = False)
    pd.DataFrame({"LipidName": regulated_lipids}).to_excel(writer, sheet_name = "Regulated lipids", index = False)
    writer._save()

    return "", dcc.send_bytes(output.getvalue(), "GO_lipidomics_results.xlsx")



@callback(
    Output("load_examples_modal", "opened", allow_duplicate = True),
    Input("load_examples_button", "n_clicks"),
    prevent_initial_call = True,
)
def open_load_examples_modal(n_clicks):
    return True



@callback(
    Output("load_examples_modal", "opened", allow_duplicate = True),
    Output("textarea_all_lipids", "value", allow_duplicate = True),
    Output("textarea_regulated_lipids", "value", allow_duplicate = True),
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
        "\n".join(examples[index]["bg"]),
        "\n".join(examples[index]["reg"]),
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
    Input("disclaimer", "n_clicks"),
    prevent_initial_call = True,
)
def disclaimer_clicked(n_clicks):
    return True



@callback(
    Output("term_lipids_modal", "opened", allow_duplicate = True),
    Output("term_lipids_modal", "title", allow_duplicate = True),
    Output("term_lipids_modal_grid", "rowData", allow_duplicate = True),
    Input("graph_enrichment_results", "cellRendererData"),
    State("sessionid", "children"),
    State("regulated_lipids", "children"),
    prevent_initial_call = True,
)
def open_term_window(row_data, session_id, regulated_lipids):
    if session_id not in sessions or "rowId" not in row_data:
        raise exceptions.PreventUpdate

    term_id = row_data["rowId"]
    if "data" not in sessions[session_id] or term_id not in sessions[session_id]["data"]:
        raise exceptions.PreventUpdate

    result = sessions[session_id]["data"][term_id]
    for k, v in result.term.term_paths.items():
        print(k, v)
    lipids = sorted(list(result.term.term_paths.keys()))
    regulated_lipids = set(regulated_lipids.split("|"))

    return (
        True,
        f"Lipids for '{result.term.name}'",
        [{"lipid": lipid, "regulated": ("X" if lipid in regulated_lipids else "")} for lipid in lipids],
    )




