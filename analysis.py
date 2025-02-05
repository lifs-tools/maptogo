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

current_path = pathlib.Path(__file__).parent.resolve()

logger = logging.getLogger(__name__)
logging.basicConfig(
    level = os.environ.get("LOGLEVEL", "INFO"),
    format = "%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
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

examples = {}
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
                                    style = {"color": "#2980B9"},
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
logger.info("loaded")



app.layout = html.Div([
    dcc.Download(id = "download_data"),
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
    html.Div(
        html.H1("Gene Ontology (GO) lipidomics enrichment analysis"),
        className = "mantine-InputWrapper-label mantine-MultiSelect-label mantine-ittua2",
    ),
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
    dmc.SimpleGrid(
        cols = 2,
        children = [
            html.Div(
                dmc.SimpleGrid(
                    cols = 2,
                    children = [
                        html.Div(
                            [
                                html.H4("All lipid names in experiment (background)", style = {"marginBottom": "2px"}),
                            ],
                            className = "mantine-InputWrapper-label mantine-MultiSelect-label mantine-ittua2",
                        ),
                        html.Div(
                            [
                                html.H4("All regulated lipid names in experiment", style = {"marginBottom": "2px"}),
                            ],
                            className = "mantine-InputWrapper-label mantine-MultiSelect-label mantine-ittua2",
                        ),
                    ],
                ),
            ),
            html.Div(
                dmc.SimpleGrid([
                    html.H3("Results", style = {"marginBottom": "2px"}),
                    html.Div(
                        html.Span(
                            dmc.ActionIcon(
                                DashIconify(icon="material-symbols:download-rounded", width = 20),
                                id = "icon_download_results",
                                title = "Download table",
                            ),
                        ),
                        style = {"height": "100%", "display": "flex", "alignItems": "flex-end", "justifyContent": "right"},
                    ),
                ], cols = 2),
            ),
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
                        minRows = 20,
                    ),
                    dmc.Textarea(
                        id = "textarea_regulated_lipids",
                        minRows = 20,
                    ),
                    dmc.Button(
                        "Load example datasets",
                        id = "load_examples_button",
                    ),
                ], cols = 2),
                html.Div(style = {"height": "20px"}),
                html.H3("Analysis parameters", style = {"marginBottom": "2px"}),
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
                                        #{"value": "9606", "label": "Homo Sapiens (hsa)"},
                                        #{"value": "10090", "label": "Mus Musculus (mmu)"},
                                    ],
                                    value = "9606",
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
                                    data = sorted(list(enrichment_ontologies["9606"].domains)),
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
                        },
                        {
                            'field': "termid",
                            "headerName": "Term ID",
                            "cellRenderer": "TermIDRenderer",
                        },
                        {
                            'field': "term",
                            "headerName": "Term",
                        },
                        {
                            'field': "pvalue",
                            "headerName": "pValue",
                        },
                    ],
                    rowData = [],
                    columnSize = "responsiveSizeToFit",
                    defaultColDef = {
                        "resizable": True,
                        "suppressMovable": True,
                    },
                    style = {
                        "height": "100%",
                    },
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




@callback(
    Output("loading_output", "children"),
    Output("graph_enrichment_results", "rowData"),
    Output("alert_enrichment", "hide"),
    Output("alert_enrichment", "children"),
    Output("background_lipids", "children"),
    Output("regulated_lipids", "children"),
    Input("button_run_enrichment", "n_clicks"),
    State("textarea_all_lipids", "value"),
    State("textarea_regulated_lipids", "value"),
    State("select_organism", "value"),
    State("select_domains", "value"),
    State("switch_ignore_unparsable_lipids", "checked"),
    State("select_test_method", "value"),
    State("switch_ignore_unknown_regulated_lipids", "checked"),
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

    results = ontology.enrichment_analysis(regulated_lipids, domains)
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

    return "", data, True, "", "|".join(lipidome), "|".join(regulated_lipids)



@callback(
    Output("download_data", "data"),
    Input("icon_download_results", "n_clicks"),
    State("graph_enrichment_results", "rowData"),
    State("background_lipids", "children"),
    State("regulated_lipids", "children"),
    prevent_initial_call = True,
)
def download_table(_, graph_enrichment_results, background_lipids, regulated_lipids):
    domains = []
    term_ids = []
    terms = []
    pvalues = []

    if len(graph_enrichment_results) == 0:
        raise exceptions.PreventUpdate

    for row in graph_enrichment_results:
        domains.append(row["domain"])
        term_ids.append(row["termid"])
        terms.append(row["term"])
        pvalues.append(row["pvalue"])

    df = pd.DataFrame({"Domain": domains, "Term ID": term_ids, "Term": terms, "pValue": pvalues})

    background_lipids = sorted(background_lipids.split("|"))
    regulated_lipids = sorted(regulated_lipids.split("|"))

    output = io.BytesIO()
    writer = pd.ExcelWriter(output, engine = 'xlsxwriter')
    df.to_excel(writer, sheet_name = "GO lipidomics results", index = False)
    pd.DataFrame({"LipidName": background_lipids}).to_excel(writer, sheet_name = "Background lipids", index = False)
    pd.DataFrame({"LipidName": regulated_lipids}).to_excel(writer, sheet_name = "Regulated lipids", index = False)
    writer._save()

    return dcc.send_bytes(output.getvalue(), "GO_lipidomics_results.xlsx")



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
    Output({"type": "checkbox_type", "index": ALL}, "checked"),
    Input({"type": "checkbox_type", "index": ALL}, "checked"),
    State({"type": "checkbox_type", "index": ALL}, "id"),
    prevent_initial_call = True,
)
def checkbox_checks(_, checkbox_ids):
    index = json.loads(callback_context.triggered[0]["prop_id"].split(".")[0])["index"]
    return [checkbox_id["index"] == index for checkbox_id in checkbox_ids]
