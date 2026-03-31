var dagcomponentfuncs = (window.dashAgGridComponentFunctions = window.dashAgGridComponentFunctions || {});
var dagfuncs = window.dashAgGridFunctions = window.dashAgGridFunctions || {};
window.dash_clientside = window.dash_clientside || {};

var LINK_COLOR = "#2980B9";

dagcomponentfuncs.CustomNoRowsOverlay = function (props) {
    return React.createElement(
        'div',
        {
            style: {
                color: 'grey',
                padding: 10,
                fontSize: props.fontSize
            },
        },
        props.message
    );
};



function get_url(term_id){
    if (term_id.startsWith("GO:")){
        return "https://www.ebi.ac.uk/QuickGO/term/" + term_id;
    }
    else if (term_id.startsWith("SMP")){
        return "https://pathbank.org/view/" + term_id;
    }
    else if (term_id.startsWith("LION:") || term_id.startsWith("CAT:")){
        return "https://bioportal.bioontology.org/ontologies/LION?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F" + term_id.replace(":", "_");
    }
    else if (term_id.startsWith("DOID:")){
        return "https://disease-ontology.org/?id=" + term_id;
    }
    else if (term_id.startsWith("MONDO:")){
        return "https://monarchinitiative.org/" + term_id;
    }
    else if (term_id.startsWith("HGNC:")){
        return "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/" + term_id;
    }
    else if (term_id.startsWith("HP:")){
        return "https://hpo.jax.org/browse/term/" + term_id;
    }
    else if (term_id.startsWith("NCBI:")){
        return "https://www.ncbi.nlm.nih.gov/gene/" + term_id.split(":")[1];
    }
    else if (term_id.startsWith("ENS") || term_id.startsWith("WBGene") || term_id.startsWith("FBgn")){
        return "https://www.ensembl.org/id/" + term_id;
    }
    else if (term_id.startsWith("R-")){
        return "https://reactome.org/PathwayBrowser/#/" + term_id;
    }
    else if (term_id.startsWith("LM")){
        return "https://www.lipidmaps.org/databases/lmsd/" + term_id;
    }
    else if (term_id.startsWith("UNIPROT:")){
        return "https://www.uniprot.org/uniprotkb/" + term_id.split(":")[1];
    }
    else if (term_id.startsWith("CHEBI:")){
        return "https://www.ebi.ac.uk/chebi/" + term_id;
    }
    return "";
}



dagcomponentfuncs.TermIDRenderer = function (props) {
    const {setData, data} = props;

    var term_ids = data["termid"].split("|");
    var reacts = [];
    var unique_key = 0;
    for (var term_id of term_ids){
        var href = get_url(term_id);
        if (term_id.startsWith("UNIPROT:")){
            term_id = term_id.split(":")[1];
        }

        if (reacts.length > 0){
            reacts.push(
                React.createElement("span", {key: "k" + unique_key}, " | ")
            );
            unique_key += 1;
        }

        if (href !== ""){
            reacts.push(React.createElement(
                "a",
                {
                    "href": href,
                    "style": {"width": "100%", "text-align": "center", "color": LINK_COLOR},
                    "target": "_blank",
                    key: "k" + unique_key,
                },
                term_id
            ));
        }
        else {
            reacts.push(React.createElement(
                "span",
                {key: "k" + unique_key},
                term_id
            ));
        }
        unique_key += 1;
    }
    return React.createElement(React.Fragment, null, reacts);
}



dagcomponentfuncs.TermRenderer = function (props) {
    const {setData, data} = props;

    function onTextClicked() {
        props.api.startEditingCell({
            rowIndex: props.rowIndex,
            colKey: props.column.getId(),
        });
        setData();
    }

    return React.createElement(
        "div",
        {
            onClick: onTextClicked,
            "style": {"cursor": "pointer", "color": LINK_COLOR, "textDecorationLine": "underline"},
        },
        data["term"]
    );
}



dagcomponentfuncs.MoleculeRenderer = function (props) {
    const {setData, data} = props;

    function onTextClicked() {
        props.api.startEditingCell({
            rowIndex: props.rowIndex,
            colKey: props.column.getId(),
        });
        setData();
    }


    return React.createElement(
        "div",
        {
            onClick: onTextClicked,
            "style": {"cursor": "pointer", "color": LINK_COLOR, "textDecorationLine": "underline"},
        },
        data["molecule"]
    );
}


document.addEventListener("DOMContentLoaded", function () {
    // Wait for Dash to render the plot
    const observer = new MutationObserver(function () {
        const plot = document.querySelector("#barplot_terms .js-plotly-plot");
        if (plot && !plot.__sunburstClickHooked) {
            plot.on('plotly_sunburstclick', function (event) {
                const pt = event.points?.[0];
                if (pt) {
                    for (var term_id of pt.id.split("|")){
                        var url = get_url(term_id);
                        if (url.length > 0) window.open(url, '_blank');
                    }
                }

                // Prevent zoom/expansion
                return false;
            });
            plot.__sunburstClickHooked = true;
        }
    });

    observer.observe(document.body, { childList: true, subtree: true });
});

