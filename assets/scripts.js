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


dagcomponentfuncs.TermIDRenderer = function (props) {
    const {setData, data} = props;

    var href = "";
    var term_ids = data["termid"].split("|");
    var reacts = [];
    var unique_key = 0;
    for (var term_id of term_ids){
        if (term_id.startsWith("GO:")){
            href = "https://amigo.geneontology.org/amigo/term/" + term_id;
        }
        else if (term_id.startsWith("SMP")){
            href = "https://pathbank.org/view/" + term_id;
        }
        else if (term_id.startsWith("LION:") || term_id.startsWith("CAT:")){
            href = "https://bioportal.bioontology.org/ontologies/LION?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F" + term_id.replace(":", "_");
        }
        else if (term_id.startsWith("DOID:")){
            href = "https://disease-ontology.org/?id=" + term_id;
        }
        else if (term_id.startsWith("MONDO:")){
            href = "https://monarchinitiative.org/" + term_id;
        }
        else if (term_id.startsWith("HGNC:")){
            href = "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/" + term_id;
        }
        else if (term_id.startsWith("HP:")){
            href = "https://hpo.jax.org/browse/term/" + term_id;
        }
        else if (term_id.startsWith("NCBI:")){
            href = "https://www.ncbi.nlm.nih.gov/gene/" + term_id.split(":")[1];
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

