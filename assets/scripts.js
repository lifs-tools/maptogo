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
    var term_id = data["termid"];
    if (term_id.startsWith("GO:")){
        href = "https://amigo.geneontology.org/amigo/term/" + term_id;
    }
    else if (term_id.startsWith("SMP")){
        href = "https://pathbank.org/view/" + term_id;
    }
    else if (term_id.startsWith("LION:") || term_id.startsWith("CAT:")){
        href = "https://bioportal.bioontology.org/ontologies/LION?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F" + term_id.replace(":", "_");
    }

    if (href !== ""){
        return React.createElement(
            "a",
            {
                "href": href,
                "style": {"width": "100%", "text-align": "center", "color": LINK_COLOR},
                "target": "_blank",
            },
            term_id
        );
    }
    return React.createElement(
        "div",
        {},
        term_id
    );
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

