var dagcomponentfuncs = (window.dashAgGridComponentFunctions = window.dashAgGridComponentFunctions || {});
var dagfuncs = window.dashAgGridFunctions = window.dashAgGridFunctions || {};
window.dash_clientside = window.dash_clientside || {};

dagcomponentfuncs.TermIDRenderer = function (props) {
    const {setData, data} = props;

    function onClickDel(){
        setData();
    }

    var href = "";
    var term = data["termid"];
    if (term.startsWith("GO:")){
        href = "https://amigo.geneontology.org/amigo/term/" + term;
    }
    else if (term.startsWith("SMP")){
        href = "https://pathbank.org/view/" + term;
    }
    else if (term.startsWith("LION:") || term.startsWith("CAT:")){
        href = "https://bioportal.bioontology.org/ontologies/LION?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2F" + term.replace(":", "_");
    }

    if (href !== ""){
        return React.createElement(
            "a",
            {
                "href": href,
                "style": {"width": "100%", "text-align": "center", "color": "#2980B9"},
                "target": "_blank",
            },
            term
        );
    }
    return React.createElement(
        "div",
        {},
        term
    );
}
