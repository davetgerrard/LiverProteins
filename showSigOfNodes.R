
#########################################
#					#
#	Dave Gerrard 			#
#	University of Manchester	#
#	2011				#
#					#
#########################################

function (GOdata, termsP.value, firstSigNodes = 10, reverse = TRUE, 
    sigForAll = TRUE, wantedNodes = NULL, putWN = TRUE, putCL = 0, 
    type = NULL, showEdges = TRUE, swPlot = TRUE, useFullNames = TRUE, 
    oldSigNodes = NULL, useInfo = c("none", "pval", "counts", 
        "def", "np", "all")[1], plotFunction = GOplot, .NO.CHAR = 20) 
{
    require("Rgraphviz") || stop("package Rgraphviz is required")
    if (!is.null(firstSigNodes)) 
        sigTerms <- sort(termsP.value)[1:firstSigNodes]
    else sigTerms <- numeric(0)
    if (putWN && !is.null(wantedNodes)) 
        baseNodes <- union(names(sigTerms), wantedNodes)
    else baseNodes <- names(sigTerms)
    if (length(baseNodes) == 0) 
        stop("No nodes were selected")
    if (putCL) {
        goDAG.r2l <- reverseArch(graph(GOdata))
        for (i in 1:putCL) {
            newNodes <- unique(unlist(adj(goDAG.r2l, baseNodes)))
            baseNodes <- union(newNodes, baseNodes)
        }
    }
    dag <- inducedGraph(graph(GOdata), baseNodes)
    if (reverse) 
        dag <- reverseArch(dag)
    termCounts <- termStat(GOdata, nodes(dag))
    if (!is.null(type)) {
        if (swPlot) 
            GOplot.counts(dag, wantedNodes = wantedNodes, nodeCounts = termCounts, 
                showEdges = showEdges)
        return(dag)
    }
    pval.info <- function(whichNodes) {
        ret.val <- format.pval(termsP.value[whichNodes], digits = 3, 
            eps = 1e-20)
        names(ret.val) <- whichNodes
        return(ret.val)
    }
    .pval = pval.info(nodes(dag))
    .def = .getTermsDefinition(whichTerms = nodes(dag), ontology(GOdata), 
        numChar = .NO.CHAR)
    .counts = apply(termCounts[, c("Significant", "Annotated")], 
        1, paste, collapse = "/")
    nodeInfo <- switch(useInfo, none = NULL, pval = .pval, def = .def, 
        counts = .counts, np = paste(.def, .pval, sep = "\\\n"), 
        all = paste(.def, .pval, .counts, sep = "\\\n"))
    if (sigForAll) 
        sigNodes <- termsP.value[nodes(dag)]
    else sigNodes <- sigTerms
    if (is.null(wantedNodes)) 
        wantedNodes <- names(sigTerms)
    complete.dag <- plotFunction(dag, sigNodes = sigNodes, genNodes = names(sigTerms), 
        wantedNodes = wantedNodes, showEdges = showEdges, useFullNames = useFullNames, 
        oldSigNodes = oldSigNodes, nodeInfo = nodeInfo)
    if (swPlot && !is.null(complete.dag)) 
        plot(complete.dag)
    return(list(dag = dag, complete.dag = complete.dag))
}
<environment: namespace:topGO>