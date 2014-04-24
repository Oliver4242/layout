#################################################
#
#' Wrapper for prefuse library 
#'
#' 
#' @param g the igraph graph object
#' @param verb if \code{TRUE} prints out debug information
#' @param iter the number of iterations of the inner algorithm
#' @return A layouts 
#' 
#' @references
#' \url{http://arxiv.org/abs/1207.6282}
#' 
#' @author Oliver D\"urr \email{oliver.duerr@@zhaw.ch}

#' 
#' @examples
#' 
#' # Simple force-directed layout 
#' g <- barabasi.game(100, directed=FALSE)
#' l <- layout.prefuse(g, iter=100)
#' # plot(g, layout=l)
layout.prefuse <- function(graph, pos=NULL, iter, verb=TRUE, graConst=-5e-6, steps=10, parameters="-Xmx1500m"){
  el <- get.edgelist(graph)
  nV <- vcount(graph)
  if (verb) print('----------------------- start lib ------------------------')
  .jinit('JavaLib.jar', force.init=FALSE, parameters=parameters) # this starts the JVM, note
  obj=.jnew("prefuseExtension/PrefuseWrapper", check=TRUE)
  if (verb) print(paste("Java Code version :", .jcall(obj, "S", method="version")))
  if (is.null(pos)) {
    print("No positions given, starting with random positions...") 
    pos = matrix(as.double(rnorm(n=2*nV, mean=0, sd=0.1)), ncol=2)
  }
  el.i <- .jarray(matrix(as.integer(el), nrow=nrow(el)),dispatch=TRUE) #TODO better way than recreation
  pos.d <- .jarray(pos,dispatch=TRUE)
  if (!is.null(E(graph)$weight) ) {
    eWeights.d <- .jarray(E(graph)$weight,dispatch=TRUE)
  } else {
    eWeights.d <- .jarray(c(0),dispatch=TRUE);
  }
  if (!is.null(V(graph)$weight) ) {
    vWeights.d <- .jarray(V(graph)$weight,dispatch=TRUE)
    if (verb) print("Using vertex weights")
  } else {
    vWeights.d <- .jarray(c(0),dispatch=TRUE);
  }
  #vWeights.d <- .jarray(c(0),dispatch=TRUE);
  #eWeights.d <- .jarray(c(0),dispatch=TRUE);
  
  if (verb) print('----------------------- Start Task Java ------------------------')
  result=.jcall(obj, returnSig="[[D", method="doLayout", el.i, pos.d, eWeights.d, vWeights.d,
                as.integer(nV), as.integer(iter), as.double(graConst), as.integer(steps))
  if (verb) print('----------------------- Finished Task Java ------------------------')
  return(t(sapply(result,.jevalArray)));
}
