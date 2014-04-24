\name{doMultilevelLayout}
\alias{doMultilevelLayout}
\title{Creates a hierarchy of layouted graphs}
\usage{
doMultilevelLayout(g, file = NULL, verb = FALSE, iter = 1000,
  RMSMIN = NA)
}
\arguments{
  \item{g}{the igraph graph object}

  \item{verb}{if \code{TRUE} prints out debug information}

  \item{iter}{the number of iterations of the inner
  algorithm}

  \item{RMSMIN}{stops the layout in a certain level if the
  RMS differences are below the threshold given by
  \code{RMSMIN}}

  \item{file}{if not \code{NULL} a pngs of are created
  during the layout-process.}
}
\value{
A result object including the graphs and the layouts for
the different levels
}
\description{
Creates a 2-dimensional layout (placement of the nodes) of
a graph using a multilevel approach (see details below).
}
\details{
In the first stage the graph is decomposed into a series of
smaller graphs using the a multilevel community detection
algorithm \code{\link{multilevel.community}}. Then the a
layout is done for the lowest level (usally the smallest
graph).
}
\examples{
# Simple example of a barabasi game
g <- barabasi.game(100, directed=FALSE)
res <- doMultilevelLayout(g)
plot(g, layout=res$layouts[[1]])
#TODO create example with a valiud file
}
\author{
Oliver D\"urr \email{oliver.duerr@zhaw.ch}
}
\references{
\url{http://arxiv.org/abs/1207.6282}
}

