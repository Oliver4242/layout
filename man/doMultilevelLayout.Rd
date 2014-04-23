\name{doMultilevelLayout}
\alias{doMultilevelLayout}
\title{Creates a hierarchy of layouted graphs}
\usage{
doMultilevelLayout(g, file = NULL, verb = TRUE, iter = 1000,
  RMSMIN = NA)
}
\arguments{
  \item{g}{the graph}

  \item{verb}{if \code{TRUE} prints out debug information}

  \item{iter}{the number of iterations of the inner
  algorithm}

  \item{RMSMIN}{stops the layout in a certain levelif

  Hallo Gallo}
}
\description{
Creates a hierarchy of layouted graphs. If file is not NULL
pngs are created during the layout-process If RMS is set,
the layout stops and pops up to the next level if RMS is
less than RMSMIN
}

