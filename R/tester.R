testPrefuse <- function() {
  require(igraph)
  g <- barabasi.game(500, directed=FALSE)
  l <- layout.prefuse(g, iter=1000, verb=TRUE)
  wt <- multilevel.community(g)
  V(g)$color <- wt$membership
  plot(g, layout=l, vertex.size=3, vertex.label=NA)
} 

testMLPrefuse <- function() {
  require(igraph)
  g <- barabasi.game(500, directed=FALSE)
  l <- layout.multi.level(graph=g, iter=1000, verb=TRUE, method="prefuse")
  #l <- layout.fruchterman.reingold(graph=g)
  wt <- multilevel.community(g)
  V(g)$color <- wt$membership
  plot(g, layout=l, vertex.size=3, vertex.label=NA)
} 

testMLNonPrefuse <- function() {
  g <- barabasi.game(500, directed=FALSE)
  l <- layout.multi.level(graph=g, iter=1000, verb=FALSE, method="fr")
  wt <- multilevel.community(g)
  V(g)$color <- wt$membership
  plot(g, layout=l, vertex.size=3, vertex.label=NA)
} 

testWrongMethod <- function() {
  g <- barabasi.game(500, directed=FALSE)
  l <- layout.multi.level(graph=g, iter=1000, verb=FALSE, method="egal")
  wt <- multilevel.community(g)
  V(g)$color <- wt$membership
  plot(g, layout=l, vertex.size=3, vertex.label=NA)
} 
