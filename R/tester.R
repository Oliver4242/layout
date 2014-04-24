testPrefuse <- function() {
  g <- barabasi.game(500, directed=FALSE)
  l <- layout.prefuse(g, iter=1000, verb=FALSE)
  wt <- multilevel.community(g)
  V(g)$color <- wt$membership
  plot(g, layout=l, vertex.size=3, vertex.label=NA)
} 

testML <- function() {
  g <- barabasi.game(500, directed=FALSE)
  res <- doMultilevelLayout(g, iter=1000, verb=FALSE, method="prefuse")
  l <- res$layouts[[1]]
  wt <- multilevel.community(g)
  V(g)$color <- wt$membership
  plot(g, layout=l, vertex.size=3, vertex.label=NA)
} 
