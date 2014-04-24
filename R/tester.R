testPrefuse <- function() {
  g <- barabasi.game(10000, directed=FALSE)
  l <- layout.prefuse(g, iter=100)
  plot(g, layout=l)
} 