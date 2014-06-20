HEIGHT = 900
WIDTH  = 1600 

#################################################
#
#' Creates a hierarchy of layouted graphs
#'
#' Creates a 2-dimensional layout (placement of the nodes) of a graph using a multilevel approach (see details below).
#' 
#' @param g the igraph graph object
#' @param verb if \code{TRUE} prints out debug information
#' @param iter the number of iterations of the inner algorithm
#' @param RMSMIN stops the layout in a certain level if the RMS differences are below the threshold given by \code{RMSMIN}
#' @param file if not \code{NULL} a pngs of are created during the layout-process.
#' @param the method using the inner layout either "prefuse" or "fr" 
#' @param mlInfo if true information such as the intermediate levels of the layout are returned
#' @return A result object including the graphs and the layouts for the different levels
#' 
#' @details In the first stage the graph is decomposed into a series of smaller graphs using the
#' a multilevel community detection algorithm \code{\link{multilevel.community}}. Then the a layout is
#' done for the lowest level (usally the smallest graph).  
#' 
#' @references
#' \url{http://arxiv.org/abs/1207.6282}
#' 
#' @author Oliver D\"urr \email{oliver.duerr@@zhaw.ch} 
 
#' 
#' @examples
#' # Simple example of a barabasi game
#' g <- barabasi.game(100, directed=FALSE)
#' res <- layout.multi.level(g)
#' plot(g, layout=res, vertex.size=3, vertex.label=NA)
#' layout.multi.level(g, file='tmp/dumm04%d.png')
layout.multi.level <- function(graph, file=NULL, verb=FALSE, iter=1000, RMSMIN=NA, method="prefuse", mlInfo=FALSE) {
  if (method != "prefuse" && method != "fr") {
    print(paste("Warning: invalide method ", method, " using prefuse."));
    method = "prefuse"
  }
  set.seed(1) 
  tic <- proc.time()[3]
  if (verb) print("Multi Level Layouter:: Detecting communities ...")
  wt <- multilevel.community(g, weights=NA) #TODO Check why it is not working with weight
  if (verb) print(paste("   found ", dim(wt$memberships)[1], "  levels of communities in ", round((proc.time()[3] - tic),5),  sep=""));
  graphs <- createMLGraphs(g, wt, verb);
  levels <- dim(wt$memberships)[1]
  layouts <- vector("list", levels + 1) 
  areas <- vector("list", levels + 1)
  RMSs  <- vector("list", levels + 1)
  if (levels == 0) levels = 1
  par(mfrow = c(1,1) )
  if (is.null(file) == FALSE) png(file=file, width=WIDTH, height=HEIGHT)
  l <- NULL;
  if (verb) {
    print("---------------- Creation of Layout ----------")
    print("Creating Layouts for the different levels  ...")
  }
  if (is.null(file) == FALSE) dev.off()
  for (k in (levels+1):1) {
      if (verb) {
        print(paste(" Layout iteration in level ", k, " : ", sep=""))
      }
      g.cur <- graphs[[k]]
      ret <- createNewLayoutLevel(g.cur=g.cur, members=wt$memberships, k=k-1, l=l, verb=verb, debug=TRUE, iter=iter, RMSMIN=RMSMIN, plot = !is.null(file), method=method)
      l <- ret$layout #Needed for the next iter
      layouts[[k]] <- l;
      areas[[k]] <- ret$area
      RMSs[[k]]  <- ret$RMS
  }
  if (mlInfo == FALSE) {
    return(layouts[[1]])
  } else{
    return (list(layouts = layouts, graphs=graphs, members=wt$memberships, areas = areas, RMSs = RMSs));
  }
}

# Creates a hierachy of graphs
# 1 The original
# After the first level
createMLGraphs <- function(g, wt, verb=TRUE) {
  #createNewGraphLevel.cmp <- cmpfun(createNewGraphLevel)
  levels <- dim(wt$memberships)[1]
  graphs <- vector("list", levels + 1)
  graphs[[1]] <- g
  if (levels > 0) { #R sucks sequence 1:0 bullshit
    for (i in 1:levels) {
      if (verb) print(paste(" Creating graphs for Level : ", i, sep=""))
      graphs[[i + 1]] = createNewGraphLevel(k=i, members=wt$memberships, g.cur = graphs[[i]], verb=verb)
    }
  }
  return(graphs);
}

# Caching the sums between the communities. 
# Takes the diagonal sparsely coded dgCMatrix (adjec) and converts it into a dense num.com x num.com Matrix
# containing the sums of the adjecency Matrix.
condensedSum <- function(adjec, com, num.com) {
  sums <- matrix(0, nrow=num.com, ncol=num.com)
  adjec <- as(adjec, "dgTMatrix") #adjec is a compressed sparse Matrix (dgCMatrix) a dgTMatrix is a "normal" sparse Matrix. 
  vals <- adjec@x
  ii <- adjec@i + 1
  jj <- adjec@j + 1
  for (nn in 1:length(ii)) {
    i <- ii[nn]    
    j <- jj[nn]
    comI = com[i];
    comJ = com[j];
    val = vals[nn];
    if (i <= j) {
      if  (comI > comJ) {
        sums[comI,comJ] = sums[comI,comJ] + val;
      } else{
        sums[comJ,comI] = sums[comJ,comI] + val;
      }
    }
  }
 return(sums);
}

######################################################
# Creates a new graph 
createNewGraphLevel <- function(k=1, members = wt$memberships, g.cur, verb=FALSE) {
  tic <- proc.time()[3]
  # Communities in the old graph
  com <- members[k,]        # The communities in the current level, indices refer to the orignal graph
  num.com <- max(com)       # Number of communities
  if (k > 1) {
    com <- rep(NA, max(members[k-1,]))
    for (i in 1:num.com) {
      com[unique(members[k-1, members[k,] == i])] <- i
    }
  } 
  #Creating a weighted graph for the level k
  graph.k <- graph.empty(n=num.com, directed=FALSE)
  vertex.weights <- rep(0, num.com)
  e.cur <- E(g.cur) 
  
  # Caching the edges between the communities
  el <- get.edgelist(g.cur)
  e.weights <- rep(0, num.com)
  for (i in 1:dim(el)[1]){
    if (com[el[i,1]] == com[el[i,2]]) {
      com.cur <- com[el[i,1]]
      if (k > 1 | is.weighted(g)) {  #TODO put into outer loop
        e.weights[com.cur] <- e.weights[com.cur] + e.cur[i]$weight
      } else{
        e.weights[com.cur] <- e.weights[com.cur] + 1  
      }
    }
  }  
  if (verb) print(paste("    Cummulative-Time for creating graph (subtask I)", round(proc.time()[3] - tic,5),  sep=""))
  
  # Caching the indices of the communities
  com.idxs <- list()
  for (i in 1:num.com) {
    com.idxs[[i]] <- which(com == i);
  }
  if (verb) print(paste("    Cummulative-Time for creating graph (subtask II)", round(proc.time()[3] - tic,5),  sep=""))
  
  # Caching the sums of the adjacency matrix (needs to be ported into C, Java)
  adjec <- g.cur[]; #Note adjec is a sparse matrix
  sums <- condensedSum(adjec, com, num.com) 
  #sums <- condensedSum(adjec=as.matrix(adjec), com=com, verb=verb)
  if (verb) print(paste("    Cummulative-Time for creating graph (subtask III)", round(proc.time()[3] - tic,5),  sep=""))
  if (verb) print(paste("    Number of communities ", num.com, sep=""))
  if (verb) pb <- myTxtProgressBar(min = 0.999, max = max(num.com,2), style = 3,  file=stderr())
  for (i in 1:num.com) {  
      if (verb) {setTxtProgressBar(pb, i)}
      com.i <- com.idxs[[i]]; #which(com == i);             Indices der Vertices, die zur gewuenschten Community gehoeren.
      weights  <- 0;
      if (k > 1 | is.weighted(g)) { 
        weights  <- sum(V(g.cur)[com.i]$weight)
      }
      vertex.weights[i] <- weights + 2 * e.weights[i] ### The weight of these meta-vertices is given by the sum of the weights of the previous graph plus twice the sum of the edge weights inside the community
      for (j in (i+1):num.com) { 
        if (j > num.com || i == num.com) next #Stupid R
        com.j <- com.idxs[[j]]  #com.j <- which(com == j) 
# Uncomment for testing
#         e.weight <- sum(adjec[com.j, com.i])  #dog-slow!!!!  Edges between the community i and j
#         if (abs(sums1[j, i] - sums[j, i]) > 1e-8) {
#           print(paste0(sums1[j,i], " ", sums[j, i])) #Lower triangle
#           scheisse
#         }
        e.weight <- sums[j, i]
        if(e.weight > 0) {
          graph.k <- add.edges(graph.k, c(i,j), attr=list(weight=e.weight)) 
        }
      }
    }
  #V(graph.k)$weight <- vertex.weights;
  if (verb) close(pb)
  if (verb) print(paste("    Time for creating graph ", round(proc.time()[3] - tic,5),  sep=""))
  graph.k <- set.vertex.attribute(graph=graph.k, name="weight", index=1:num.com, value=vertex.weights)
  return (graph.k);
}

ck.colors <- function(n)
  ###
  ### create a vector of n colors
  ### 20-4-2010; ch
{
  t.hues <- seq( 240, -30, length = n)/360
  t.hues[t.hues < 0] <- t.hues[t.hues < 0] + 1
  t.hues <- hsv( h = t.hues, s = 0.8, v = 1)
  return(t.hues)
}

##### Plots a graph
plotGraph <- function(g, l, com, k, rms) {
  minSize = sqrt(WIDTH * HEIGHT * 0.01 / length(l))
  V(g)$size = min(minSize, 10)
  V(g)$label.cex <- 0.001
  if (is.null(com) == FALSE) {
    V(g)$color <- ck.colors(max(com))[com]
  }
  text = paste("Level ", k, " Dim ", " x[", round(min(l[,1]),2), " ", round(max(l[,1]),2), 
               "] y[", round(min(l[,2]),2), " ", round(max(l[,2],2)), "]", 
               "|V| =", length(V(g)), "|E| = ", length(E(g)), sep="")
  plot.igraph(g, layout=l, main=text)
}

# g.cur   the current graph which is used for the layout
# members the member after the community detection
# k       the current level which is to be laid to
# l       the layout of the previous level
# iter    the number of iterations (maximum)
# RMSMIN  the root-mean-square, if below this value the algorithm stops earlier
# plot    if true a plot is made
createNewLayoutLevel <- function(g.cur, members, k, l, verb=verb, debug=FALSE, iter=100, RMSMIN=NA, plot=FALSE, method) {
  ############     Remove
  blocksize = iter
  #iter = 3000
  java = FALSE
  if (method == "prefuse") {
    java = TRUE
    if (verb) {
      print("Using prefuse library")
    }
  }
  #if (k > 3) {iter = iter * 3}
  ############     Remove (End)
  
  if (iter >= blocksize) {
    blocks = as.numeric(iter / blocksize);
  } else {
    blocks = 2;
  }
  area <- rep(NA, blocks)
  RMS <- rep(NA, blocks)
  if (verb) {
    print(paste(" Doing layout for level ", k, " iter", iter, " blocks", blocks, " Graph ", sep=""))
    tic <- proc.time()[3]
    print(g.cur)
  }
  com <- NULL
  if (!is.null(l)) { # We have a graph from the previous level, which we use as a starting point
    if (k > 0) {      
      com <- rep(NA, max(members[k,]))
      for (i in 1:length(com)) {
        com[unique(members[k, members[k+1,] == i])] <- i
      }
      #print(com)
    } else {
      com <- rep(NA, length(members[k+1,]))
      for (i in 1:length(com)) {
        com[members[k+1,] == i] <- i     
      }
      #print(com)
    }
    l.init <- matrix(0, nrow=length(com), ncol=2);
    # The scaling 
    nold <- length(unique(members[k + 1,]))
    nnew <- NA
    if (k > 0) {
      nnew <- length(unique(members[k,]))
    } else {
      nnew <- length(V(g.cur))
    }
    fact <- nnew / nold #Scaling, radius assumed to be proportional to number of vetrices this seems to be valid for baramabasi
    for (i in 1:length(com)) { #Setting the starting point
      l.init[i,]  <- fact * l[com[i],];# + rnorm(n=2, mean=0, sd=0.1)
    }
    
    if (debug) {
      if (verb) print(paste("   Using  layout (from level ) ", k , " x[", min(l[,1]), " ", max(l[,1]), "] y[", min(l[,1]), " ", max(l[,2]), "]", sep = ""))
      if (verb) print(paste("   Graph  is weighted?", is.weighted(g), sep=""))
      if (verb) print(paste("   Number of iteration blocks ", blocks, sep=""))
      if (verb) pb <- myTxtProgressBar(min = 0.999, max = blocks, style = 3, file=stderr())
      for (i in 1:blocks) {
          if (verb) setTxtProgressBar(pb, i)
          if (i == 1) {
            start = l.init;
          } else {
            start = l;
          }
          lold <- l;
          if (!java) {
            if (is.weighted(g.cur)) {
                l <- layout.fruchterman.reingold(graph=g.cur, params=list(niter=blocksize, weights=E(g.cur)$weight, start=start, maxdelta=1));
            } else {
                l <- layout.fruchterman.reingold(graph=g.cur, params=list(niter=blocksize, start=start, maxdelta=1));
            }
          } else {
            l <- layout.prefuse(graph=g.cur, pos=start, iter=blocksize, verb=verb)
          }
          area[i] <- max(l[,1]) - min(l[,1]) * max(l[,2]) - min(l[, 2])
          RMS[i]  <- sqrt(sum((l[,1] - start[,1])^2 * (l[,2] - start[,2])^2) / length(l))
          if (plot) {
            plotGraph(g.cur, l, com, RMS[i]);
          }
          if (is.na(RMSMIN) == FALSE && RMS[i] < RMSMIN) {
            if (verb) print(paste("   Aborting early Below minimal RMS ", RMS[i], " after ", i*blocksize, " iterations", sep=""))
            break #Cancel calculation, if RMS is below threshold
          }
      }
      if (verb) close(pb)
    } else { #The all in one go is not implement yet
      print("------------      --------------------    Spring Layout ")
      l <- layout.spring(g.cur, params=list(niter=100,start=l.init))
      #l <- layout.kamada.kawai(graph=g.cur, params=list(niter=500));
    } 
    } else { #No intial Layout given, starting from Scratch
      if (verb) print(paste("  Creating the initial layout (Lowest Level in the Graph) ", is.null(l), sep=""))
      if (verb) pb <- myTxtProgressBar(min = 0.999, max = blocks, style = 3, file=stderr())
      for (i in 1:blocks) {
        if (verb) setTxtProgressBar(pb, i)
        start = l;
        if (!java) {
          if (is.weighted(g.cur)) {
              l <- layout.fruchterman.reingold(g.cur, params=list(start=start, weights=E(g.cur)$weight, niter=blocksize, maxdelta=1))
          } else {
              l <- layout.fruchterman.reingold(g.cur, params=list(start=start, niter=blocksize, maxdelta=1))
          }
        } else {
          l <- layout.prefuse(graph=g.cur, pos=start,iter=blocksize, verb=verb)
        } 
        area[i] <- max(l[,1]) - min(l[,1]) * max(l[,2]) - min(l[, 2])
        if (is.null(start) == FALSE)  RMS[i]  <- sqrt(sum((l[,1] - start[,1])^2 * (l[,2] - start[,2])^2) / length(l))
        if (plot) {
          plotGraph(g.cur, l, com, -1, RMS[i]);
        }
        if (is.na(RMS[i]) == FALSE && is.na(RMSMIN) == FALSE && RMS[i] < RMSMIN) {
          if (verb) print(paste("Below minimal RMS ", RMS[i], " after ", i*blocksize, " iterations", sep=""))
          break #Cancel calculation, if RMS is below threshold
        }
      }
      if (verb) close(pb);
  }
  if (verb) print(paste("Time for layout ", round(proc.time()[3] - tic,5),  sep=""))
  return (list(layout=l, area=area, RMS=RMS));
}

###############
# Own version of progressbar not spoiling the log done
myTxtProgressBar <- function (min = 0, max = 1, initial = 0, char = "=", width = NA, 
                              title, label, style = 1, file = "") 
{
  if (!identical(file, "") && !(inherits(file, "connection") && 
                                  isOpen(file))) 
    stop("'file' must be \"\" or an open connection object")
  if (!style %in% 1L:3L) 
    style <- 1
  .val <- initial
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L
  nw <- nchar(char, "w")
  if (is.na(width)) {
    width <- getOption("width")
    if (style == 3L) 
      width <- width - 10L
    width <- trunc(width/nw)
  }
  if (max <= min) 
    stop("must have 'max' > 'min'")
  up1 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb < nb) {
      cat(paste(rep.int(char, nb - .nb), collapse = ""), 
          file = file)
      flush.console()
    }
    else if (.nb > nb) {
      cat("\r", paste(rep.int(" ", .nb * nw), collapse = ""), 
          "\r", paste(rep.int(char, nb), collapse = ""), 
          sep = "", file = file)
      flush.console()
    }
    .nb <<- nb
  }
  up2 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb <= nb) {
      cat("\r", paste(rep.int(char, nb), collapse = ""), 
          sep = "", file = file)
      flush.console()
    }
    else {
      cat("\r", paste(rep.int(" ", .nb * nw), collapse = ""), 
          "\r", paste(rep.int(char, nb), collapse = ""), 
          sep = "", file = file)
      flush.console()
    }
    .nb <<- nb
  }
  up3 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc) 
      return()
    #     cat(paste(c("\r  |", rep.int(" ", nw * width + 6)), collapse = ""), file = file)
    #     cat(paste(c("\r  |", rep.int(char, nb), rep.int(" ",                       nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""), 
    #         file = file)
    cat(".");
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  getVal <- function() .val
  kill <- function() if (!.killed) {
    cat("\n", file = file)
    flush.console()
    .killed <<- TRUE
  }
  up <- switch(style, up1, up2, up3)
  up(initial)
  structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}


  
