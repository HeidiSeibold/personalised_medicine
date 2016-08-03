
##' For tree b in forest do the pruning
##' 
##' Trick with pruning as long as no pruning neccessary anymore via attr(., "PruningDone")
##' Find a nicer solution???
##' 
##' @param b index of tree
##' @param forest forest
##' @param endpoint character string classifying the endpoint variable (e.g. "survival")
prune_forest.tree <- function(b, forest, endpoint = "survival") {
  
  alldata <- forest$data
  stopifnot(nrow(alldata) > 0)
  
  # get tree and training data
  tree.node <- forest$nodes[[b]]
  traindat <- alldata[forest$weights[[b]] == 1, ]
  tree <- party(tree.node, data = traindat)
  
  # nodeids
  terminals <- nodeids(tree, terminal = TRUE)
  allids <- nodeids(tree)
  if(length(terminals) == 1) return(tree.node) # stop if stump
  
  # data in terminal nodes
  data.terminals <- data_party(tree, terminals)
  
  # which terminal node should be cut off?
  if(endpoint == "survival") {
    needs.pruning <- function(d) {
      tb <- table(d[, c("cens", "Riluzole")])
      (!("1" %in% rownames(tb)) || any(tb["1", ] <  5) || sum(tb["1", ]) < 10) 
    }
  } else {
    needs.pruning <- function(d) {
      tb <- table(d[, "Riluzole"])
      any(tb < 5)
    }
  }
  
  toprune <- sapply(data.terminals, needs.pruning)
  if(sum(toprune) == 0) return(tree.node) # stop if no pruning needed
  terminals.toprune <- terminals[toprune]
  
  # prune at parent of the terminal that should be cut off
  prune_here <- function(id) {
    kids <- sapply(kids_node(nodeapply(tree, id)[[1]]), function(x) x$id)
    any(kids %in% terminals.toprune)
  }
  prune.here <- allids[sapply(allids, prune_here)]
  # nodeprune can not deal with first parent then child 
  # Error: inherits(node, "partynode") is not TRUE
  prune.here <- prune.here[order(prune.here, decreasing = TRUE)] 
 
 
  pruned.tree <- nodeprune(tree.node, prune.here)
  return(pruned.tree)
  
}




##' Prune model-based forest
##' 
##' @param forest forest
##' @param endpoint character string classifying the endpoint variable (e.g. "survival")
prune_forest <- function(forest, endpoint = "survival") {
  
  forest.pruned <- lapply(1:length(forest$nodes), prune_forest.tree, 
                          forest = forest, endpoint = endpoint)
  
  #copied from cforest
  ret <- partykit:::constparties(nodes = forest.pruned, 
                                 data = forest$data, weights = forest$weights,
                                 fitted = forest$fitted, terms = forest$terms, 
                                 info = forest$info)
  class(ret) <- c("cforest", class(ret))
  
  return(ret)
}








