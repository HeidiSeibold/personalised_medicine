
##' For tree number b in given forest compute log-likelihood
##' 
##' @param b index of the tree
##' @param forest forest
##' @param formula formula for the model 
##' @param response character string giving the name of the response variable
##' @param basemod character string giving the name of the model used (e.g. "glm")
##' @param loglik a function that computes the log-likelihood when given a model (e.g. comp_loglik.lm)
##' @param OOB should the log-likelihood be computed out-of-bag?
##' @param perm name of the variable which to permute (useful for variable importance)
##' @param ... further details passed on to the model (e.g. family)
##' 
##' Beware!!! 
##' If log-likelihood is computed out-of-bag it can happen that the model is not computable,
##' because there may be no observations of one treatment group in a terminal node for the
##' out-of-bag data.
##' The same problem can occur, if perm != NULL
get_ll <- function(b, forest, formula, response, basemod, loglik, OOB, perm = NULL, ...) {
  
  alldata <- forest$data
  stopifnot(nrow(alldata) > 0)
  alldata$id <- 1:nrow(alldata)
  
  ### get tree and training data
  tree <- forest$nodes[[b]]
  traindat <- alldata[forest$weights[[b]] == 1, ]
  
  ### which data should be used
  if(OOB) {
    samp <- alldata[forest$weights[[b]] == 0, ]
  } else {
    samp <- traindat
  }
  
  if(length(tree) == 0) {
    
    samp$pnode <- 1
    terminals <- 1
    data.terminals <- list(only = traindat)
    
  } else {
    
    ### make tree as constparty object for predicting nodes
    fitted <- data.frame("(fitted)" = fitted_node(tree, data = traindat),
                         "(response)" = traindat[ , response])
    names(fitted) <- c("(fitted)", "(response)")
    ctree <- party(tree, data = traindat, 
                   fitted = fitted)
    
    # cat(paste0(b, ", "))
    # which terminal node is each observation of samp in when variable perm is permuted
    
    
    samp$pnode <- predict(ctree, newdata = samp, type = "node", perm = perm)
    #   samp$node <- predict(ctree, newdata = samp, type = "node")
    
    # ids terminal nodes
    terminals <- nodeids(ctree, terminal = TRUE)
    
    # data in terminal nodes
    data.terminals <- data_party(ctree, terminals)
  }
  
  
  ### compute models in every terminal node
  get.model <- function(data, basemod, formula, ...) {
    if(basemod == "survreg") {
      require("eha")
      m <- aftreg(formula = formula, data = data, ...)    
    } else {
      m <- do.call(basemod, args = list(formula = formula, data = data, ...))
    }
    return(m)
  }
  models.terminals <- lapply(data.terminals, get.model, basemod, formula, ...)
  names(models.terminals) <- terminals
  
  
  
  llfun <- function(x) {
    node <- as.character(unique(x$pnode))
    stopifnot(length(node) == 1)
    sum(loglik(mod = models.terminals[[node]], 
                       ndat = x, response = response))
    
  }
  LL <- sum(daply(samp, .(pnode), .fun = llfun))
  
  return(LL)

}



##' Get forest log-likelihood based on tree log-likelihoods within the forest
##' 
##' @param forest forest
##' @param basemod character string giving the name of the model used (e.g. "glm")
##' @param perm name of the variable which to permute (useful for variable importance)
##' @param OOB should the log-likelihood be computed out-of-bag?
##' @param loglik a function that computes the log-likelihood when given a model (e.g. comp_loglik.lm)
##' @param offset character string that can be included into formula as offset (e.g. "log(ALSFRS.Start)")
##' @param ... further details passed on to the model (e.g. family)
get_loglik <- function(forest, basemod = "lm", perm, OOB = TRUE, loglik, offset, ...) {
    
  terms <- as.character(forest$terms[[2]])[2]
  terms <- strsplit(terms, split = " ")[[1]]
  response <- terms[1]
  

  if(is.null(offset)) {
    
    if(basemod %in% c("survreg", "coxph")) {
      formula <- as.formula(paste("Surv(", terms[1], ",", terms[3], ") ~", paste(terms[seq(5, length(terms), 2)], collapse = "+")))   
      response <- c(terms[1], terms[3])
    } else {
      formula <- as.formula(paste(terms[1], "~", paste(terms[seq(3, length(terms), 2)], collapse = "+")))      
    }
  } else {
    offset <- paste("offset(", offset, ")")
    formula <- as.formula(paste(terms[1], "~", paste(terms[seq(3, length(terms)-1, 2)], collapse = "+"), "+", offset))
  }
#   print(formula)

  logliks <- sapply(1:length(forest$nodes), get_ll, 
                    forest = forest, formula = formula, response = response, 
                    basemod = basemod, loglik = loglik, OOB = OOB, perm = perm, ...)
  
  return(logliks)
}




##' Compute variable importance for a model-based forest
##' 
##' @param forest forest
##' @param basemod character string giving the name of the model used (e.g. "glm")
##' @param OOB should the log-likelihood be computed out-of-bag?
##' @param loglik a function that computes the log-likelihood when given a model (e.g. comp_loglik.lm)
##' @param parallel Should computations be done in parallel using the parallel-package?
##' @param offset character string that can be included into formula as offset (e.g. "log(ALSFRS.Start)")
##' @param ... further details passed on to the model (e.g. family)
varimp <- function(forest, basemod = "lm", loglik = comp_loglik.lm, 
                   OOB = TRUE, parallel = FALSE, offset = NULL, ...) {
  
  zvars <- as.character(forest$terms[[2]][[3]])
  zvars <- gsub("\n    ", "", zvars)
  zvars <- unlist(strsplit(zvars, split = " \\+ ")) 
  zvars <- zvars[zvars %in% attr(forest$terms, "term.labels")]
  print(zvars)
  
  ## get forest log-likelihood
  loglik_base <- get_loglik(forest, loglik = loglik, basemod = basemod, 
                            perm = NULL, OOB = OOB, offset = offset, ...)
  #cat(paste("loglik_base", loglik_base, "\n"))
  
  ## get forest log-likelihood where every time one partitioning variable z is permuted
  if(parallel) {
    require("parallel")
    loglik_perm <- mcmapply(zvars, FUN = function(z) get_loglik(forest = forest,
                                                                loglik = loglik,
                                                                basemod = basemod,
                                                                perm = z,
                                                                OOB = OOB,
                                                                offset = offset, ...),
                            mc.set.seed = FALSE)
  } else {
    loglik_perm <- sapply(zvars, function(z) get_loglik(forest = forest,
                                                        loglik = loglik,
                                                        basemod = basemod,
                                                        perm = z,
                                                        OOB = OOB,
                                                        offset = offset, ...))
  }
  
  ## Variable importances are the means of the tree-wise differences between the original
  ## log-likelihood and the log-likelihood with one z permuted
  diff_logliks <- loglik_base - loglik_perm
  VI <- colMeans(diff_logliks)
  
  return(data.frame(variable = names(VI), VI = VI))
}


