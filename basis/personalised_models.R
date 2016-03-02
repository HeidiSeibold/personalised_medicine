##' compute model
##' 
##' @param w vector of weights
##' @param formula formula for model
##' @param basemod type of model to compute (lm, survreg,...). In case "survreg" is given
##'        actually aftreg form package eha is used 
##' @param dat data 
##' @param y logical. If TRUE response is returned with model fit (see e.g. ?lm)
##' @param returncoefs logical. Should instead of the model fit only the coefficients be
##'        returned?
comp_mod <- function(w, formula, basemod, dat, y, returncoefs = FALSE, ...) {
  if(basemod == "coxph") library("survival")
  
  if(basemod == "survreg") {
    library("eha")
    rrep <- rep(1:nrow(dat), w)
    m <- aftreg(formula = formula, data = dat[rrep, 1:3], y = y)    
  } else {
    m <- do.call(basemod, args = list(formula = formula, data = dat, weights = w, 
                                      y = y, subset = (w > 0), ...))
  }
  
  if(any(is.na(coef(m)))) browser()
  if(returncoefs) {
    return(coef(m))
  } else {
    return(m)    
  }
}


##' personalised models from model-based cforest
##' @param object An object as returned by cforest.
##' @param basemod parametric model used.
##' @param newdata An optional data frame containing test data.
##' @param OOB a logical defining out-of-bag predictions (only if newdata = NULL).
##' @param offset offset set as "variable" or "log(variable)" respectively. 
##'        If not NULL, the last part of the formula stored in object is
##'        assumed to be the offset.
##' @param parallel should the computation be conducted parallel?
##' @param ncores number of cores used for parallel computation
##' @param ... additional arguments passed on to basemod.
person_mods <- function(object, basemod = "lm", newdata = NULL, OOB = TRUE, offset = NULL, 
                        parallel = FALSE, ncores = NULL, returncoefs = FALSE, ...) {
  weights <- predict(object, type = "weights", newdata = newdata, OOB = OOB)
  # message(list(fivenum = fivenum(weights), mean = mean(weights)))
  
  terms <- as.character(object$terms[[2]])[2]
  terms <- strsplit(terms, split = " ")[[1]]
  
  if(!is.character(offset)) {
    if(basemod %in% c("survreg", "coxph")) {
      formula <- as.formula(paste("Surv(", terms[1], ",", terms[3], ") ~", paste(terms[seq(5, length(terms), 2)], collapse = "+")))   
    } else {
      formula <- as.formula(paste(terms[1], "~", paste(terms[seq(3, length(terms), 2)], collapse = "+")))      
    }
  } else {
    offset <- paste("offset(", offset, ")")
    formula <- as.formula(paste(terms[1], "~", paste(terms[seq(3, length(terms)-1, 2)], collapse = "+"), "+", offset))
  }
  
  if(parallel) {
    library("parallel")
    if(is.null(ncores)) ncores <- min(detectCores(), ncol(weights))
    cl <- makeCluster(ncores)
    
    message(paste("Computing on", ncores, "cores."))
    mods <- parApply(cl = cl, weights, 2, comp_mod, formula = formula, 
                     basemod = basemod, dat = object$data, y = TRUE, returncoefs = returncoefs, ...)
    
    stopCluster(cl)
  } else {
    mods <- apply(weights, 2, comp_mod, formula = formula, 
                  basemod = basemod, dat = object$data, y = TRUE, returncoefs = returncoefs, ...)
  }

  return(mods)
}




##' compute log-likelihood contribution(s) from a linear model
##' @param mod object of class lm. 
##' @param ndat data
##' @param response name of the response variable
comp_loglik.lm <- function(mod , ndat, response = "y") {
  y <- ndat[, response]
  yhat <- predict(mod, newdata = ndat)

  - (yhat - y)^2
}


##' compute log-likelihood contribution(s) from a AFT model (aftreg)
##' @param mod object of class lm
##' @param ndat data
comp_loglik.wb <- function(mod, ndat, response = c("survival.time", "cens")) {
  y <- log(ndat[ , response[1]])
  delta <- ndat[ , response[2]]
  
  x <- ndat$Riluzole == "Yes"
  cf <- coef(mod)[1:2]
  cf["RiluzoleYes"] <- - cf["RiluzoleYes"]
  xb <- rowSums(t(rbind(x, 1) * cf))
  
  sc <- exp(-coef(mod)[3])
  
  (- (delta * log(sc)) + (delta * ((y - xb)/sc)) - ( exp((y - xb)/sc) ) )
  
}



##' compute log-likelihood contribution(s) from a gaussian model with log-link
##' @param mod object of class lm
##' @param ndat data
comp_loglik.ALSFRS <- function(mod, ndat, response = "ALSFRS.halfYearAfter") {
  y <- ndat[, response]
  y0 <- ndat$ALSFRS.Start
  
  yhat <- exp(colSums(rbind(1, ndat$Riluzole, log(y0)) * c(coef(mod), 1)))
  val <- - (y - yhat)^2
}




##' compute log-likelihood of i-th observation (i-th personalised model)
##' @param i index of the observation
##' @param mods list of personalised models
##' @param dat data.frame 
##' @param loglik function to compute special loglik. Not needed if basemod = "lm".
comp_loglik <- function(i, mods, dat, basemod = "lm", loglik) {
  ndat <- dat[i, ]
  mod <- mods[[i]]
  if(basemod == "lm" & is.null(loglik)) {
    comp_ll <- comp_loglik.lm
  } else {
    comp_ll <- loglik
  }
  
  comp_ll(mod, ndat = ndat)
}


##' compute log-likelihood contributions from cforest based on personalised models
##' @inheritParams person_mods
##' @param loglik function to compute special loglik. Not needed if basemod = "lm".
loglik_pmods <- function(object, basemod = "lm", newdata = NULL, OOB = TRUE, loglik = NULL, parallel = FALSE, ...) {

  pmods <- person_mods(object, basemod = basemod, newdata = newdata, OOB = OOB, parallel = parallel, ...)
  
  if(is.null(newdata)) newdata <- object$data

  llc <- sapply(1:nrow(newdata), comp_loglik, mods = pmods, dat = newdata, basemod = basemod, loglik = loglik)

  return(llc)
}

##' Plot the distribution of a personalised coefficient
##' @param parameter name of the parameter to plot
##' @param data set with coefficients. Variable names are the parameters.
plot_param <- function(parameter, dat) { 
  dat <- as.data.frame(dat)
  npar <- gsub("\\(|\\)", "", parameter)
  names(dat)[names(dat) == parameter] <- npar
  ggplot(dat, aes_string(x = npar)) + geom_line(stat = "density") + # geom_density() + 
	ylab("")
}
