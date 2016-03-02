##' fitting function for glm with offset and log-link
##' @param data ALSFRS data
##' @param weights weights
##' @param parm which parameters are we interested in. c(1,2) corresponds to intercept and Riluzole parameter.
my.lmlog <- function(data, weights, parm = c(1,2)) {
  
  tb <- table(data[["Riluzole"]][weights > 0])
  ## only one treatment arm left; we don't want to split further...
  # if (any(tb == 0)) return(matrix(0, nrow = nrow(data), ncol = length(parm)))
  if (any(tb < 5)) return(matrix(0, nrow = nrow(data), ncol = length(parm)))
  
  mod <- glm(ALSFRS.halfYearAfter ~ Riluzole + offset(log(ALSFRS.Start)), 
             weights = weights, data = data, subset = weights > 0, 
             family = gaussian(link = "log"), start = c(-0.159, 0.009)) # start from base model
  ef <- as.matrix(estfun(mod)[, parm])
  ret <- matrix(0, nrow = nrow(data), ncol = ncol(ef))
  ret[weights > 0,] <- ef
  ret
}


# ##' compute Log-Likelihoods from personal model (random forest) from a given patient
# ##' @param row row of given patient in data frame
# ##' @param data ALSFRS data
# ##' @param mods list of personal models
# ##' @param bmod base model (parametric model)
# compLogLik <- function(row, data, mods, bmod = NULL) {
#   m <- mods[[row]]
#   d <- data[row, ]
#   y <- d$ALSFRS.halfYearAfter
#   y0 <- d$ALSFRS.Start
#   
#   # rf model
#   yhat <- exp(sum(c(1, d$Riluzole, log(y0)) * c(coef(m), 1)))
#   val <- - log((y - yhat)^2)
#   
#   if(!is.null(bmod)) {
#     # base model
#     byhat <- exp(sum(c(1, d$Riluzole, log(y0)) * c(coef(bmod), 1)))
#     bval <- - log((y - byhat)^2)
#     
#     return(data.frame(logLik = val, base_logLik = bval, resid = yhat - y, base_resid = byhat - y))
#   } else {
#     return(data.frame(logLik = val, resid = yhat - y))
#   }
#   
# }
