
##' make a dependence plot for both continuous and categorical variables
##' from data frame containing the individualised treatment effects
##' and the given variable
##' 
##' @param variable name of the variable in pd data.frame
##' @param treatment name of the treatment effect variable in pd data.frame
##' @param pd data frame containing the info on the personalised model and the covariates
##' @param nmean logical. TRUE if for a numeric variable the mean should be computed if there
##' are several values of the variable. Defaults to FALSE. Has no impact if variable is a 
##' factor.
##' @param treatment_symbol symbol to be used on the y-axis. If NULL this is "bar(beta)" 
##' or "beta" depending on nmean.
dependenceplot <- function(variable, treatment = "RiluzoleYes", pd, nmean = FALSE, treatment_symbol = NULL) {
  
  betarange <- range(pd[,treatment])
  if(is.null(treatment_symbol)) treatment_symbol <- ifelse(nmean, "bar(beta)", "beta")  
  
  if(is.factor(pd[,variable])) {
    
    p <- ggplot(data = pd, aes_string(x = variable, y = treatment)) +
      geom_boxplot() + ylim(betarange) +
      stat_summary(fun.y = mean, colour = "darkred", geom = "point", 
                   shape = 18, size = 7, show.legend = FALSE) +
      theme(panel.grid.major.x = element_blank()) + 
      ylab(parse(text = paste(treatment_symbol, "(", variable, ")")))
      # ylab(bquote(bar(beta)( .(variable) )))
    
  } else {
    
    if(nmean == TRUE) {
      pd1 <- ddply(pd, as.quoted(variable), 
                   function(x) cbind(treatment = mean(x[,treatment]), nobs = nrow(x)))
    } else {
      pd1 <- pd[ ,c(variable, treatment)]
      names(pd1) <- c(variable, "treatment")
      pd1$nobs <- 0.1
    }
    if(all(pd1$nobs == 0.1)) {
      bk <- waiver()
      rg <- c(1)
    } else {
      maxnobs <- max(pd1$nobs[!is.na(pd1[,variable])])
      bk <- c(1, round(maxnobs/3, 0), round(2*maxnobs/3, 0), maxnobs)
      rg <- c(1, 6)
    }
    
    p <- ggplot(data = pd1, aes_string(x = variable, y = "treatment")) + 
      geom_point(alpha = 0.4, aes(size = nobs)) + 
      geom_smooth(fill = NA) + ylim(betarange) + 
      scale_size_continuous(range = rg, breaks = bk, name = "number of \nobservations") + 
      ylab(parse(text = paste(treatment_symbol, "(", variable, ")")))
      # ylab(bquote(bar(beta)( .(variable) )))
    
    if(nmean == FALSE) p <- p + theme(legend.position = "none") #+ ylab(bquote(beta( .(variable) )))
    
  }
return(p)
}


plot_pda <- function(variable) {
  p <- ggplot(data = pd, aes_string(x = variable, y = "RiluzoleYes"))
  if(is.factor(pd[,variable])) {
    p <- p +
      geom_boxplot() + ylim(betarange) +
      theme(panel.grid.major.x = element_blank())
  } else {
    p <- p + 
      geom_point(alpha = 0.5) + 
      geom_smooth(fill = NA) + ylim(betarange)
  }
  print(p)
}


plot_pdb <- function(variable) {
  pd1 <- ddply(pd, as.quoted(variable), function(x) cbind(RiluzoleYes = mean(x$RiluzoleYes), nobs = nrow(x)))
  
  p <- ggplot(data = pd1, aes_string(x = variable, y = "RiluzoleYes")) + 
    geom_point(alpha = 0.5, aes(size = nobs)) + 
    geom_smooth(fill = NA) + ylim(betarange) 
  
  maxnobs <- max(pd1$nobs[!is.na(pd1[,variable])])
  if(is.factor(pd[,variable])) {
    breaks <- c(min(pd1$nobs), maxnobs)
    p <- p + scale_size_continuous(breaks = breaks, name = "number of \nobservations")
  } else {
    breaks <- c(1, round(maxnobs/3, 0), round(2*maxnobs/3, 0), maxnobs)
    p <- p + scale_size_continuous(breaks = breaks, name = "number of \nobservations")
  }
  print(p)
}
