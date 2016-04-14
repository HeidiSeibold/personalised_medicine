
## ---- surv_data --------------------------------------------------------------
load("../data/RALSfinal.rda")  

del <- grepl(".halfYearAfter", names(data))
ALSsurvdata <- data[ , -which(del)]

### delete all rows where there is no survival info
delete <- sapply(ALSsurvdata$survival.time, is.na)
table(delete)
ALSsurvdata <- ALSsurvdata[ -which(delete) , ]


### delete columns with more than 50% NAs
### except scores, t.onsettrt, Riluzole
keepvarnames <- c("survival", "t.onsettrt", "Riluzole")
keepvars <- grepl(paste0(keepvarnames, collapse = "|"), names(ALSsurvdata))


pNA <- 0.5 * nrow(ALSsurvdata)
ALSsurvdata <- ALSsurvdata[ , keepvars | (colSums(is.na(ALSsurvdata)) < pNA)]
names(ALSsurvdata)

keepvars1 <- grepl(paste0(keepvarnames, collapse = "|"), names(ALSsurvdata))
dim(na.omit(ALSsurvdata[ , !keepvars1]))

### delete all of non-usefull variables
delete <- grepl("Delta|delta|SubjectID|Unit|Onset|treatment.group", names(ALSsurvdata))
ALSsurvdata <- ALSsurvdata[ , -which(delete)]

ALSsurvdata <- ALSsurvdata[complete.cases(ALSsurvdata[, c("survival.time", "cens", "Riluzole")]),]
ALSsurvdata$survival.time[ALSsurvdata$survival.time == 0] <- 0.1

# ### Rename variables for plotting
names(ALSsurvdata)[names(ALSsurvdata) == "t.onsettrt"] <- "time_onset_treatment"
Z <- names(ALSsurvdata)[!(names(ALSsurvdata) %in% c("cens", "survival.time", "Riluzole"))]
names(ALSsurvdata)[names(ALSsurvdata) %in% Z] <- tolower(names(ALSsurvdata)[names(ALSsurvdata) %in% Z])
names(ALSsurvdata) <- gsub("fam.hist.", "family_history_", names(ALSsurvdata))


rm(data)
save(ALSsurvdata, file = "../data/ALSsurvdata.rda")




## ---- surv_fittingfunction --------------------------------------------------------------
my.wb <- function(data, weights, parm = c(1, 2, 3)) {
  
  tb <- table(data[weights > 0, c("cens", "Riluzole")])
  ## only one treatment arm left; we don't want to split further...
  if (!("1" %in% rownames(tb)) || any(tb["1", ] <  5) || sum(tb["1", ]) < 10) 
    return(matrix(0, nrow = nrow(data), ncol = length(parm)))
  
  
  mod <- survreg(Surv(survival.time, cens) ~ Riluzole, data = data, 
                 weights = weights, subset = weights > 0, dist = "weibull",
                 init = c(6.7, 0)) # inits based on base model
  ef <- as.matrix(estfun(mod)[,parm])
  ret <- matrix(0, nrow = nrow(data), ncol = ncol(ef))
  ret[weights > 0,] <- ef
  ret
}


## ---- surv_forest --------------------------------------------------------------
set.seed(1235)

### formula
Z <- names(ALSsurvdata)[!(names(ALSsurvdata) %in% c("cens", "survival.time", "Riluzole"))]
fm <- as.formula(paste("survival.time + cens + Riluzole ~ ", paste(Z, collapse = "+")))

bmod_surv <- survreg(Surv(survival.time, cens) ~ Riluzole, data = ALSsurvdata, 
                dist = "weibull")

bmod <- aftreg(Surv(survival.time, cens) ~ Riluzole, data = ALSsurvdata, 
                dist = "weibull")

summary(bmod)
save(bmod, bmod_surv, file = "ALSsurv_bmod.rda")


### forest
## with cores != NULL, not reproducible (seed depends on parallel)
survforest <-  cforest(fm, data = ALSsurvdata, ytrafo = my.wb, ntree = 100, cores = NULL, 
                   perturb = list(replace = FALSE),
                   control = ctree_control(teststat = "max", testtype = "Univ",
                                           mincriterion = 0.95, minsplit = 40, minbucket = 30))

survforest <- prune_forest(survforest)
# save(Z, fm, survforest, file = "ALSsurv_forest.rda")



## ---- surv_pm --------------------------------------------------------------
mods <- person_mods(survforest, basemod = "survreg", newdata = NULL, OOB = TRUE, 
                    parallel = TRUE, init = c(6.7, 0))

cf_surv <- t(sapply(mods, coef))
summary(cf_surv)
colnames(cf_surv) <- gsub("Yes", "", colnames(cf_surv))

save(cf_surv, file = "ALSsurv_personalModels.rda")


cf_surv[,"Riluzole"] <- - cf_surv[,"Riluzole"]
cf_surv[,"log(shape)"] <- - cf_surv[,"log(shape)"]



cf_surv <- cbind(cf_surv, ALSsurvdata)
names(cf_surv) <- gsub("\\(|\\)", "", names(cf_surv))

ggplot(data = cf_surv, aes(x = logscale, y = logshape, color = Riluzole)) + 
  geom_point(size = I(1.2)) + 
  labs(x = bquote(log(alpha[2])), y = bquote(alpha[1]), color = bquote(beta)) +
  scale_colour_gradient(limits=c(min(cf_surv$Riluzole), max(cf_surv$Riluzole)), 
                        low="green", high="black", space="Lab")



## ----surv_varimp--------------------------------------------------------------
set.seed(444)
VI <- varimp(forest = survforest, basemod = "survreg", loglik = comp_loglik.wb, 
             OOB = TRUE, parallel = TRUE)
VI$variable <- factor(VI$variable, levels = VI$variable[order(VI$VI)])
VI[order(VI$VI), ]


save(VI, file = "ALSsurv_varimp.rda")


ggplot(VI, aes(y = VI, x = variable)) + 
  geom_bar(stat = "identity", width = .1) + 
  coord_flip() +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 13))

rm(survforest)


## ----surv_pdplot--------------------------------------------------------------
### partial dependency plots
load("ALSsurv_personalModels.rda")
cf_surv <- as.data.frame(cf_surv)

smed <- ldply(1:nrow(cf_surv), function(i) {
  
  x <- cf_surv[i,]
  
  p = 0.5
  s0 <- qweibull(p = p, shape = 1/exp(-x[,"log(shape)"]), 
                 scale = exp(x[,"log(scale)"]))
  s1 <- qweibull(p = p, shape = 1/exp(-x[,"log(shape)"]), 
                 scale = exp(x[,"log(scale)"] - x[,"Riluzole"]))
  
  data.frame(id = i, x, s0, s1, sdiff = s1 - s0,
             row.names = NULL)
  
})


pd <- cbind(smed, ALSsurvdata)

a <- lapply(Z, dependenceplot, treatment = "sdiff", pd = pd, treatment_symbol = "Delta[.5]")
print(a)

b <- lapply(Z, dependenceplot, treatment = "Riluzole", pd = pd, nmean = TRUE, treatment_symbol = "bar(Delta)[.5]")
print(b)

rm(pd)
rm(a)
rm(b)

## ----surv_logLiks-------------------------------------------------------------
set.seed(5)

## forest
logLiks <- sapply(1:nrow(ALSsurvdata), comp_loglik, mods = mods, dat = ALSsurvdata, 
                  basemod = "survreg", loglik = comp_loglik.wb)

## base model
logLik_bmod_surv <- sum(comp_loglik.wb(mod = bmod, ndat = ALSsurvdata))
(logLik_rf_surv <- sum(logLiks))


## forest with splits in alpha_1 and alpha_2
my.wb_alpha <- function(data, weights, parm = c(1, 3)) {
  my.wb(data, weights, parm)
}
survforest_alpha <-  cforest(fm, data = ALSsurvdata, ytrafo = my.wb_alpha, ntree = 100, cores = NULL, 
                             perturb = list(replace = FALSE),
                             control = ctree_control(teststat = "max", testtype = "Univ",
                                                     mincriterion = 0.95, minsplit = 40, minbucket = 30))
survforest_alpha <- prune_forest(survforest_alpha)
mods_alpha <- person_mods(survforest_alpha, basemod = "survreg", newdata = NULL, OOB = TRUE, 
                          parallel = TRUE, init = c(6.7, 0))

logLiks_alpha <- sapply(1:nrow(ALSsurvdata), comp_loglik, mods = mods_alpha, dat = ALSsurvdata, 
                        basemod = "survreg", loglik = comp_loglik.wb)
(logLik_rf_surv_alpha <- sum(logLiks_alpha))

## forest with splits in beta
my.wb_beta <- function(data, weights, parm = c(2)) {
  my.wb(data, weights, parm)
}
survforest_beta <-  cforest(fm, data = ALSsurvdata, ytrafo = my.wb_beta, ntree = 100, cores = NULL, 
                             perturb = list(replace = FALSE),
                             control = ctree_control(teststat = "max", testtype = "Univ",
                                                     mincriterion = 0.95, minsplit = 40, minbucket = 30))
survforest_beta <- prune_forest(survforest_beta)
mods_beta <- person_mods(survforest_beta, basemod = "survreg", newdata = NULL, OOB = TRUE, 
                          parallel = TRUE, init = c(6.7, 0))

logLiks_beta <- sapply(1:nrow(ALSsurvdata), comp_loglik, mods = mods_beta, dat = ALSsurvdata, 
                        basemod = "survreg", loglik = comp_loglik.wb)
(logLik_rf_surv_beta <- sum(logLiks_beta))

save(logLik_bmod_surv, logLik_rf_surv, logLik_rf_surv_alpha, logLik_rf_surv_beta,
     file = "ALSsurv_logLiks.rda")


rm(mods)
rm(bmod)
rm(mods_alpha)
rm(mods_beta)



## ----surv_bootstrapLogliks ------------------------------------------------------------
set.seed(12)

# number of bootstrap samples
B <- 50

# get info to get parametric bootstrap samples
mmbmod <- model.matrix(bmod_surv)
rweib_scale <- exp(mmbmod %*% coef(bmod_surv))
rweib_shape <- 1/bmod_surv$scale
maxt <- max(ALSsurvdata$survival.time)


## i: patient i
## b: bootstrapsample b
all_ynew <- ldply(1:nrow(ALSsurvdata), function(i) {
  
  sc <- rweib_scale[i]
  stime <- rweibull(n = B, shape = rweib_shape, scale = sc)
  
  if(ALSsurvdata$cens[i] == 0) {
    censored <- stime > ALSsurvdata$survival.time[i]
    stime[censored] <- ALSsurvdata$survival.time[i]
  } else {
    censored <- stime > maxt
    stime[censored] <- maxt
  }
  
  cens <- as.numeric(!censored)
  stime < maxt
  
  data.frame(i = rep(i, B), b = 1:B, survival.time = stime, cens = cens)
})

ggplot(ALSsurvdata, aes(survival.time, color = factor(cens))) + 
  geom_line(stat = "density") + geom_rug()

fit <- survfit(Surv(survival.time, cens) ~ 1, data = ALSsurvdata) 
plot(fit) 

# check censoring percentage
n <- nrow(ALSsurvdata)
pc <- round(sum(ALSsurvdata$cens == 0) / n * 100, digits = 2)
pc_bs <- round(daply(all_ynew, .(b), function(yn) sum(yn$cens == 0) / n * 100), digits = 2)
message(
  paste("Original data has",
        pc, 
        "% censoring. \nSimulated has", 
        paste(pc_bs, collapse = "; "))
)


get_bslogliks_surv <- function(ynew) {
  
  bsdata <- ALSsurvdata
  bsdata$survival.time <- ynew$survival.time
  bsdata$cens <- ynew$cens
  
  ## compute forest log-likelihood on bootstrap sample
  bssurvforest <-  cforest(fm, data = bsdata, ytrafo = my.wb, ntree = 100, cores = NULL, 
                           perturb = list(replace = FALSE),
                           control = ctree_control(teststat = "max", testtype = "Univ",
                                                   mincriterion = 0.95, minsplit = 40, minbucket = 30))
  bssurvforest <- prune_forest(bssurvforest) 
  
  bsmods <- person_mods(bssurvforest, basemod = "survreg", newdata = NULL, OOB = TRUE, 
                        parallel = TRUE, ncores = 75, init = c(6.7, 0))
  rm(bssurvforest)
  bootstrapped.logliks_pm <- sum(sapply(1:nrow(bsdata), comp_loglik, mods = bsmods, 
                                                dat = bsdata, basemod = "survreg", 
                                                loglik = comp_loglik.wb))
  
  ## compute base model log-likelihood on bootstrap sample
  bsbmod <- aftreg(Surv(survival.time, cens) ~ Riluzole, data = bsdata, 
                   dist = "weibull")
  bootstrapped.logliks_bm <- sum(sum(comp_loglik.wb(mod = bsbmod, ndat = bsdata)))
  
  c(base_model = bootstrapped.logliks_bm,
    forest = bootstrapped.logliks_pm)
}


bootstrapped.logliks_surv <- ddply(all_ynew, .(b), get_bslogliks_surv, .progress = "text")
message("BS done")

save(bootstrapped.logliks_surv, file = "ALSsurv_bootstrapLogLiks.rda")

rm(bootstrapped.logliks_surv)






## ----surv_rankplot--------------------------------------------------------------
### Median Survival improvement ###

load("ALSsurv_personalModels.rda")
cf_surv <- as.data.frame(cf_surv)

smed <- ldply(1:nrow(cf_surv), function(i) {
  
  x <- cf_surv[i,]
  
  p = 0.5
  s0 <- qweibull(p = p, shape = 1/exp(-x[,"log(shape)"]), 
                 scale = exp(x[,"log(scale)"]))
  s1 <- qweibull(p = p, shape = 1/exp(-x[,"log(shape)"]), 
                 scale = exp(x[,"log(scale)"] - x[,"Riluzole"]))
  
  data.frame(id = i, x, s0, s1, sdiff = s1 - s0,
             row.names = NULL)
  
})

## get data
V <- as.character(VI$variable[order(VI$VI, decreasing = TRUE)][1:5])
rk <- cbind(smed, ALSsurvdata[, V])
rk$Rank <- rank(rk$sdiff)

## Plot median survival against rank
p_ril <- ggplotGrob(
  ggplot(aes(x = Rank, y = sdiff), data = rk) + 
    geom_point() + theme_bw() + 
    ylab("Difference in \n median survival")
)

## Plot variables with highest VI against rank
make_ggplotgrob <- function(z) {
  
  p <- ggplot(aes_string(x = "Rank", y = z), data = rk) +
    theme_bw() + theme(axis.ticks.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       legend.position = "none")
  
  if(is.numeric(rk[,z])) {
    p <- p + geom_point(alpha = I(0.2))
  } else {
    p <- p + geom_point(alpha = I(0.2), aes_string(colour = z))
  }
  
  ggplotGrob(p)
  
}

p_z <- lapply(V, make_ggplotgrob)
names(p_z) <- V


## arrange plots in list, make sure they align
p_all <- list()
p_all[[1]] <- p_ril
p_all[2:(length(p_z) + 1)] <- p_z


maxWidth <- do.call(grid:::unit.pmax, lapply(p_all, function(p) p$widths[2:3]))

p_all <- lapply(p_all, function(p) {
  p$widths[2:3] <- maxWidth
  p
})

## Plot
do.call(grid.arrange, c(p_all, ncol = 1))







## ---- surv_newmodel --------------------------------------------------------------
ALSsurvdata$totr <- 1 
ALSsurvdata$totr[ALSsurvdata$time_onset_treatment < 600] <- 0
ALSsurvdata$totr <- factor(ALSsurvdata$totr)

nbmod <- coxph(Surv(survival.time, cens) ~ 
                Riluzole * age + 
                Riluzole * time_onset_treatment : totr, 
              data = ALSsurvdata)
summary(nbmod)
