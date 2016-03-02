## ----data ------------------------------------------------------------

### data preparation ##########################
load("../data/RALSfinal.rda")

ALSFRSdata <- data

hya <- grep(".halfYearAfter", names(ALSFRSdata), value = TRUE)
del <- hya[-grep("ALSFRS", hya)]
st <- grep(".Start", names(ALSFRSdata), value = TRUE)
del <- c(del, st[-grep("ALSFRS", hya)])
del <- c(del, "survival.time", "cens")
ALSFRSdata <- ALSFRSdata[ , !(names(ALSFRSdata) %in% del)]
ALSFRSdata <- ALSFRSdata[complete.cases(ALSFRSdata[, c("ALSFRS.halfYearAfter", 
                                                       "ALSFRS.Start", 
                                                       "Riluzole")]),]


### delete all of non-usefull variables
delete <- grepl("Delta|delta|SubjectID|Unit|Onset|treatment.group", names(ALSFRSdata))
ALSFRSdata <- ALSFRSdata[ , -which(delete)]

### delete Basophil count and Hematocrit (non-explicable values)
ALSFRSdata$Value_Absolute_Basophil_Count <- NULL
# ALSFRSdata$Value_Hematocrit <- NULL
sm1 <- ALSFRSdata$Value_Hematocrit < 1 & !is.na(ALSFRSdata$Value_Hematocrit)
ALSFRSdata$Value_Hematocrit[sm1] <- ALSFRSdata$Value_Hematocrit[sm1] * 100

### delete extremely high Phosphorus value
ALSFRSdata$Value_Phosphorus[!is.na(ALSFRSdata$Value_Phosphorus) & 
                              ALSFRSdata$Value_Phosphorus > 5] <- NA



### delete columns with more than 50% NAs
### except scores, t.onsettrt, Riluzole
keepvarnames <- c("ALSFRS", "t.onsettrt", "Riluzole")
keepvars <- grepl(paste0(keepvarnames, collapse = "|"), names(ALSFRSdata))

pNA <- 0.5 * nrow(ALSFRSdata)
ALSFRSdata <- ALSFRSdata[ , keepvars | (colSums(is.na(ALSFRSdata)) < pNA)]


# ### Rename variables for plotting
names(ALSFRSdata)[names(ALSFRSdata) == "t.onsettrt"] <- "time_onset_treatment"
# names(ALSFRSdata)[names(ALSFRSdata) == "SubjectLiters_fvc"] <- "FVC"
# names(ALSFRSdata)[names(ALSFRSdata) == "Value_Phosphorus"] <- "phosphorus"
Z <- names(ALSFRSdata)[!(names(ALSFRSdata) %in% c("ALSFRS.halfYearAfter", 
                                                  "Riluzole", 
                                                  "ALSFRS.Start"))]
names(ALSFRSdata)[names(ALSFRSdata) %in% Z] <- 
  tolower(names(ALSFRSdata)[names(ALSFRSdata) %in% Z])
names(ALSFRSdata) <- gsub("fam.hist.", "family_history_", names(ALSFRSdata))


dim(ALSFRSdata)
save(ALSFRSdata, file = "../data/ALSFRSdata.rda")



## ----forest ------------------------------------------------------------
message("forest")
set.seed(1234)

### formula
Z <- names(ALSFRSdata)[!(names(ALSFRSdata) %in% c("ALSFRS.halfYearAfter", "Riluzole", "ALSFRS.Start"))]
fm <- as.formula(paste("ALSFRS.halfYearAfter + Riluzole + ALSFRS.Start ~ ", paste(Z, collapse = "+")))

# trace("cforest", quote( message(print(which(sapply(forest, function(x) inherits(x, "partynode")) == FALSE))) ), at = 27, print = TRUE)

### forest
## with cores != NULL, not reproducible (seed depends on parallel)
forest <-  cforest(fm, data = ALSFRSdata, ytrafo = my.lmlog, 
                   ntree = 100, cores = NULL, 
                   perturb = list(replace = FALSE),
                   control = ctree_control(teststat = "max", testtype = "Univ",
                                           mincriterion = 0.95, minsplit = 40, minbucket = 30))
forest <- prune_forest(forest, endpoint = "numeric")

### base model 
bmod <- glm(ALSFRS.halfYearAfter ~ Riluzole + offset(log(ALSFRS.Start)),
            data = ALSFRSdata, family = gaussian(link = "log"))
save(bmod, file = "ALSFRS_bmod.rda")


## ---- pm ----
message("personalised models")
### personalized models
mods <- person_mods(forest, basemod = "glm", newdata = NULL, OOB = TRUE, 
                     offset = "log(ALSFRS.Start)", family = gaussian(link = "log"),
                     parallel = TRUE)


cf <- t(sapply(mods, coef))
summary(cf)
colnames(cf) <- gsub("\\(|\\)|Yes", "", colnames(cf))

#densp <- lapply(colnames(cf), plot_param, dat = cf)

#grid.arrange(densp[[1]] + xlab(bquote(alpha)) + annotate("text", x = coef(bmod)["(Intercept)"], y = 0, label = "alpha", parse = TRUE),
#             densp[[2]] + xlab(bquote(beta)) + annotate("text", x = coef(bmod)["RiluzoleYes"], y = 0, label = "beta", parse = TRUE),
#             nrow = 1, left = "density")


save(cf, file = "ALSFRS_personalModels.rda")


cf <- cbind(cf, ALSFRSdata)
ggplot(cf, aes(x = Intercept, y = Riluzole, color = log(time_onset_treatment))) + geom_point()
ggplot(cf, aes(x = Intercept, y = Riluzole, color = speech)) + geom_point()
ggplot(cf, aes(x = Intercept, y = Riluzole, color = log(value_phosphorus))) + geom_point()
ggplot(cf, aes(x = Intercept, y = Riluzole, color = subjectliters_fvc)) + geom_point()
ggplot(cf, aes(x = Intercept, y = Riluzole, color = age)) + geom_point()
ggplot(cf, aes(x = Intercept, y = Riluzole, color = height)) + geom_point()
ggplot(cf, aes(x = Intercept, y = Riluzole, color = weakness)) + geom_point()


## ----logLiks------------------------------------------------------------
message("logliks")
### logLik
logLiks <- sapply(1:nrow(ALSFRSdata), comp_loglik, mods = mods, dat = ALSFRSdata, 
                  basemod = "glm", loglik = comp_loglik.ALSFRS)

logLik_bmod <- sum(comp_loglik.ALSFRS(mod = bmod, ndat = ALSFRSdata))
(logLik_rf <- sum(logLiks))


save(logLik_bmod, logLik_rf, file = "ALSFRS_logLiks.rda")

rm(mods)



## ----bootstrapLogliks ------------------------------------------------------------
message("bootstrapLogliks")
set.seed(12)

# number of bootstrap samples
B <- 50

# get info to get parametric bootstrap samples
cfbmod <- coef(bmod)
sdbmod <- sqrt(summary(bmod)$dispersion)
mmbmod <- model.matrix(bmod)
y0 <- ALSFRSdata$ALSFRS.Start
yhat <- exp(mmbmod %*% cfbmod) * y0
# all_ynew <- t(sapply(yhat, rnorm, n = B, sd = sdbmod))

### bootstrap sample without negative values
all_ynew <- t(sapply(yhat, function(yh) {
  bss <- rnorm(n = B, mean = yh, sd = sdbmod)
  while(any(bss < 0)) bss[bss < 0] <- rnorm(n = sum(bss < 0), mean = yh, sd =sdbmod)
  return(bss)
}))

# how many observations have a probability of more than 1% to have a negative bootstrap
# sample
table(pnorm(0, mean = yhat, sd = sdbmod) > 0.01)

ggplot(ALSFRSdata, aes(x = Riluzole)) + 
  geom_jitter(aes(y = ALSFRS.halfYearAfter), color = 2, alpha = 0.3, position=position_jitter(height = 0)) +
  geom_jitter(aes(y = all_ynew[,1]), alpha = 0.3, position=position_jitter(height = 0))



get_bslogliks <- function(ynew, start) {
  
  bsdata <- ALSFRSdata
  bsdata$ALSFRS.halfYearAfter <- ynew
  bsforest <-  cforest(fm, data = bsdata, ytrafo = my.lmlog, 
                       ntree = 100, cores = NULL, 
                       perturb = list(replace = FALSE),
                       control = ctree_control(teststat = "max", testtype = "Univ",
                                               mincriterion = 0.95, minsplit = 40, minbucket = 30))
  bsforest <- prune_forest(bsforest, endpoint = "numeric")
  
  bsmods <- person_mods(bsforest, basemod = "glm", newdata = NULL, OOB = TRUE, 
                        offset = "log(ALSFRS.Start)", family = gaussian(link = "log"),
                        parallel = TRUE, start = start)
  bootstrapped.logliks_pm <- sum(sapply(1:nrow(bsdata), comp_loglik, mods = bsmods, dat = bsdata, 
                                           basemod = "glm", loglik = comp_loglik.ALSFRS))
  
  ## compute base model log-likelihood on bootstrap sample
  bsbmod <- glm(ALSFRS.halfYearAfter ~ Riluzole + offset(log(ALSFRS.Start)),
                data = bsdata, family = gaussian(link = "log"), start = cfbmod)
  bootstrapped.logliks_bm <- sum(comp_loglik.ALSFRS(mod = bsbmod, ndat = bsdata))
  
  c(base_model = bootstrapped.logliks_bm,
    forest = bootstrapped.logliks_pm)
}

# bootstrapped.logliks <- as.data.frame(t(apply(all_ynew, 2, get_bslogliks)))

bootstrapped.logliks <- adply(all_ynew, 2, get_bslogliks, start = cfbmod,
                              .progress = "text")

ggplot(bootstrapped.logliks) + 
  geom_line(aes(forest), stat = "density") + 
  geom_line(aes(base_model), stat = "density", linetype = 2) +
  geom_rug(aes(x = logLik_bmod), linetype = 2) +
  annotate("text", x = logLik_bmod, y = 0, 
           label = "base model") +
  geom_rug(aes(x = logLik_rf)) +
  annotate("text", x = logLik_rf, y = 0, 
           label = "forest") +
  xlab("log-likelihood")
  
  

# for(b in 1:B) {
#   message(b)
#   ## get bootstrap sample
#   bssample <- sample(1:nrow(ALSFRSdata), replace = TRUE)
#   bsdata <- ALSFRSdata[bssample, ]
#   
#   
#   ## compute forest log-likelihood on bootstrap sample
#   message("compute forest")
#   bsforest <-  cforest(fm, data = bsdata, ytrafo = my.lmlog, 
#                        ntree = 100, cores = NULL, 
#                        perturb = list(replace = FALSE),
#                        control = ctree_control(teststat = "max", testtype = "Univ",
#                                            mincriterion = 0.95, minsplit = 40, minbucket = 30))
#   bsforest <- prune_forest(bsforest, endpoint = "numeric")
#   
#   message("compute personalized models")
#   bsmods <- person_mods(bsforest, basemod = "glm", newdata = NULL, OOB = TRUE, 
#                         offset = "log(ALSFRS.Start)", family = gaussian(link = "log"),
#                         parallel = TRUE)
#   message("compute logliks")
#   bootstrapped.logliks[b, 2] <- sum(sapply(1:nrow(bsdata), comp_loglik, mods = bsmods, dat = bsdata, 
#                                            basemod = "glm", loglik = comp_loglik.ALSFRS))
#   
#   ## compute base model log-likelihood on bootstrap sample
#   bsbmod <- glm(ALSFRS.halfYearAfter ~ Riluzole + offset(log(ALSFRS.Start)),
#                 data = bsdata, family = gaussian(link = "log"))
#   bootstrapped.logliks[b, 1] <- sum(comp_loglik.ALSFRS(mod = bsbmod, ndat = bsdata))
#   
# }

save(bootstrapped.logliks, file = "ALSFRS_bootstrapLogLiks.rda")


## ----varimp ------------------------------------------------------------
message("varimp")

set.seed(122)
VI <- varimp(forest = forest, basemod = "glm", loglik = comp_loglik.ALSFRS, 
             OOB = TRUE, parallel = TRUE, 
             offset = "log(ALSFRS.Start)", family = gaussian(link = "log"))
VI$variable <- factor(VI$variable, levels = VI$variable[order(VI$VI)])
VI[order(VI$VI), ]

ggplot(VI, aes(y = VI, x = variable)) + geom_bar(stat = "identity", width = .1) + 
  coord_flip() + theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 13))

save(VI, file = "ALSFRS_varimp.rda")



## ----rankplot--------------------------------------------------------------
## get data
V <- as.character(VI$variable[order(VI$VI, decreasing = TRUE)][1:5])
rk <- cbind(cf, ALSFRSdata[, V])
rk$Rank <- rank(rk$Riluzole)

## Plot treatment effect against rank
p_ril <- ggplotGrob(
  ggplot(aes(x = Rank, y = Riluzole), data = rk) + 
    geom_point() + theme_bw() + 
    ylab(bquote(hat(beta)))
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




## ----pdplot ------------------------------------------------------------
library("plyr")
#load(file = "../data/ALSFRSdata.rda")
#load(file = "ALSFRS_forest.rda")
#load(file = "ALSFRS_personalModels.rda")
source("basis/dependence_plots.R")

### partial dependency plots
pd <- cbind(cf, ALSFRSdata)
library("ggplot2")
a <- lapply(Z, dependenceplot, treatment = "Riluzole", pd = pd)
print(a)

b <- lapply(Z, dependenceplot, treatment = "Riluzole", pd = pd, nmean = TRUE)
print(b)
