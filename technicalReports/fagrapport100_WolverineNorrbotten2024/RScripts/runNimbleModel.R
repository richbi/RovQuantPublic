library(nimble)
library(nimbleSCR)
##----  1. OPSCR models ----
rm(list=ls())
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport100_WolverineNorrbotten2024/input/SCR")
# load the sex-specific input data
# Females
load("Snap54.FCleaned2024NorrbottenCovWillNewSeason2023_1.RData")
# Males
#load("Snap54.MCleaned2024NorrbottenCovWillNewSeason2023_1.RData")


## load a custom nimble function located in fagrapport89_Wolverine2014-2023/RScripts
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport100_WolverineNorrbotten2024/RScripts")
source("dbin_LESS_Cached_MultipleCovResponse.R")

model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
model$calculate()
cmodel <- compileNimble(model)
cmodel$calculate()
MCMCconf <- configureMCMC(model = model,
                          monitors = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1)

MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
# RUN THE MCMC 
MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
                                                      nburnin = 1,
                                                      niter = 10,
                                                      nchains = 1,
                                                      samplesAsCodaMCMC = TRUE))
overallRuntimeEnd <- proc.time()

