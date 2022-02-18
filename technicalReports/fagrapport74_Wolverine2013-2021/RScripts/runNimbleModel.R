library(nimble)
library(nimbleSCR)
# 1. SCR models 
# input file to run the single season SCR models are located here:
# fagrapport74_Wolverine2013-2021/input/singleSeasonSCR
# "Fa" stands for females, and "Ma" for males. One input file is given for each monitoring seasons (9 in total)
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport74_Wolverine2013-2021/input/singleSeasonSCR")
load("25.J_FaSnap2012_1.RData")# load the female input file for the 2012-2013 monitoring season.

## load a custom nimble function located in fagrapport74_Wolverine2013-2021/RScripts
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport74_Wolverine2013-2021/RScripts")
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
#plot check 
overallRuntimeEnd <- proc.time()

# 2. OPSCR models 
rm(list=ls())
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport74_Wolverine2013-2021/input/OPSCR")
load("25.J_Fa1.RData")

## load a custom nimble function located in fagrapport74_Wolverine2013-2021/RScripts
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport74_Wolverine2013-2021/RScripts")
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
#plot check 
overallRuntimeEnd <- proc.time()

