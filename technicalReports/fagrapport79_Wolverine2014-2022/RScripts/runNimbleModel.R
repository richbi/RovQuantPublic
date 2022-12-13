library(nimble)
library(nimbleSCR)
##----  1. OPSCR models ----
rm(list=ls())
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport79_Wolverine2014-2022/input/OPSCR")

# load the sex-specific input data
# Females
load("47.aJ_FaCleaned1.RData")
# Males
#load("28.M_2021_INPUTChain1.RData")


## load a custom nimble function located in fagrapport74_Wolverine2013-2021/RScripts
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport79_Wolverine2014-2022/RScripts")
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

