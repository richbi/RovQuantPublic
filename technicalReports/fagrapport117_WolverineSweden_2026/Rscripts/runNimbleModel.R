library(nimble)
library(nimbleSCR)
##----  1. OPSCR models ----
rm(list=ls())
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport117_WolverineSweden_2026/input/SCR")

# load the sex-specific input data
# Females
load("SCRinput2025_2026-10-07_female_1.RData")
# Males
#load("SCRinput2025_2026-10-07_male_1.RData")


## load a custom nimble function located in fagrapport89_Wolverine2014-2023/RScripts
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport117_WolverineSweden_2026/RScripts")
source("dbinomLocal_normalWolverine.R")

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