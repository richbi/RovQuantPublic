library(nimble)
library(nimbleSCR)
##----  1. OPSCR models ----
rm(list=ls())
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport103_Wolves2015-2025/input/OPSCR")

# load the sex-specific input data
# Females
load("39.F_2024_sf_INPUTChain1.RData")
# Males
#load("39.M_2024_sf_INPUTChain1.RData")


## load a custom nimble function located in fagrapport74_Wolverine2013-2023/RScripts
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport97_Wolves2014-2024/RScripts")
source("dbinomLocal_normalWolf.R")

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

