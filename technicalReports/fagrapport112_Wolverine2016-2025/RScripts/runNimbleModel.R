library(nimble)
library(nimbleSCR)
##----  1. OPSCR models ----
rm(list=ls())
#FOR OPSCR MODELS 
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport112_Wolverine2016-2025/input/OPSCR")
# load the sex-specific input data
# Females
load("54.Cleaned2025Hunn_Chain1.RData")
# Males
#load("54.Cleaned2025Hann_Chain1.RData")

#FOR SCR MODELS 
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport112_Wolverine2016-2025/input/SCR")
# load the sex-specific input data
# Females
# load("Snap2024_54.Cleaned2025Hunn_Chain1.RData")
# Males
#load("Snap2024_54.Cleaned2025Hann_Chain1.RData")

## load a custom nimble function located in fagrapport89_Wolverine2014-2023/RScripts
setwd("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport112_Wolverine2016-2025/RScripts")
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

