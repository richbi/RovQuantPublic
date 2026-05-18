##----- Wolf OPSCR model
##-- This script loads the data and runs the OPSCR model used in the MINA fagrapport 113 
##-- "Estimates of wolf density, abundance, and population dynamics in Sweden and Norway, 2016–2026".

##-- Clear session and load required libraries
rm(list=ls())
library(nimble)
library(nimbleSCR)


##-- load custom nimble functions and sampler needed to run the model
##-- (located in technicalReports/fagrapport113_Wolves2016-2026/RScripts)
source("RScripts/dbinomLocal_normalWolf.R")


##-- Load Male or Female data input
##-- (located in technicalReports/fagrapport113_Wolves2016-2026/input)

##-- Females
load("input/40.F_2026.1.3_INPUTChain1.RData")

##-- Males
# load("input/40.M_2026.1.3_INPUTChain1.RData")


##-- Fit the OPSCR model with nimbleSCR
##-- Build and compile the OPSCR model
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
cmodel <- compileNimble(model)
cmodel$calculate()


##-- Configure MCMC samplers
MCMCconf <- configureMCMC(model = model,
                          monitors = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1)


##-- Build and run MCMC chains
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
                                                      nburnin = 1,
                                                      niter = 10,
                                                      nchains = 1,
                                                      samplesAsCodaMCMC = TRUE))

