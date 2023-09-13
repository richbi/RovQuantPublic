##----- Norwegian Bear OPSCR model
##-- This script loads the data and runs the OPSCR model used in the 
##-- MINA fagrapport 86 - "Estimates of brown bear density, abundance and
##-- population dynamics in Norway 2012 - 2022".

##-- Clear session and load required libraries
rm(list=ls())
library(nimble)
library(nimbleSCR)


##-- load custom nimble functions and sampler needed to run the model
##-- (located in fagrapport86_Bear.NOR2012-2022/RScripts)
source("technicalReports/fagrapport86_Bear_NOR_2012-2022/RScripts/dbinomLocal_normalCovs.R")
source("technicalReports/fagrapport86_Bear_NOR_2012-2022/RScripts/dcatHR.R")
source("technicalReports/fagrapport86_Bear_NOR_2012-2022/RScripts/sampler_categorical_general1.R")


##-- Load Male or Female data input
##-- (located in fagrapport86_Bear.NOR2012-2022/input)

##-- Females
load("input/OPSCR.NOR.F_input.RData")


##-- Males
#load("input/OPSCR.NOR.M_input.RData")


##-- Fit the OPSCR model with nimbleSCR
##-- Build and compile the OPSCR model
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
model$calculate()
cmodel <- compileNimble(model)
cmodel$calculate()


##-- Configure MCMC samplers
MCMCconf <- configureMCMC(model = model,
                          monitors = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1)


##-- Update default sampler with custom sampler for "z" nodes
zIndex <- which(is.na(nimData$z), arr.ind = T)
zIndex <- zIndex[zIndex[ ,2] > 1, ]
zNodes <- paste0("z[",zIndex[,1],",",zIndex[,2],"]")
MCMCconf$removeSamplers(zNodes)
for(i in 1:length(zNodes)){
  MCMCconf$addSampler(target = zNodes[i],
                      type = 'sampler_categorical_general1',
                      control = list("numCategories"= 4))
}#i


##-- Build and run MCMC chains
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
                                                      nburnin = 1,
                                                      niter = 10,
                                                      nchains = 1,
                                                      samplesAsCodaMCMC = TRUE))
