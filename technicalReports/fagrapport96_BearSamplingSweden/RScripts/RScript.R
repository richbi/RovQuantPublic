rm(list=ls())
#LOAD LIBRARIES
library(nimble)
library(nimbleSCR)
library(raster)
library(sf)

#LOAD REQUIRED DATA AND FUNCTIONS
load("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport94_BearSamplingSweden/input/simulationSCR/RData.RData")
source("C:/My_documents/rovquant/analyses/Rgit/RovQuantPublic/technicalReports/fagrapport94_BearSamplingSweden/RScripts/dbinomLocal_normalCovs.R")

## IN THIS SIMULATION SET UP, DENSITY AND DETECTIONS ARE PREDICTED FOR THE TWO LAST YEARS (2020-2021) USING RESULTS FROM THE OPSCR MODEL (Dupont et al 2024)
## ---- 1. SIMULATE DENSITY  ----
## ----   1.1 MODEL CODE  ----
modelCodeDensity <- nimbleCode(
  {
    habIntensity[1:n.habwindows] <- exp(betaDens[1] * habDens[1:n.habwindows,1] +
                                          betaDens[2] * habDens[1:n.habwindows, 2])
    sumHabIntensity <- sum(habIntensity[1:n.habwindows])
    logHabIntensity[1:n.habwindows] <- log(habIntensity[1:n.habwindows])
    logSumHabIntensity <- log(sumHabIntensity)
    for (t in 1:n.years) {
      for (i in 1:n.individuals) {
        sxy[i, 1:2,t] ~ dbernppAC(lowerCoords = lowerHabCoords[1:n.habwindows,1:2],
                                  upperCoords = upperHabCoords[1:n.habwindows,1:2],
                                  logIntensities = logHabIntensity[1:n.habwindows],
                                  logSumIntensity = logSumHabIntensity,
                                  habitatGrid = habitatGrid[1:y.max,1:x.max],
                                  numGridRows = y.max,
                                  numGridCols = x.max)
        
      }
    }
  }
)

## ----   1.2 BUILD MODEL  ----
modelsxy <- nimbleModel( code = modelCodeDensity,
                         constants = nimConstantsDens,
                         data = nimDataDens,
                         check = F,       
                         calculate = F)  

## -----  1.3 SIMULATE ------
nodesToSimsxy <- modelsxy$getDependencies(c("sxy","habIntensity","sumHabIntensity",
                                            "logHabIntensity","logSumHabIntensity"),
                                          self = T,
                                          downstream = F,
                                          returnScalarComponents = TRUE)

#USE POSTERIOR DISTRIBUTION FROM OPSCR (DUPONT ET AL 2024)
i=1 #we chose the first posterior
modelsxy$betaDens <- betaDensPosterior[i,]


#SIMULATE
set.seed(i)
modelsxy$simulate(nodes = nodesToSimsxy,
                  includeData = F)

## -----  1.4 CHECK AND PLOT PREDICTED DENSITY ------
t=2 #we plot the last year
densPred <- calculateDensity(modelsxy$sxy[,,t],
                             habitatGrid = nimData$habitatGrid,
                             numWindows = nimConstants$n.habwindows,
                             nIndividuals = dim(zPosterior)[2],
                             indicator = zPosterior[i,,t])
#FILL IN PREDICTED DENSITY IN RASTERS
r <- habitat.r
habitat.r[habitat.r[]%in% 1] <- densPred
#PLOT
plot(habitat.r)

## ---- 2. SIMULATE DETECTIONS  ----
## ----   2.1 MODEL CODE  ----
modelCode <- nimbleCode(
  {
    sigma ~ dunif(0, 5)
    
    for (d in 1:n.detCovs) {
      betaDet[d] ~ dunif(-5, 5)
    }
    
    
    for (c in 1:n.counties) {
      p0[c] ~ dunif(0, 1)
    }
    
    for (i in 1:n.individuals){
      for (t in 1:n.years) {
        y.alive[i, 1:lengthYCombined, t] ~ dbinomLocal_normalCovs(s = sxy[i,1:2, t],
                                                                  size = size[1:n.detectors], 
                                                                  p0Traps = p0[1:n.counties],
                                                                  sigma = sigma, 
                                                                  trapCoords = detCoords[1:n.detectors, 1:2], 
                                                                  localTrapsIndices = localDetIndices[1:n.habwindowsLocal,1:localDetNumMax],
                                                                  localTrapsNum = localDetNum[1:n.habwindowsLocal],
                                                                  resizeFactor = resizeFactor,
                                                                  habitatGrid = habitatGridLocal[1:y.maxLocal,1:x.maxLocal], 
                                                                  indicator = (z[i, t]), 
                                                                  lengthYCombined = lengthYCombined, 
                                                                  allowNoLocal = 1, 
                                                                  trapCovsIntercept = county[1:n.detectors],
                                                                  trapCovs = detCovs[1:n.detectors,1:n.detCovs], 
                                                                  trapBetas = betaDet[1:n.detCovs])
      }
    }
  }
)

## ----   2.2 BUILD NIMBLE MODEL  ----
# nimData$size <- nimData$size[,1]

# plot(nimData$detCoords[,2,2]~nimData$detCoords[,1,2])

# BUILD NIMBLE MODEL 
modelSims <- nimbleModel( code = modelCode,
                          constants = nimConstants,
                          data = nimData,
                          check = F,       
                          calculate = F)  

## ----   2.3 SIMULATE DETECTION DATA ----
nodesToSim <- modelSims$getDependencies(c("z"),
                                        self = T,
                                        downstream = F,
                                        returnScalarComponents = TRUE)

#SET P0 VALUES TO GET THE DIFFERENT INTENSITIES 
# SCENARIO 1
modelSims$p0 <- c(0.0027,0.0027,0.0046,0.0075,0.0050,0.0065)
# #SCENARIO 2
# model$p0 <- c(0.0027,0.0027,0.0046,0.0075,0.0050,0.0065)/2
# #SCENARIO 3
# model$p0 <- c(0.0027,0.0027,0.0046,0.0075,0.0050,0.0065)/4
# #SCENARIO 4
# model$p0 <- c(0.0027,0.0027,0.0046,0.0075,0.0050,0.0065)/8
#USE THE PREDICTED SXY OBTAINED IN  1.
modelSims$sxy <- modelsxy$sxy #Scaledsxy[itera[ite],,,(n.yearsTotal-n.years+1):n.yearsTotal] 
#Z 
modelSims$z <- zPosterior[i,,]
#sigma
modelSims$sigma <- sigmaPosterior[i]


#BETA DETS
modelSims$betaDet[1] <- 1 #EFFECT OF THE OBS COVARIATE (CONTINUOUS)
modelSims$betaDet[2] <- 0 #EMPTY COLUMN FOR THE EFFECT OF CATEGORIES (FOR FITTING)
modelSims$betaDet[3] <- betaDetPosterior[i,1]#roads  

#SIMULATE
modelSims$simulate(nodes = nodesToSim,
                   includeData = F)


## ---- 3. FIT MODEL STRUCTURED SAMPLING ----
## ----   3.1 MODEL CODE ----
modelCodeFit <- nimbleCode({
  
  betaDens[1] ~ dnorm(0, 0.01)
  betaDens[2] ~ dnorm(0, 0.01)
  habIntensity[1:n.habwindows] <- exp(betaDens[1] * habDens[1:n.habwindows,1] +
                                        betaDens[2] * habDens[1:n.habwindows, 2])
  sumHabIntensity <- sum(habIntensity[1:n.habwindows])
  logHabIntensity[1:n.habwindows] <- log(habIntensity[1:n.habwindows])
  logSumHabIntensity <- log(sumHabIntensity)
  for (i in 1:n.individuals) {
    sxy[i, 1:2] ~ dbernppAC(lowerCoords = lowerHabCoords[1:n.habwindows,1:2],
                            upperCoords = upperHabCoords[1:n.habwindows,1:2],
                            logIntensities = logHabIntensity[1:n.habwindows],
                            logSumIntensity = logSumHabIntensity,
                            habitatGrid = habitatGrid[1:y.max,1:x.max],
                            numGridRows = y.max,
                            numGridCols = x.max)
    
  }
  
  psi ~ dunif(0,1)
  for (i in 1:n.individuals) {
    z[i] ~ dbern(psi)
  }
  
  sigma ~ dunif(0, 5)
  for (d in 1:n.detCovs) {
    betaDet[d] ~ dunif(-5, 5)
  }
  # for (t in 1:n.years) {
  for (c in 1:n.counties) {
    p0[c] ~ dunif(0, 1)
  }
  for (i in 1:n.individuals) {
    y.alive[i, 1:lengthYCombined] ~ dbinomLocal_normalCovs(s = sxy[i,1:2],
                                                           size = size[1:n.detectors], 
                                                           p0Traps = p0[1:n.counties],
                                                           sigma = sigma, 
                                                           trapCoords = detCoords[1:n.detectors, 1:2], 
                                                           localTrapsIndices = localDetIndices[1:n.habwindowsLocal,1:localDetNumMax],
                                                           localTrapsNum = localDetNum[1:n.habwindowsLocal],
                                                           resizeFactor = resizeFactor,
                                                           habitatGrid = habitatGridLocal[1:y.maxLocal,1:x.maxLocal], 
                                                           indicator = (z[i] == 1), 
                                                           lengthYCombined = lengthYCombined, 
                                                           allowNoLocal = 1, 
                                                           trapCovsIntercept = county[1:n.detectors],
                                                           trapCovs = detCovs[1:n.detectors,1:n.detCovs], 
                                                           trapBetas = betaDet[1:n.detCovs])
  }
  N <- sum(z[1:n.individuals])
  
})

## ----   3.2 SET UP THE SIMUALTED DATA FOR FITTING ----
#We simulate two years but fit only the last one. 
t=2
ySparse <- modelSims$y.alive
#reduce the size of the dataset 
whichDets <- which(ySparse[,1,t]>0)
whichNOTDets <- which(ySparse[,1,t]%in%0)

#AUGMENT DATA 
#set M to 4000
M <- 4000
#sample a X number of ids 
whichNOTDets  <-  sample(whichNOTDets, size = M - length(whichDets), replace = T)
whichidToKeep <- c(whichDets, whichNOTDets)#

#NIMDATA 
nimData1 <- nimData
nimData1$y.alive <- ySparse[whichidToKeep,,t]
nimData1$lowerHabCoords <- ScaledLowUpCoords$lowerHabCoords
nimData1$upperHabCoords <- ScaledLowUpCoords$upperHabCoords

#PARAMETERS TO BE ESTIMATED
nimData1$z <- NULL
nimData1$sxy <- NULL
nimData1$p0 <- NULL
nimData1$betaDet <- NULL
nimData1$sigma <- NULL
nimData1$betaDens <- NULL


#NIMCONSTANTS 
nimConstants1 <- nimConstants
nimConstants1$n.individuals <- length(whichidToKeep)
nimConstants1$n.habwindows <- dim(nimData1$habDens)[1]
nimConstants1$habitatGrid <- ScaledLowUpCoords$habitatGrid
nimConstants1$x.max <- dim(ScaledLowUpCoords$habitatGrid)[2]
nimConstants1$y.max <- dim(ScaledLowUpCoords$habitatGrid)[1]
#nimConstants1$lengthYCombined <- ySparse$lengthYCombined
nimConstants1$period <- NULL
nimConstants1$localHabNumMax <- NULL
nimConstants1$n.years <- NULL
nimConstants1$n.regions <- NULL


# INITIAL VALUES 
initsz <- modelSims$z[whichidToKeep,t]

nimInits1 <- list(betaDens = runif(2,-0.1,0.1),
                  sigma = runif(1,0.5,1),
                  p0 = runif(nimConstants$n.counties,0.1,0.2),
                  betaDet =runif(nimConstants1$n.detCovs,-0.1,0.1),
                  psi = runif(1,0.4,0.5),
                  z = initsz,
                  sxy = modelSims$sxy[whichidToKeep,,t])

#PARAMETERS TO MONITOR 
nimParams <- c("betaDens","sigma","p0","betaDet","psi","N")
nimParams2 <- c("z","sxy")

#SAVE TRUE PARAMETER VALUES 

## ----   3.3 FITTING THE MODEL ----
modelStructured <- nimbleModel( code = modelCodeFit,
                                constants = nimConstants1,
                                data = nimData1,
                                inits = nimInits1,
                                check = F,       
                                calculate = F)  
modelStructured$calculate()
cmodelStructured <- compileNimble(modelStructured)
MCMCconfStructured <- configureMCMC(model = modelStructured,
                                    monitors  = c(names(nimInits1),"N"),
                                    control = list(reflective = TRUE),
                                    thin = 1)
MCMCStructured <- buildMCMC(MCMCconfStructured)
cMCMCStructured <- compileNimble(MCMCStructured, project = modelStructured, resetFunctions = TRUE)
## Run MCMC (RUN FOR LONGER)
MCMCRuntimeStructured <- system.time(samplesStructured <- runMCMC( mcmc = cMCMCStructured,
                                                                   nburnin = 100,
                                                                   niter = 500,
                                                                   nchains = 1,
                                                                   samplesAsCodaMCMC = TRUE))

## ---- 4. FIT MODEL OPPORTUNISTIC SAMPLING ----
## ----   4.1 CREATE PROXY IN SPATIAL EFFORT ----
# WE INTRODUCE SOME ERRORS IN THE COVARIATES OF SPATIAL EFFORT TO SIMULATE WHAT COULD OCCUR IF A PROXY WAS USED 
eps <- 0.15
error.r[] <- rnorm(length(error.r), 0, eps)
#------ ADD AUTOCORRELATION TO ERROR
error.r<-
  focal(
    error.r,
    w = matrix(1, nrow = 5, ncol = 5),
    pad = TRUE,
    padValue = 0
  )
#------ ADD AUTOCORRELATED ERROR TO DET PROB RASTER
with.error.r <- orig.r + error.r

#------ SCALE AGAIN
with.error.r <- scale(with.error.r)

#------  CATEGORIZE THE p0 COVARIATE
splits.orig <- quantile(orig.r, p=c(0.45, 0.9), na.rm=TRUE)
with.error.categ.r <- orig.categ.r <- orig.r
with.error.categ.r[] <- cut(with.error.r[], 
                            breaks = c(min(with.error.r[], na.rm=TRUE)*1.1,
                                       splits.orig, max(with.error.r[], na.rm=TRUE))
)

orig.categ.r[]<- cut(orig.r[], 
                     breaks = c(min(with.error.r[], na.rm=TRUE)*1.1,
                                splits.orig, max(orig.r[], na.rm=TRUE))
)


#CATEGORICAL USED FOR FITTING
detCat <- raster::extract(with.error.categ.r,detectors)
#deal with NAs
isna <- which(is.na(detCat))
tmp <- raster::extract(with.error.categ.r,detectors[isna, ],
                       buffer = 15000, fun = min, na.rm = T)
detCat[isna] <- tmp
##turn to category
detCovsCat <- matrix(0,nrow=length(detCat),ncol=2)
detCovsCat[,1] <- ifelse(detCat %in% 2, 1, 0)
detCovsCat[,2] <- ifelse(detCat %in% 3, 1, 0)
#quick check
colSums(detCovsCat)

#DETECTOR COVARIATE 
#CATEGORICAL COV + ROAD 
nimData1$detCovs <- cbind(detCovsCat, nimData$detCovs[,3]) #cbind(nimData1$detCovs[,1],detOtherTrueBinary )
nimConstants1$n.detCovs <- dim(nimData1$detCovs)[2]

## ----   4.2 FIT MODEL----
modelOpportunistic <- nimbleModel( code = modelCodeFit,
                                   constants = nimConstants1,
                                   data = nimData1,
                                   inits = nimInits1,
                                   check = F,       
                                   calculate = F)  
modelOpportunistic$calculate()
cmodelOpportunistic <- compileNimble(modelOpportunistic)
MCMCconfOpportunistic <- configureMCMC(model = modelOpportunistic,
                                       monitors  = c(names(nimInits1),"N"),
                                       control = list(reflective = TRUE),
                                       thin = 1)
MCMCOpportunistic <- buildMCMC(MCMCconfOpportunistic)

cMCMCOpportunistic <- compileNimble(MCMCOpportunistic, project = modelOpportunistic, resetFunctions = TRUE)
gc()
## Run MCMC (RUN FOR LONGER)
MCMCRuntimeOpportunistic <- system.time(samplesOpportunistic <- runMCMC( mcmc = cMCMCOpportunistic,
                                                                         nburnin = 100,
                                                                         niter = 500,
                                                                         nchains = 1,
                                                                         samplesAsCodaMCMC = TRUE))


##PLOT RESULTS 
library(basicMCMCplots)
chainsPlot(samplesOpportunistic,var = c("N","betaDens","sigma","p0","psi"))#,line = c(mean(model$sigma),nimData$p0,sum(nimData$z[,t])/nimConstants1$n.individuals))
chainsPlot(samplesStructured,var = c("N","betaDens","sigma","p0","psi"))#,line = c(mean(model$sigma),nimData$p0,sum(nimData$z[,t])/nimConstants1$n.individuals))


