
## -------------------------------------------------------------------------- --
##
## Script name: Spatial variation in individual detectability 
##
## Purpose of script: 
## This R script is meant to demonstrate how to calculate individual detectability
## (defined as the probability to be detected at least once) for any OPSCR model
##  and arbitrary habitat raster.
##
## Author: Pierre Dupont
## Email: pierre.dupont@nmbu.no
##
## Date Created: 2025-01-21
##
## -------------------------------------------------------------------------- --
##
## Notes: For this script to work, you will have to adjust the names and paths 
## to the correct OPSCR model input and processed output. You may also need to
## adjust :
##    1 - The names of the OPSCR model output (e.g. "nimOutput","myNimbleOutput")
##    2 - The dimensions of the parameter sims.list (e.g. "p0", "sigma")
##    3 - 
##
## -------------------------------------------------------------------------- --

rm(list = ls())
gc()


##-- Identify user and set corresponding DropBox and Git/Rovquant directory
if(Sys.info()['user'] == 'pidu') {
  dir.git <- "C:/My_documents/RovQuant"
  dir.dropbox <- "C:/Users/pidu/AQEG Dropbox/AQEG Team Folder/RovQuant"
}
if(Sys.info()['user'] == 'pierredupont') {
  dir.git <- "/Users/pierredupont/Documents/RovQuant"
  dir.dropbox <- "/Users/pierredupont/Dropbox (AQEG)/AQEG Team Folder/RovQuant/"
}
if(Sys.info()['user'] == 'cymi') {
  dir.git <- "C:/My_documents/rovquant/"
  dir.dropbox <- "C:/Users/cymi/Dropbox (Old)/AQEG Dropbox/AQEG Team Folder/RovQuant" 
}
if(Sys.info()['user'] == 'richbi') {
  dir.git <- "C:/Users/richbi/OneDrive - Norwegian University of Life Sciences/PROJECTS/RovQuant"
  dir.dropbox <- "C:/Users/richbi/AQEG Dropbox/AQEG Team Folder/RovQuant"
}
if(Sys.info()['user'] == 'seasunci') {
  dir.git <- "C:/Users/seasunci/02_RovQuant/RovQuant"
  dir.dropbox <- "C:/Users/seasunci/AQEG Dropbox/AQEG Team Folder/RovQuant"
}


##-- Libraries
library(raster)
library(sf)
library(Rcpp)
library(nimble)
library(dplyr)


##-- Source C++ function
sourceCpp(file.path(dir.dropbox,"detection_prob/Individual_detec/FINAL_data_GitHub/cppFunction/GetDetectability.cpp"))

##-- Output directory (to print out results)
dir.out <- file.path(dir.dropbox, "detection_prob/Individual_detec/FINAL_data_GitHub/output/Bear")



#################################################################
############ ---- Example for BEAR - Males ---- #################
#################################################################

## ------ 1. LOAD NECESSARY OBJECTS ------ 

##-- load necessary input (Habitat, Detectors, habitatRasterResolution, country polygon, output OPSCR model)
load(file.path(dir.dropbox, "detection_prob/Individual_detec/FINAL_data_GitHub/input/NecessaryInput_BEAR.RData"))

## ------ 2. INITIAL SET-UP ------

##-- Get detectors coordinates (original)
detectors.xy <- as.matrix(st_coordinates(myDetectors$main.detector.sp))

##-- Get number of years  -> below

##-- Create raster for extraction 
##-- 5km in this case - as for density extraction - but could be something else
regions.r <- habitatRasterResolution$`5km`[["Regions"]]
regions.r <- crop(regions.r, myHabitat$habitat.r)

##-- Set all habitat cells outside the area of interest to NA
regions.r[regions.r[] < 23] <- NA
plot(regions.r)

##-- Identify habitat cells
isHabitat <- which(!is.na(regions.r[]))

##-- Get habitat cell coordinates
regions.xy <- coordinates(regions.r)[isHabitat, ]

##-- Get region name for each habitat cell (format to run the C++ function)
regionsNames <- sort(unique(na.omit(regions.r[])))
regions.rgmx <- do.call(rbind, lapply(regionsNames, function(x)regions.r[] == x))
regions.rgmx[is.na(regions.rgmx)] <- 0
row.names(regions.rgmx) <- factorValues(regions.r, regionsNames)[,1]
regions.rgmx <- regions.rgmx[ ,isHabitat]



## ------ 3. CALCULATE DETECTABILITY -----
## ------   3.1. MALE -----

##-- load input data for males
load(file.path(dir.dropbox, "detection_prob/Individual_detec/FINAL_data_GitHub/input/NOR.2023.7_M_input_1.RData"))

##-- Get number of years  
years <- as.numeric(dimnames(nimData$z)[[2]]) 
n.years <- length(years)

##-- Calculate individual detectability
DetectabilityRegionsM <- list()
for(t in 1:n.years){
  ## Calculate p0 at each detector location
  p0_M <- rep(NA,nimConstants$n.detectors)
  for(j in 1:nimConstants$n.detectors){
    ##--- REPLACE CALCULATION OF p0 DEPENDING ON THE MODEL
    ##--- FOR EXAMPLE THIS WILL HAVE TO BE DUPLICATED FOR THE WOLVERINE
    ##--- ONE CALCULATION FOR SYSTEMATIC AND ONE FOR OPPORTUNISTIC DATA
    p0_M[j] <-
      ilogit(logit(myResults_M$mean$p0[nimConstants$county[j],t]) +    # Intercept, a different one for each "region" (South, Center and North)
               myResults_M$mean$betaDet[1] * nimData$detCovs[j,1,t] +  # Cov1: Distance to roads
               myResults_M$mean$betaDet[2] * nimData$detCovs[j,2,t])   # Cov2: Detection of carnivores (Rovbase & Skandovs) 
  }#j
  
  ## CALCULATE DETECTABILITY USING C++ FUNCTION
  system.time(
    DetectabilityRegionsM[[t]] <- GetDetectability_mean(
      p0 = p0_M,
      sigma = myResults_M$mean$sigma*myHabitat$resolution,
      habitatxy = regions.xy,
      detectorxy = detectors.xy,
      size = nimData$size)         # number of subdetectors per detector
  )
  #print(DetectabilityRegionsM[[t]]$summary)
}#t

## PLOT DETECTABILITY RASTERS FOR EACH YEAR
par(mfrow = c(2,6))
for(t in 1:n.years){
  detectab.r <- regions.r
  detectab.r[isHabitat] <- DetectabilityRegionsM[[t]]$MeanCell
  plot(detectab.r)
}#t



## ------   3.3. SAVE DETECTABILITY OBJECTS ------

save( DetectabilityRegionsM, file = file.path(dir.out, "Detectability_5km_Males.RData"))

# load(file.path(dir.out, "Detectability_5km_Males.RData"))



## ------   3.4. PLOT DETECTABILITY MAPS -----

##-- MALES
pdf(file = file.path(dir.out, "DetectabilityMaps_5km_Males.pdf"),
    width = 20, height = 12)

##-- Set color scale
max <- max(c(unlist(lapply( DetectabilityRegionsM,
                            function(x) max(x$MeanCell[], na.rm = T)))))
cuts <- seq(0, max+0.01*max, length.out = 101) ##-- set breaks
col <- hcl.colors(length(cuts)-1, "YlOrRd", rev = TRUE)

##-- layout
mx <- rbind(c(1,rep(1:6, each = 2)),
            c(rep(1:6, each = 2), 6))
mx <- rbind(mx, mx + 6)
nf <- layout(mx,
             widths = c(rep(1,ncol(mx))),
             heights = rep(1,2))
#layout.show(nf)
par(mar = c(0,0,0,0))

for(t in 1:n.years){
  plot(st_geometry(COUNTRIES), border = NA, col = "gray80")
  detectab.r <- regions.r
  detectab.r[isHabitat] <- DetectabilityRegionsM[[t]]$MeanCell
  image(detectab.r, add=TRUE, breaks = cuts, col = col, legend=FALSE)
  plot(st_geometry(COUNTRIES), border = grey(0.4), col = NA, add=TRUE)
  mtext(text = years[t], side = 1, -8, adj = 0.2, cex = 1.2, font = 2)
  
  if(t == n.years){
    segments(x0 = 900000, x1 = 900000,
             y0 = 6600000, y1 = 6600000 + 500000,
             col = grey(0.3), lwd = 4, lend = 2)
    text(860000, 6600000+500000/2, labels = "500 km", srt = 90, cex = 1.4)
    plot( detectab.r,
          legend.only = T,
          breaks = cuts,
          col = col,
          legend.width = 2,
          axis.args = list( at = round(seq(0,max,length.out = 6),2),
                            labels = round(seq(0,max,length.out = 6),2),
                            cex.axis = 1.2),
          smallplot = c(0.85, 0.9, 0.2, 0.6),
          legend.args = list(text = "Detection prob.",
                             side = 2, font = 1, line = 0.5, cex = 1))
  }#if
}#t
dev.off()

