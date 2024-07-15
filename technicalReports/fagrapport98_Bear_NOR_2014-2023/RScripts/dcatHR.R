
#### 1.Density function ####
dcatHR <- nimbleFunction(run = function( 
    x = double(0),
    z = double(0),
    gamma = double(0),
    mhH = double(0),
    mhW = double(0),
    log = integer(0, default = 0)
){
  # Return type declaration
  returnType(double(0))
  
  if(z == 1){
    logLikelihood <- dcat(x, prob = c(1 - gamma, gamma), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 2){
    mhH1 <- exp(mhH)
    mhW1 <- exp(mhW)
    
    h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
    w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
    phi <- 1-h-w
    
    logLikelihood <- dcat(x, prob = c(0, phi, h, w), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
  
  if(z == 3 | z == 4){
    logLikelihood <- dcat(x, prob = c(0, 0, 0, 1), log=1)
    if(log == 1){return(logLikelihood)}else{return(exp(logLikelihood))}
  }
})

#### 2.Sampling function ####
rcatHR <- nimbleFunction(run = function( 
    n = integer(0),
    z = double(0),
    gamma = double(0),
    mhH = double(0),
    mhW = double(0)
){
  # Return type declaration
  returnType(double(0))

  if(z == 1){
    state <- rcat(1, prob = c(1 - gamma, gamma))
    return(state)
  }
  
  if(z == 2){
    mhH1 <- exp(mhH)
    mhW1 <- exp(mhW)
    
    h <- (1-exp(-(mhH1+mhW1)))* (mhH1/(mhH1+mhW1))
    w <- (1-exp(-(mhH1+mhW1)))* (mhW1/(mhH1+mhW1))
    phi <- 1-h-w
    
    state <- rcat(1, prob = c(0, phi, h, w))
    return(state)
  }
  
  if(z == 3 | z == 4){
    state <- 4
    return(state)
  }
  
})

#### 3.Registration ####
registerDistributions(list(
  dcatHR = list(
    BUGSdist = "dcatHR(z, gamma, mhH, mhW)",
    types = c( "value = double(0)",
               "z = double(0)",
               "gamma = double(0)",
               "mhH = double(0)",
               "mhW = double(0)"
    ),
    pqAvail = FALSE,
    discrete = TRUE
  )))

