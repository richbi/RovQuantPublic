// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <unordered_map>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
// QUANTILE FUNCTION 
NumericVector quantileCpp( NumericVector x,
                           NumericVector q
) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(q - 0.000000001)];
}


// [[Rcpp::export]]
List  GetDetectability_normal( NumericMatrix p0,                 // detector-specific values of p0   
                               NumericVector sigma,              // SIGMA VALUES 
                               NumericMatrix habitatxy,          // HABITAT COORDINATES
                               NumericMatrix detectorxy,         // DETECTOR COORDINATES
                               NumericVector size,               // NUMBER OF SUB-DETECTORS 
                               NumericMatrix regionID,           // MATRIX WITH REGION ID, ONE ROW PER REGION WITH 1 AND 0 WHETHER IT BELONGS TO THE REGION OR NOT NEED TO PROVIDE ROWNAMES TO IDENTIFY REGIONS.
                               NumericVector probs = NumericVector::create(0.025,0.975), //CI TO BE RETURNED
                               double localDist = 100,
                               bool display_progress = true,     // DISPLAY PROGRESS BAR
                               bool returnPosteriorCells = false // IF POSTERIORS SHOULD BE RETURNED
){
  // INITIALIZE OBJECTS 
  int niter = p0.nrow(); 
  int ncells = habitatxy.nrow();
  int ndets = detectorxy.nrow();
  int nregions = regionID.nrow();
  NumericMatrix pTot(ncells, niter);
  double maxD2 = localDist*localDist;
  
  // CALCULATE DISTANCE^2
  Progress prog1(ncells, display_progress);
  NumericMatrix dist2(ncells,ndets);
  for (int i = 0; i < ncells; i++) {
    if(Progress::check_abort())return -1.0;
    prog1.increment(); 
    for (int j = 0; j < ndets; j++) {
      dist2(i,j) = pow(habitatxy(i,0) - detectorxy(j,0),2) + pow(habitatxy(i,1) - detectorxy(j,1),2);
    }
  }
  
  
  // INITIATE PROGRESS BAR
  Progress prog2(niter, display_progress);
  
  // LOOP FOR EACH MCMC ITERATION  
  for(int ite = 0; ite < niter; ite++){
    // UPDATE PROGRESS BAR
    if(Progress::check_abort())return -1.0;
    prog2.increment(); 
    
    // PREPARE SIGMA
    double sigma1 = (2*sigma(ite)*sigma(ite));
    
    // LOOP FOR EACH HABITAT CELL
    for (int i = 0; i < ncells; i++){
      double prod = 1;
      // LOOP FOR EACH DETECTOR
      for (int j = 0; j < ndets; j++) {
        if(dist2(i,j) < maxD2){
          
          // CALCULATE 1-DETECTION PROBABILITY
          double OneMinusP = pow(1 - (p0(ite,j) * exp(-dist2(i,j)/sigma1)),size(j));
          // TAKE THE PRODUCT
          prod *= OneMinusP;
      }
      }
      pTot(i,ite) = 1-prod;
    }
  }
  
  // CELL STATISTICS 
  NumericVector outMean(ncells);
  NumericVector outMedian(ncells);
  NumericVector outsd(ncells);
  NumericVector outCV(ncells);
  NumericVector outCIL(ncells);
  NumericVector outCIH(ncells);
  
  for (int i = 0; i < ncells; i++){
    NumericVector tmpMean = pTot(i,_);
    outMean(i) = mean(tmpMean);
    outMedian(i) = median(tmpMean);
    outsd(i) = sd(tmpMean);
    outCV(i) = (100*outsd(i))/outMean(i);
    NumericVector outCILr= quantileCpp(tmpMean ,probs);
    outCIL(i) = outCILr(0);
    outCIH(i) = outCILr(1);
  }
  
  //================
  //==== RETURN ====
  //================
  // INITIALIZE OBJECTS
  NumericMatrix subsetSum(nregions, niter);
  NumericMatrix summary(nregions+1,4);
  // ATTRIBUTE NAMES FOR THE SUMMARY TABLE 
  CharacterVector Names = rownames(regionID);
  // ADD A TOTAL ROW
  CharacterVector Names1 = CharacterVector::create("Total");
  CharacterVector Names2(Names.size() + Names1.size());
  std::copy(Names.begin(), Names.end(), Names2.begin());
  std::copy(Names1.begin(), Names1.end(), Names2.begin() + Names.size());
  rownames(summary) = Names2;
  colnames(summary) = CharacterVector::create("mean", "median", "95%CILow","95%CIHigh");
  
  // SUMMARY AND POSTERIOR FOR EACH REGIONS 
  for (int r = 0; r < nregions; r++){
    // SUBSET TO CELL WITHIN REGIONS 
    NumericVector T = regionID(r,_);
    colvec tIdx(T.begin(), T.size(), false); 
    mat Xmat(pTot.begin(), pTot.nrow(), pTot.ncol(), false);
    mat subMat = Xmat.rows(find(tIdx == 1));
    int nsub = sum(T);
    // AVERAGE THE DETECTION PROBABILITY OVER ALL CELLS OF THE REGION
    for (int ite = 0; ite < niter; ite++){
      double total = 0;
      for (int j = 0; j < nsub; j++) {
       double meanP = subMat(j,ite)/nsub;
        total += meanP;
      }
      subsetSum(r,ite) = total;
    }
    // FILL IN THE SUMMARY TABLE
    NumericVector tmp = subsetSum(r,_);
    summary(r,0) = mean(tmp);
    summary(r,1) = median(tmp);
    NumericVector outCILr = quantileCpp(tmp, probs);
    summary(r,2) = outCILr(0);
    summary(r,3) = outCILr(1);
  }
  
  // SUMMARY AND POSTERIOR FOR ALL REGIONS
  // SUBSET TO CELLS WITHIN ALL REGIONS 
  NumericVector AllRegionsID = colSums(regionID);
  mat Xmat1(pTot.begin(), pTot.nrow(), pTot.ncol(), false);
  colvec tIdx1(AllRegionsID.begin(), AllRegionsID.size(), false); 
  mat subMat1 = Xmat1.rows(find(tIdx1 > 0));
  int nsub1 = sum(AllRegionsID);
  
  // AVERAGE THE DETECTION PROBABILITY OVER ALL CELLS
  NumericVector subsetSum1(niter) ;
  for (int ite = 0; ite < niter; ite++){
    double total = 0;
    for (int j = 0; j < nsub1; j++) {
      total += subMat1(j, ite)/nsub1;
    }
    subsetSum1(ite) = total;
  }
  
  // FILL IN THE "Total" SUMMARY TABLE
  summary(nregions,0) = mean(subsetSum1);
  summary(nregions,1) = median(subsetSum1);
  NumericVector outCILr1 = quantileCpp(subsetSum1, probs);
  summary(nregions,2) = outCILr1(0);
  summary(nregions,3) = outCILr1(1);
  
  // OUTPUT
  if(returnPosteriorCells){
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum,
                        Named("PosteriorCells") = pTot);
  } else {
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum);
  }
}


// [[Rcpp::export]]
List  GetDetectability_mean( NumericVector p0,             // detector-specific values of p0   
                             double sigma,                 // SIGMA VALUES 
                             NumericMatrix habitatxy,      // HABITAT COORDINATES
                             NumericMatrix detectorxy,     // DETECTOR COORDINATES
                             NumericVector size,           // NUMBER OF SUB-DETECTORS
                             bool display_progress = true  // DISPLAY PROGRESS BAR
){
  // INITIALIZE OBJECTS 
  int ncells = habitatxy.nrow();
  int ndets = detectorxy.nrow();
  NumericVector pTot(ncells);
  double sigma1 = 2*sigma*sigma;
  
  // INITIATE PROGRESS BAR
  Progress prog2(ncells, display_progress);
  
  // LOOP FOR EACH HABITAT CELL
  for (int i = 0; i < ncells; i++){
    // UPDATE PROGRESS BAR
    if(Progress::check_abort())return -1.0;
    prog2.increment(); 
    
    double prod = 1;
    
    // LOOP FOR EACH DETECTOR
    for (int j = 0; j < ndets; j++) {
      // CALCULATE SQUARED DISTANCE
      double d2 = pow(habitatxy(i,0) - detectorxy(j,0),2) + pow(habitatxy(i,1) - detectorxy(j,1),2);
      
      // CALCULATE 1-DETECTION PROBABILITY
      double OneMinusP = pow(1 - (p0(j) * exp(-d2/sigma1)),size(j));
      
      // TAKE THE PRODUCT
      prod *= OneMinusP;
    }
    pTot(i) = 1-prod;
  }
  
  // OUTPUT
  return List::create(Named("MeanCell") = pTot);
}