// Title:    Gibbs Sampler for MIBEN & MIBL
// Author:   Kyle M. Lang
// Created:  2014-AUG-20
// Modified: 2016-MAY-11
// Purpose:  This code is part of the R package mibrr.
//           This function will do the Gibbs sampling for Multiple Imputation
//           with the Bayesian Elastic Net (MIBEN) and Multiple Impution with
//           the Bayesian LASSO (MIBL).

//--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------//
//  Copyright (C) 2016 Kyle M. Lang <kyle.lang@ttu.edu>                        //  
//                                                                             //
//  This file is part of mibrr.                                                //
//                                                                             //
//  This program is free software: you can redistribute it and/or modify it    //
//  under the terms of the GNU Lesser General Public License as published by   //
//  the Free Software Foundation, either version 3 of the License, or          //
//  (at you option) any later version.                                         //
//                                                                             //
//  This program is distributed in the hope that it will be useful, but        //
//  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY //
//  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public    //
//  License for more details.                                                  //
//                                                                             //
//  You should have received a copy of the GNU Lesser General Public License   //
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.      //
//-----------------------------------------------------------------------------//

#include <RcppEigen.h>
#include "MibrrDefs.hpp"
#include "MibrrData.hpp"
#include "MibrrGibbs.hpp"

// [[Rcpp::export]]
Rcpp::List runGibbs(Eigen::MatrixXd inData,
		    Eigen::VectorXd dataScales,
		    int             nTargets,  
		    Eigen::VectorXd lambda1Starts,
		    Eigen::VectorXd lambda2Starts,
		    Eigen::VectorXd sigmaStarts,
		    Eigen::MatrixXd tauStarts,
		    Eigen::MatrixXd betaStarts,
		    double          missCode,
		    int             nMcemApproxIters,
		    int             nMcemTuneIters,
		    int             nMcemApproxBurn,
		    int             nMcemApproxGibbs,
		    int             nMcemTuneBurn,		   
		    int             nMcemTuneGibbs,
		    int             nMcemPostBurn,
		    int             nMcemPostGibbs,
		    int             lambdaWindow,
		    double          emConvTol,
		    bool            verboseIters,
		    bool            verboseErrors,
		    bool            doMibl,
		    bool            regIntercept)
{
  // Initialize the various classes needed below:
  MibrrData mibrrData(inData, 
		      dataScales,
		      missCode);
  MibrrGibbs *mibrrGibbs = new MibrrGibbs[nTargets];
  
  // Specify some useful constants:
  int nPreds = mibrrData.nPreds();
  int nObs = mibrrData.nObs();
  int nTotalEmIters = nMcemApproxIters + nMcemTuneIters;
   
  // Initialize all parameters, setup the Gibbs sampler and the EM Optimization:
  mibrrData.fillMissing(nTargets);
  for(int j = 0; j < nTargets; j++) {
    Eigen::VectorXd betaStartVec = betaStarts.col(j);
    Eigen::ArrayXd tauStartArray = tauStarts.col(j).array();

    if(doMibl) {
      mibrrGibbs[j].doMibl();
      emConvTol = 0.0; // Set dummy value for EM convergence criterion
    }
    
    // Must call setupOptimizer() before startParameters()!
    mibrrGibbs[j].setupOptimizer(nTotalEmIters,
				 lambdaWindow,
				 emConvTol);
    
    mibrrGibbs[j].startParameters(betaStartVec,
				  tauStartArray,
				  sigmaStarts[j],
				  lambda1Starts[j],
				  lambda2Starts[j]);
    
    mibrrGibbs[j].setVerbosity(verboseIters, verboseErrors);
    mibrrGibbs[j].setTargetIndex(j);
    mibrrGibbs[j].setRegInterceptFlag(regIntercept);
  }
  
  // Specify containers for the parameters' starting values:
  Eigen::MatrixXd dvStartMat = MatrixXd(nObs, nTargets);

  for (int k = 0; k < (nTotalEmIters + 1); k++) {// LOOP over MCEM iterations
    int emIterNum = k + 1;
    int nGibbsIters, nBurnIns, nPosteriorBurnIns;
    if(k < nMcemApproxIters) {
      // Small gibbs samples for EM burn in
      nGibbsIters = nMcemApproxGibbs + nMcemApproxBurn;
      nBurnIns = nMcemApproxBurn;
    } 
    else if((k >= nMcemApproxIters) & (k < nTotalEmIters)) {
      // Large gibbs samples for EM tuning phase
      nGibbsIters = nMcemTuneGibbs + nMcemTuneBurn;
      nBurnIns = nMcemTuneBurn;
    } 
    else {
      // Large final gibbs sample from the convergent model
      nPosteriorBurnIns = 
	nMcemPostBurn < 0 ? nMcemTuneBurn : nMcemPostBurn;
      nGibbsIters = nPosteriorBurnIns + nMcemPostGibbs;
      nBurnIns = nPosteriorBurnIns;
    }
    
    if(verboseIters) {
      if(k < nMcemApproxIters) {
	Rcpp::Rcout << "\nDoing MCEM approximation iteration " << emIterNum;
	Rcpp::Rcout << " of " << nMcemApproxIters << "\n" << endl;
      } 
      else if(k == nTotalEmIters) {
	Rcpp::Rcout << "\nSampling from the stationary posterior\n" << endl;
      }
      else {
	Rcpp::Rcout << "\nDoing MCEM tuning iteration ";
	Rcpp::Rcout << emIterNum - nMcemApproxIters;
	Rcpp::Rcout << " of " << nMcemTuneIters << "\n" << endl;
      }
    }
    
    for(int i = 0; i < nGibbsIters; i++) {// LOOP over Gibbs iterations
      if(verboseIters & (i % (nGibbsIters / 10) == 0)) {
	int iterOut = i + 1;
	Rcpp::Rcout << "Doing Gibbs iteration " << (i + 1);
	Rcpp::Rcout << " of " << nGibbsIters << endl;
	// Improve the output's aesthetics:
	if(i == nGibbsIters - (nGibbsIters / 10)) Rcpp::Rcout << "\n";
      }
      
      for(int j = 0; j < nTargets; j++) {// LOOP over target variables
	bool lastEmApprox = k == nMcemApproxIters;
	if(lastEmApprox) mibrrGibbs[j].setLambdas();
	
	bool changeNDraws = 
	  (i == 0) & ((k == 0) | lastEmApprox | (k == nTotalEmIters));
	if(changeNDraws) mibrrGibbs[j].setNDraws(nGibbsIters - nBurnIns);
	
	// Update the Gibbs samples:
	mibrrGibbs[j].doGibbsIteration(mibrrData);
	
	if ((i + 1) == nBurnIns) mibrrGibbs[j].startGibbsSampling(mibrrData);
       	
	if((k < nTotalEmIters) & ((i + 1) == nGibbsIters)) {
	  mibrrGibbs[j].updateLambdas(); // Optimize the penalty parameters
	  mibrrGibbs[j].restartParameters(mibrrData);
	}
      }// CLOSE for (int j = 0; j < nTargets, j++)  
    }// CLOSE for (int i = 0; i < nGibbsIters; i++)
  }// END for(int k = 0; k < nTotalEmIters; k++)
  
  RList outList(nTargets);
  for(int j = 0; j < nTargets; j++) {
    outList[j] = 
      RList::create(Rcpp::Named("imps"         ) = mibrrGibbs[j].getImpSam(), 
		    Rcpp::Named("beta"         ) = mibrrGibbs[j].getBetaSam(),
		    Rcpp::Named("tau"          ) = mibrrGibbs[j].getTauSam(),
		    Rcpp::Named("sigma"        ) = mibrrGibbs[j].getSigmaSam(),
		    Rcpp::Named("lambdaHistory") =
		    mibrrGibbs[j].getLambdaHistory()
		    );    
  }
  return outList;
}// END runGibbs() 


