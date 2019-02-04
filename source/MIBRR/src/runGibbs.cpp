// Title:    Gibbs Sampler for MIBEN & MIBL
// Author:   Kyle M. Lang
// Created:  2014-AUG-20
// Modified: 2019-JAN-18
// Purpose:  This function will do the Gibbs sampling for the Bayesian Elastic
//           Net and Bayesian LASSO models that underlie MIBRR's core functions.

//--------------------- COPYRIGHT & LICENSING INFORMATION --------------------//
//  Copyright (C) 2019 Kyle M. Lang <k.m.lang@uvt.nl>                         //
//                                                                            //
//  This file is part of MIBRR.                                               //
//                                                                            //
//  This program is free software: you can redistribute it and/or modify it   //
//  under the terms of the GNU General Public License as published by the     //
//  Free Software Foundation, either version 3 of the License, or (at you     //
//  option) any later version.                                                //
//                                                                            //
//  This program is distributed in the hope that it will be useful, but       //
//  WITHOUT ANY WARRANTY; without even the implied warranty of                //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General  //
//  Public License for more details.                                          //
//                                                                            //
//  You should have received a copy of the GNU General Public License along   //
//  with this program. If not, see <http://www.gnu.org/licenses/>.            //
//----------------------------------------------------------------------------//

#include "MibrrData.h"
#include "MibrrGibbs.h"

// The following plugin should allow us to use C++11:
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List runGibbs(Eigen::MatrixXd           data,
		    int                       nTargets,  
		    Rcpp::List                missList,
		    Eigen::VectorXi           respCounts,
		    Eigen::VectorXd           lambda1,
		    Eigen::VectorXd           lambda2,
		    Eigen::VectorXd           l1Parms,
		    Eigen::VectorXd           l2Parms,
		    Eigen::VectorXd           sigmaStarts,
		    Eigen::MatrixXd           tauStarts,
		    Eigen::MatrixXd           betaStarts,
		    int                       burnSams,
		    int                       totalSams,
		    int                       penType,
		    double                    ridge,
		    bool                      verbose,
		    bool                      fullBayes,
		    bool                      noMiss,
		    bool                      savePpSams,
		    std::vector<unsigned int> seeds)
{
  // Unpack the list of missing row indices:
  std::vector< std::vector<int> > missIndices;
  for(int v = 0; v < nTargets; v++) missIndices.push_back(missList[v]);
  
  // Initialize the various classes needed below:
  MibrrData  mibrrData(data, missIndices, respCounts, noMiss);
  MibrrGibbs *mibrrGibbs = new MibrrGibbs[nTargets];
  
  // Initialize all parameters and setup the Gibbs sampler:
  for(int j = 0; j < nTargets; j++) {
    Eigen::VectorXd betaStartVec  = betaStarts.col(j);
    Eigen::ArrayXd  tauStartArray = tauStarts.col(j).array();
  
    // Define the type of regularization:
    mibrrGibbs[j].setPenType(penType);
  
    if(penType == 0)// Doing basic ridge?
      mibrrGibbs[j].setRidge(ridge);
 
    if(fullBayes) {// Fully Bayesian estimation?
      mibrrGibbs[j].doFullBayes(); 
      mibrrGibbs[j].setLam1Parms(l1Parms);
  
      if(penType == 2)// Doing MIBEN?
	mibrrGibbs[j].setLam2Parms(l2Parms);
    }
    
    if(savePpSams) mibrrGibbs[j].savePpSams();
    
    mibrrGibbs[j].seedRng(seeds[j]);
 
    mibrrGibbs[j].startParameters(betaStartVec,
				  tauStartArray,
				  sigmaStarts[j],
				  lambda1[j],
				  lambda2[j]);
  
    mibrrGibbs[j].setTargetIndex(j);
    mibrrGibbs[j].setDoImp(!noMiss);
    mibrrGibbs[j].setNDraws(totalSams - burnSams);
     
    if(!verbose) mibrrGibbs[j].beQuiet(); 
  }// CLOSE for(in j ==0; j < nTargets; j++)
  
  for(int i = 0; i < totalSams; i++) {// LOOP over Gibbs iterations
    // Print a nice progress message:
    if(verbose) {
      int marg, max;
      bool check0;
      if(i < burnSams) {
	marg   = burnSams % 5;
	max    = burnSams - marg;
	check0 = (i % (max / 5) == 0) & ((burnSams - i) > marg);
	if(check0) {
	  Rcpp::Rcout << "Doing Gibbs burn-in iteration " << (i + 1);
	  Rcpp::Rcout << " of " << burnSams << endl;
	}
      }
      else {
	marg   = (totalSams - burnSams) % 5;
	max    = (totalSams - burnSams) - marg;
	check0 = ((i - burnSams) % (max / 5) == 0) & ((totalSams - i) > marg);
	if(check0) {
	  Rcpp::Rcout <<
	    "Doing Gibbs sampling iteration " << (i + 1) - burnSams;
	  Rcpp::Rcout << " of " << totalSams - burnSams << endl;
	}
      }
    }
   
    // Improve the output's aesthetics:
    bool check1 = verbose & ((i == burnSams - 1) || (i == totalSams - 1)); 
    if(check1) Rcpp::Rcout << "\n";
    
    for(int j = 0; j < nTargets; j++) {// LOOP over target variables
      // Compute new centers and scales for the jth target's predictors:
      mibrrData.updateMoments(j);

      // Update the Gibbs samples:
      mibrrGibbs[j].doGibbsIteration(mibrrData);
      
      // Start saving iterations after burn-in:
      if((i + 1) == burnSams) mibrrGibbs[j].startGibbsSampling(mibrrData);
    }   
  }// CLOSE for (int i = 0; i < totalSams; i++)
  
  RList outList(nTargets);
  for(int j = 0; j < nTargets; j++)
    outList[j] = 
      RList::create(Rcpp::Named("imps"  ) = mibrrGibbs[j].getImpSam(),
		    Rcpp::Named("ppSams") = mibrrGibbs[j].getPpSam(),
		    Rcpp::Named("beta"  ) = mibrrGibbs[j].getBetaSam(),
		    Rcpp::Named("tau"   ) = mibrrGibbs[j].getTauSam(),
		    Rcpp::Named("sigma" ) = mibrrGibbs[j].getSigmaSam(),
		    Rcpp::Named("lambda") = mibrrGibbs[j].getLambdaSam()
		    );    

  return outList;
}// END runGibbs()
