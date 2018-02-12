// Title:    Gibbs Sampler for MIBEN & MIBL
// Author:   Kyle M. Lang
// Created:  2014-AUG-20
// Modified: 2018-FEB-12
// Purpose:  This function will do the Gibbs sampling for the Bayesian Elastic
//           Net and Bayesian LASSO models that underlie MIBRR's core functions.

//--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------//
//  Copyright (C) 2018 Kyle M. Lang <k.m.lang@uvt.nl>                          //  
//                                                                             //
//  This file is part of MIBRR.                                                //
//                                                                             //
//  This program is free software: you can redistribute it and/or modify it    //
//  under the terms of the GNU General Public License as published by the      //
//  Free Software Foundation, either version 3 of the License, or (at you      //
//  option) any later version.                                                 //
//                                                                             //
//  This program is distributed in the hope that it will be useful, but        //
//  WITHOUT ANY WARRANTY; without even the implied warranty of                 //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General   //
//  Public License for more details.                                           //
//                                                                             //
//  You should have received a copy of the GNU General Public License along    //
//  with this program. If not, see <http://www.gnu.org/licenses/>.             //
//-----------------------------------------------------------------------------//

#include "MibrrData.h"
#include "MibrrGibbs.h"

// The following plugin should allow us to use C++11:
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List runGibbs(Eigen::MatrixXd           data,
		    Eigen::VectorXd           dataScales,
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
		    bool                      verbose,
		    bool                      doBl,
		    bool                      fullBayes,
		    bool                      adaptScales,
		    bool                      simpleIntercept,
		    bool                      noMiss,
		    std::vector<unsigned int> seeds)
{
  // Unpack the list of missing row indices:
  std::vector< std::vector<int> > missIndices;
  for(int v = 0; v < nTargets; v++) missIndices.push_back(missList[v]);

  // Initialize the various classes needed below:
  MibrrData  mibrrData(data, dataScales, missIndices, respCounts, noMiss);
  MibrrGibbs *mibrrGibbs = new MibrrGibbs[nTargets];
     
  // Initialize all parameters and setup the Gibbs sampler:
  for(int j = 0; j < nTargets; j++) {
    Eigen::VectorXd betaStartVec  = betaStarts.col(j);
    Eigen::ArrayXd  tauStartArray = tauStarts.col(j).array();

    if(doBl) mibrrGibbs[j].doBl(); // Using Bayesian LASSO?

    if(fullBayes) { // Fully Bayesian estimation?
      mibrrGibbs[j].doFullBayes(); 
      mibrrGibbs[j].setLam1Parms(l1Parms);

      if(!doBl) mibrrGibbs[j].setLam2Parms(l2Parms);
    }
    
    mibrrGibbs[j].seedRng(seeds[j]);
    
    mibrrGibbs[j].startParameters(betaStartVec,
				  tauStartArray,
				  sigmaStarts[j],
				  lambda1[j],
				  lambda2[j]);
    
    mibrrGibbs[j].setTargetIndex(j);
    mibrrGibbs[j].setDoImp(!noMiss);
    mibrrGibbs[j].setNDraws(totalSams - burnSams);
    
    if(!verbose)        mibrrGibbs[j].beQuiet(); 
    if(simpleIntercept) mibrrGibbs[j].useSimpleInt();   
  }// CLOSE for(in j ==0; j < nTargets; j++)
  
  for(int i = 0; i < totalSams; i++) {// LOOP over Gibbs iterations
    bool check0 = verbose & (i % (totalSams / 10) == 0); 
    if(check0) {
      if(i < burnSams) {
	Rcpp::Rcout << "Doing Gibbs burn-in iteration " << (i + 1);
	Rcpp::Rcout << " of " << burnSams << endl;
      }
      else {
	Rcpp::Rcout << "Doing Gibbs sampling iteration " << (i + 1) - burnSams;
	Rcpp::Rcout << " of " << totalSams - burnSams << endl;
      }
    }
    // Improve the output's aesthetics:
    bool check1 = verbose & ((i == burnSams - 1) || (i == totalSams - 1)); 
    if(check1) Rcpp::Rcout << "\n";
    
    for(int j = 0; j < nTargets; j++) {// LOOP over target variables
      // Update the Gibbs samples:
      mibrrGibbs[j].doGibbsIteration(mibrrData);
      // Start saving iterations after burn-in:
      if((i + 1) == burnSams) mibrrGibbs[j].startGibbsSampling(mibrrData);
    }
    
    if(adaptScales) mibrrData.computeDataScales();
    
  }// CLOSE for (int i = 0; i < nGibbsIters; i++)
  
  RList outList(nTargets);
  for(int j = 0; j < nTargets; j++)
    outList[j] = 
      RList::create(Rcpp::Named("imps" )  = mibrrGibbs[j].getImpSam(), 
		    Rcpp::Named("beta" )  = mibrrGibbs[j].getBetaSam(),
		    Rcpp::Named("tau"  )  = mibrrGibbs[j].getTauSam(),
		    Rcpp::Named("sigma")  = mibrrGibbs[j].getSigmaSam(),
		    Rcpp::Named("lambda") = mibrrGibbs[j].getLambdaSam()
		    );    

  return outList;
}// END runGibbs()
