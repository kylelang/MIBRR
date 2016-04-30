// Title:    Header file for the MyGibbs Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2016-APR-29
// Purpose:  This class contains the Gibbs sampling-related functions
//           for the MIBRR package.

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

#ifndef MYGIBBS_H
#define MYGIBBS_H

#include "MyErrors.hpp"
#include "MyData.hpp"
#include <nlopt.hpp>

using namespace std;
using namespace Eigen;

class MyGibbs: public MyErrors {

 public:

  ///// CONSTRUCTORS / DESTRUCTOR /////

  MyGibbs();
  
  ~MyGibbs();


  ///// ACCESSORS /////

  VectorXd getBetas() const;
  // @return: a vector of regression coefficients (and the intercept term)
  //          for the elastic net model
  
  ArrayXd getTaus() const;
  // @return: an array of auxiliary penalty hyperparameters 
  //          for the elastic net model

  double getSigma() const;
  // @return: the residual variance of the elastic net model
  
  MatrixXd getBetaSam() const; 
  // @return: a matrix containing the Gibbs sample of beta
  
  ArrayXXd getTauSam() const; 
  // @return: a two-dimentional array containing the Gibbs sample of tau
  
  VectorXd getSigmaSam() const;
  // @return: a vector containing the Gibbs sample of sigma
  
  MatrixXd getImpSam() const;
  // @return: a matrix containing the Gibbs sample of the DV
  
  MatrixXd getLambdaHistory() const;
  // @return: a matrix containing the estimates of Lambda at each iteration
  //          of the MCEM algorithm
  
  VectorXd getLambdas() const;
  // @return: the current values of the penalty parameters
  
  double getLambdas(int) const;
  // @param:  integer telling which lambda to return
  //          (1 give "LASSO" version, 2 gives "ridge" version)
  // @return: the current value of either lambda1 and lambda2

  int getNDraws() const;
  // @return: the current number of retained Gibbs sampling draws
  

  ///// MUTATORS /////
  
  void setBetas(VectorXd&);
  // @param: a new vector of coefficients for the elastic net model
  
  void setTaus(ArrayXd&);
  // @param: a new array of hyperparameters for the elastic net model
  
  void setSigma(double);
  // @param: a new residual variance of the elastic net model
  
  void setLambdas(VectorXd&);
  // @param: a new value for Lambda
  
  void setLambdas(double, double);
  // @param1: a new value for lambda1
  // @param2: a new value for lambda2
  
  void setLambdas(double);
  // @param1: a new value for the LASSO lambda

  void setLambdas();
  // @effect: set lambda1 = XX, lambda2 = XX
  
  void setNDraws(int);
  //@param: a new value for the retained Gibbs sampling draws

  void setNEmIters(int);
  //@param: a new value for the number of MCEM iterations

  void setTargetIndex(int);
  // @param: a new value for the target variable's column index

  void setVerbosity(bool, bool);
  // @param1: logical switch indicating if verbose iteration history should be
  //          reported
  // @param2: logical switch indicating if verbose errors be reported
  
  void doMiben();
  // @effect: set the imputation model to the Bayesian elastic net 
  
  void doMibl();
  // @effect: set the imputation model to the Bayesian LASSO 
  
  void startGibbsSampling(MyData&);
  // @effect: start storing the parameters' Gibbs sampled values
 
  void stopGibbsSampling();
  // @effect: stop storing the parameters' Gibbs sampled values
  
  void setupOptimizer(int,
		      int,
		      double);
  // @param1: total number of iterations in the MCEM algorithm
  // @param2: number of tuning iterations in the MCEM algorithm
  // @param3: convergence criterion for the MCEM optimizers
  // @effect: parameterize the MIBEN EM optimization
  
  //void setupOptimizer(int, int);
  // @param1: total number of iterations in the MCEM algorithm
  // @param2: number of tuning iterations in the MCEM algorithm
  // @effect: parameterize the MIBL EM optimization

  void startParameters(VectorXd&,
		       ArrayXd&,
		       double,
		       double,
		       double);
  // @param1: starting values for beta
  // @param2: starting values for tau
  // @param3: starting value for sigma
  // @param4: starting value for lambda1
  // @param5: starting value for lambda2
  // @effect: provide starting values for all model parameters
  
  void startParameters(VectorXd&,
		       ArrayXd&,
		       double,
		       VectorXd&);
  // @param1: starting values for beta
  // @param2: starting values for tau
  // @param3: starting value for sigma
  // @param4: starting value for Lambda
  // @effect: provide starting values for all model parameters
  
  void restartParameters(MyData&);
  // @param:  an initialized MyData object
  // @effect: restart all model parameters at the posterior means of their
  //          respective Gibbs samples

  
  ///// PARAMETER UPDATE FUNCTIONS /////

  void updateTaus(MyData&);
  // @param:  an initialized MyData object
  // @effect: update _taus (for the Bayesian Elastic Net) based on 
  // current values of other member variables

  void updateBetas(MyData&);
  // @param:  an initialized MyData object
  // @effect: update _betas based on current values of other member variables

  void updateSigma(MyData&);
  // @param:  an initialized MyData object
  // @effect: update _sigma (for the Bayesian Elastic Net) based on 
  //          current values of other member variables
  
  void updateImputations(MyData&);
  // @param:  an initialized MyData object
  // @effect: update the imputations based on current values of member variables
  
  void doGibbsIteration(MyData&);
  // @param:  an initialized MyData object
  // @effect: run a single iteration of the Gibbs sampler
  

  ///// MCEM OPTIMIZATION FUNCTIONS //////

  double eNetLambdaObjective(const std::vector<double>&,
			     std::vector<double>&,
			     void*);
  // @param1: starting values for Lambda values
  // @param2: starting values for the gradient
  // @param3: pointer to additional data
  // @return: the evaluated EM objective function (with gradient) for the 
  //          Elastic Net penalty parameters
  
  void optimizeMibenLambdas(bool);
  // @param:  logical switch indicated whether we are pre-optimizing (true) or
  //          fully optimizing (false) Lambda
  // @effect: Numerically optimize the MIBEN penalty parameters
  
  void updateLambdas();
  // @effect: update the ENET or LASSO penalty parameters using marginal,
  //          numerical/deterministic optimization.
  
 private:
  
  VectorXd _betas;
  ArrayXd  _taus;
  double   _sigma;
  MatrixXd _betaSam;
  ArrayXXd _tauSam;
  VectorXd _sigmaSam;
  MatrixXd _impSam;
  MatrixXd _lambdaHistory;
  VectorXd _lambdas;
  double   _emConvTol;
  int      _nDraws;
  int      _nEmIters;
  int      _targetIndex;
  int      _drawNum;
  int      _emIterNum;
  int      _optIterCount;
  int      _lambdaWindow;
  bool     _verboseIters;
  bool     _storeGibbsSamples;
};

#endif
