// Title:    Header file for the MibrrGibbs Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2017-SEP-30
// Purpose:  This class contains the Gibbs sampling-related functions for the
//           MIBRR package.

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

#ifndef MIBRRGIBBS_H
#define MIBRRGIBBS_H

#include "MibrrData.hpp"

using namespace std;
using namespace Eigen;

class MibrrGibbs {

public:
  ///////////////////////// CONSTRUCTORS / DESTRUCTOR /////////////////////////////

  MibrrGibbs();

  ~MibrrGibbs();

  //////////////////////////////// ACCESSORS //////////////////////////////////////

  VectorXd getBetas() const;
  // @return: regression coefficients (and intercept) for the elastic net model
  
  ArrayXd getTaus() const;
  // @return: auxiliary penalty hyperparameters for the elastic net model

  double getSigma() const;
  // @return: residual variance of the elastic net model

  MatrixXd getBetaSam() const; 
  // @return: (burnt in) Gibbs sample of beta

  ArrayXXd getTauSam() const; 
  // @return: (burnt in) Gibbs sample of tau

  VectorXd getSigmaSam() const;
  // @return: (burnt in) Gibbs sample of sigma

  MatrixXd getImpSam() const;
  // @return: (burnt in) Gibbs sample of the DV

  //MatrixXd getLambdaHistory() const;
  // @return: estimates of Lambda at each iteration of the MCEM algorithm

  int getNDraws() const;
  // @return: current number of retained Gibbs sampling draws

  //int getNEmIters() const;
  // @return: number of MCEM iterations requested

  bool getVerbosity() const;
  // @return: flag indicating if printed output is currently verbose

  bool getElasticNetFlag() const;
  // @return: current value of the _useElasticNet switch

  bool getDoImp() const;
  // @return: current value of the flag denoting if missing data are to be
  //          imputed

  bool getSimpleIntercept() const;
  // @return: current value of the flag denoting if the unconditional mean of the
  //          target variable should be used as the mean of the intercept's
  //          posterior (as opposed to using the conditional mean of the target)

  VectorXd getLambdas() const;
  // @return: current values of the penalty parameters (Lambda)

  double getLambdas(int) const;
  // @param:  which lambda to return (1 = "LASSO", 2 = "ridge")
  // @return: current value of either lambda1 and lambda2

  
  //////////////////////////////// MUTATORS /////////////////////////////////////


  void setBetas(VectorXd&);
  // @param: new coefficients for the elastic net model

  void setTaus(ArrayXd&);
  // @param: new hyperparameters for the elastic net model

  void setSigma(double);
  // @param: new residual variance of the elastic net model

  //void setNEmIters(int);
  // @param: new value for the number of MCEM iterations

  void setTargetIndex(int);
  // @param: new value for the target variable's column index

  void setNDraws(int);
  // @param: new value for the number of retained Gibbs sampling draws

  void setDoImp(bool);
  // @param: new value for the logical switch for imputation
  
  void beQuiet();
  // @effect: turn off verbose output

  void doBl();
  // @effect: set the imputation model to the Bayesian LASSO 

  void useSimpleInt();
  // @effect: use the unconditional mean of the target as the mean of the
  //          intercept's posterior

  void setLambdas(VectorXd&);
  // @param: new value for Lambda

  void setLambdas(double, double);
  // @param1: new value for lambda1
  // @param2: new value for lambda2

  void setLambdas(double);
  // @param1: new value for the LASSO lambda

  //void setLambdas();
  // @effect: set lambdas to the average of their values within _lambdaWindow

  void startParameters(VectorXd&, ArrayXd&, double, double, double);
  // @param1: starting values for beta
  // @param2: starting values for tau
  // @param3: starting value for sigma
  // @param4: starting value for lambda1
  // @param5: starting value for lambda2
  // @effect: provide starting values for all model parameters

  void startParameters(VectorXd&, ArrayXd&, double, VectorXd&);
  // @param1: starting values for beta
  // @param2: starting values for tau
  // @param3: starting value for sigma
  // @param4: starting value for Lambda
  // @effect: provide starting values for all model parameters

  //void setupOptimizer(int, int, double, bool);
  // @param1: number of MCEM iterations
  // @param2: number of iterations in the smoothing window
  // @param3: convergence criterion for the optimization methods
  // @param4: logical switch for two-phase optimization for MIBEN's lambdas
  // @effect: parameterize the optimizers used to estimate Lambda

  void startGibbsSampling(const MibrrData&);
  // @effect: start storing the parameters' Gibbs sampled values

  void restartParameters(MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: restart all model parameters at the posterior means of their
  //          respective Gibbs samples


  ////////////////////////// RANDOM VARIATE SAMPLERS ////////////////////////////


  double drawInvGamma(double, double) const;
  // @param1: shape parameter
  // @param2: scale parameter
  // @return: random Inverse Gamma variate

  double calcIncGamma(const double, const double, const bool);
  // @param1: shape parameter of the underlying gamma distribution
  // @param2: threshold value cutting off the upper or lower tail
  // (i.e., the underlying variate whose probability or [1 - probability] 
  // is being returned)
  // @param3: true = lower incomplete gamma, false = upper incomplete gamma
  // @return: value of the incomplete gamma function (i.e., the un-normalized 
  // area under the upper or lower tail of the Gamma CDF).

  // The following function was adapted from source code originally
  // implemented by Robert E. Wheeler (2001-MAR) in the R package SuppDists:
  double drawInvGauss(const double, const double);
  // @param1: mean parameter (mu)
  // @param2: shape parameter (lambda)
  // @return: random variate from the inverse Gaussian distribution


  ///////////////////////// PARAMETER UPDATE FUNCTIONS //////////////////////////


  void updateTaus(const MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: update _taus based on current values of other member variables

  void updateBetas(const MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: update _betas based on current values of other member variables

  void updateSigma(const MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: update _sigma based on current values of other member variables
  
  void updateImputations(MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: update the imputations based on current values of member variables
  
  void doGibbsIteration(MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: run a single iteration of the Gibbs sampler


  //////////////////////// MCEM OPTIMIZATION FUNCTIONS //////////////////////////

  //// Beginning with Version 0.0.0.9001, the MCEM steps have been moved back to
  //// the R layer because installing NLopt has proven to be too much of a
  //// portability issue.
  
  //double lambdaObjective(const std::vector<double>&,
  //			 std::vector<double>&, void*);
  // @param1: starting values for Lambda
  // @param2: starting values for Lambda's gradient
  // @param3: pointer to additional data
  // @return: the BEN penalty parameters' evaluated EM objective function 

  //void optimizeMibenLambdas(const bool);
  // @param:  are pre-optimizing (true) or fully optimizing (false) Lambda?
  // @effect: numerically optimize the MIBEN penalty parameters

  //void updateLambdas();
  // @effect: update the BEN or LASSO penalty parameters using marginal,
  //          numerical/deterministic optimization.

  
  //////////////////////// EXCEPTION HANDLING FUNCTIONS /////////////////////////

  
  void tauError(int) const;
  // @param:  the orignial error code
  // @effect: dispatch an appropriate error message to stderr

  void betaError(exception&) const;
  // @param:  orignial exception object
  // @effect: dispatch an appropriate error message to stderr

  //void lambdaError() const;
  // @effect: dispatch an appropriate error message to stderr

  //void lambdaError(exception&) const;
  // @param1: exception thrown by nlopt
  // @effect: dispatch an appropriate error message to stderr
  
private:
  VectorXd _betas;
  ArrayXd  _taus;
  double   _sigma;
  VectorXd _lambdas;
  MatrixXd _betaSam;
  ArrayXXd _tauSam;
  VectorXd _sigmaSam;
  MatrixXd _impSam;
  //MatrixXd _lambdaHistory;
  //double   _emConvTol;
  int      _targetIndex;
  int      _nDraws;
  int      _drawNum;
  //int      _nEmIters;
  //int      _emIterNum;
  //int      _optIterCount;
  //int      _lambdaWindow;
  bool     _verbose;
  bool     _useElasticNet;
  bool     _storeGibbsSamples;
  bool     _doImp;
  bool     _simpleIntercept;
  //bool     _twoPhaseOpt;
  //string   _optPrefix;
  //string   _algName;
  //int      _optMethod;
};

#endif
