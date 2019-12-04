// Title:    Header file for the MibrrGibbs Class
// Author:   Kyle M. Lang
// Created:  2014-12-04
// Modified: 2019-12-04
// Purpose:  This class contains the Gibbs sampling-related functions for the
//           MIBRR package.

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

#ifndef MIBRRGIBBS_H
#define MIBRRGIBBS_H

#include "MibrrData.h"
#include "MibrrSamplers.h"

class MibrrGibbs: public MibrrSamplers {

public:
  //////////////////////// CONSTRUCTORS / DESTRUCTOR ///////////////////////////
  
  MibrrGibbs();
  
  ~MibrrGibbs();
  
  /////////////////////////////// ACCESSORS ////////////////////////////////////

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
  // @return: (burnt in) Gibbs sample of Y_miss

  MatrixXd getPpSam() const;
  // @return: (burnt in ) posterior predictive sample of Y_obs

  MatrixXd getLambdaSam() const;
  // @return: (burnt in) Gibbs sample of the penalty parameters
  
  int getNDraws() const;
  // @return: current number of retained Gibbs sampling draws

  int getPenType() const;
  // @return: an integer code indictating the type of regularization being used
  //          0 := No shrinkage priors (i.e., basic ridge if _ridge > 0.0)
  //          1 := Bayesian elastic net
  //          2 := Bayesian LASSO

  double getRidge() const;
  // @return: the current value of the ridge penalty
  
  bool getVerbosity() const;
  // @return: flag indicating if printed output is currently verbose
  
  bool getDoImp() const;
  // @return: current value of the flag denoting if missing data are to be
  //          imputed

  VectorXd getLambdas() const;
  // @return: current values of the penalty parameters (Lambda)

  double getLambdas(int) const;
  // @param:  which lambda to return (1 = "LASSO", 2 = "ridge")
  // @return: current value of either lambda1 and lambda2

  
  //////////////////////////////// MUTATORS ////////////////////////////////////


  void setBetas(VectorXd&);
  // @param: new coefficients for the elastic net model

  void setTaus(ArrayXd&);
  // @param: new hyperparameters for the elastic net model

  void setSigma(double);
  // @param: new residual variance of the elastic net model

  void setTargetIndex(int);
  // @param: new value for the target variable's column index

  void setNDraws(int);
  // @param: new value for the number of retained Gibbs sampling draws

  void setDoImp(bool);
  // @param: new value for the logical switch for imputation

  void beQuiet();
  // @effect: turn off verbose output

  void doFullBayes();
  // @effect: set the estimation method to fully Bayesian Gibbs sampling
  
  void savePpSams();
  // @effect: set flag to save posterior predictive samples

  void setPenType(int);
  // @param: a new value for the integer code of regularization type

  void setRidge(double);
  // @param: a new value for the ridge penalty
  
  void setLam1Parms(VectorXd&);
  // @param: new value of lambda1 prior parameters
  
  void setLam2Parms(VectorXd&);
  // @param: new value of lambda2^2 prior parameters
  
  void setLambdaParms(VectorXd&);
  // @param: concatenated values for lambda1 & lambda2^2 prior parameters

  void setFinalRep(bool);
  // @param: new value for the finalRep flag

  void setLambdas(VectorXd&);
  // @param: new value for Lambda

  void setLambdas(double, double);
  // @param1: new value for lambda1
  // @param2: new value for lambda2

  void setLambdas(double);
  // @param1: new value for the LASSO lambda

  void startParameters(VectorXd&, ArrayXd&, double, double, double, bool);
  // @param1: starting values for beta
  // @param2: starting values for tau
  // @param3: starting value for sigma
  // @param4: starting value for lambda1
  // @param5: starting value for lambda2
  // @param6: should we use the means of beta to update sigma?
  // @effect: provide starting values for all model parameters

  void startParameters(VectorXd&, ArrayXd&, double, VectorXd&, bool);
  // @param1: starting values for beta
  // @param2: starting values for tau
  // @param3: starting value for sigma
  // @param4: starting value for Lambda
  // @param6: should we use the means of beta to update sigma?
  // @effect: provide starting values for all model parameters

  void startGibbsSampling(const MibrrData&);
  // @effect: start storing the parameters' Gibbs sampled values

  
  ///////////////////////// PARAMETER UPDATE FUNCTIONS /////////////////////////

  void updateLambdas(const MibrrData&);
  // @effect: update _lambdas based on current values of other member variables
  
  void updateTaus(const MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: update _taus based on current values of other member variables

  void updateBetas(MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: update _betas based on current values of other member variables

  void updateSigma(MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: update _sigma based on current values of other member variables
  
  void updateImputations(MibrrData&, bool);
  // @param1: an initialized MibrrData object
  // @param2: store a posterior predictive sample (true) or update the
  //          imputations (false)?
  // @effect: update the imputations based on current values of member variables
  
  void doGibbsIteration(MibrrData&);
  // @param:  an initialized MibrrData object
  // @effect: run a single iteration of the Gibbs sampler


  //////////////////////// MCEM OPTIMIZATION FUNCTIONS /////////////////////////

  //// Beginning with Version 0.0.0.9001, the MCEM steps have been moved back to
  //// the R layer because installing NLopt has proven to be too much of a
  //// portability issue.
    
  //////////////////////// EXCEPTION HANDLING FUNCTIONS ////////////////////////

  
  void tauError(int) const;
  // @param:  the orignial error code
  // @effect: dispatch an appropriate error message to stderr

  void betaError(exception&) const;
  // @param:  orignial exception object
  // @effect: dispatch an appropriate error message to stderr
 
private:
  VectorXd _betas;
  ArrayXd  _taus;
  double   _sigma;
  double   _ridge;
  VectorXd _lambdas;
  VectorXd _l1Parms;
  VectorXd _l2Parms;
  VectorXd _betaMeans;
  MatrixXd _betaSam;
  ArrayXXd _tauSam;
  VectorXd _sigmaSam;
  MatrixXd _impSam;
  MatrixXd _ppSam;
  MatrixXd _lambdaSam;
  int      _targetIndex;
  int      _nDraws;
  int      _drawNum;
  int      _penType;
  bool     _verbose;
  bool     _storeGibbsSamples;
  bool     _doImp;
  bool     _savePpSams;
  bool     _fullBayes;
  bool     _useBetaMeans;
  bool     _finalRep;
};

#endif
