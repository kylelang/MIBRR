// Title:    Header file for the MyGibbs Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2015-FEB-25
// Purpose:  This class contains the Gibbs sampling-related functions
//           for the MIBRR package.

//--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------//
//    Copyright (C) 2015 Kyle M. Lang <kylelang@ku.edu>                        //  
//                                                                             //
//    This program is free software: you can redistribute it and/or modify     //
//    it under the terms of the GNU General Public License as published by     //
//    the Free Software Foundation, either version 3 of the License, or        //
//    (at your option) any later version.                                      //
//                                                                             //
//    This program is distributed in the hope that it will be useful,          //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of           //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            //
//    GNU General Public License for more details.                             //
//                                                                             //
//    You should have received a copy of the GNU General Public License        //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.    //
//-----------------------------------------------------------------------------//

#ifndef MYGIBBS_H
#define MYGIBBS_H


#include "MyParams.hpp"
#include "MyData.hpp"
#include <nlopt.hpp>

using namespace std;
using namespace Eigen;

class MyGibbs {

 public:

  ///// CONSTRUCTORS / DESTRUCTOR /////

  MyGibbs();
  
  ~MyGibbs();


  ///// ACCESSORS /////

  VectorXd getBetas() const;
  // @return a vector of regression coefficients 
  // (and the intercept term) for the elastic net model

  ArrayXd getTaus() const;
  // @return an array of auxiliary penalty hyperparameters 
  // for the elastic net model

  double getSigma() const;
  // @return the residual variance of the elastic net model

  MatrixXd getBetaSam() const; 

  ArrayXXd getTauSam() const; 

  VectorXd getSigmaSam() const;

  MatrixXd getImpSam() const;
  
  MatrixXd getLambdaHistory() const;
  
  VectorXd getLambdas() const;
  // @return the current values of the penalty parameters
  
  double getLambdas(int) const;
  // @param a number indicating which lambda to return (either 1 or 2)
  // @return the current value of either lambda1 and lambda2

  int getNDraws() const;
  // @return the current number of retained Gibbs sampling draws
  

  ///// MUTATORS /////
  
  void setBetas(VectorXd&);
  // @param a new vector of coefficients for the elastic net model
  
  void setTaus(ArrayXd&);
  // @param a new array of hyperparameters for the elastic net model
  
  void setSigma(double);
  // @param a new residual variance of the elastic net model
  
  void setLambdas(VectorXd&);
  // @param a new value for Lambda

  void setLambdas(double, double);
  // @param1 a new value for lambda1
  // @param2 a new value for lambda2

  void setLambdas(double);
  // @param1 a new value for the LASSO lambda

  void setLambdas();

  void setNDraws(int);
  //@param a new value for the retained Gibbs sampling draws

  void setNEmIters(int);
  //@param a new value for the number of MCEM iterations

  void setTargetIndex(int);
  // @param a new value for the target variable's column index

  void setVerbosity(bool, bool);
  // @param1 logical switch indicating if verbose iteration history
  // should be reported (true) or not (false)
  // @param2 logical switch indicating if verbose errors
  // should be reported (true) or not (false)
  
  void doMiben();
  // Call to set the imputation model to the Bayesian elastic net 
  
  void doMibl();
  // Call to set the imputation model to the Bayesian LASSO 
  
  void startGibbsSampling(MyData&);
  // Start storing the parameters' Gibbs sampled values
 
  void stopGibbsSampling();

  void setupOptimizer(int,
		      int,
		      double);
  // Parameterize the MIBEN EM optimization
  
  void setupOptimizer(int, int);
  // Parameterize the MIBL EM optimization

  void startParameters(VectorXd&,
		       ArrayXd&,
		       double,
		       double,
		       double);
  
  void startParameters(VectorXd&,
		       ArrayXd&,
		       double,
		       VectorXd&);

  void restartParameters(MyData&);
  

  ///// PARAMETER UPDATE FUNCTIONS /////

  void updateTaus(MyData&);
  // Update _taus (for the Bayesian Elastic Net) based on 
  // current values of other member variables

  void updateBetas(MyData&);
  // Update _betas based on current values of other member variables

  void updateSigma(MyData&);
  // Update _sigma (for the Bayesian Elastic Net) based on 
  // current values of other member variables
  
  void updateImputations(MyData&);
  // Update the imputations based on current values of member variables
  
  void doGibbsIteration(MyData&);
  // Run a single iteration of the Gibbs sampler


  ///// MCEM OPTIMIZATION FUNCTIONS //////

  double eNetLambdaObjective(const std::vector<double>&,
			     std::vector<double>&,
			     void*);
  // The objective function (with gradient embedded) for the 
  // Elastic Net penalty parameters

  void optimizeMibenLambdas(bool); 
  // Numerically optimize the MIBEN penalty parameters

  void updateLambdas();
  // Update the ENET or LASSO penalty parameters using
  // marginal, numerical optimization.

 private:
  
  VectorXd _betas;
  ArrayXd _taus;
  double _sigma;
  MatrixXd _betaSam;
  ArrayXXd _tauSam;
  VectorXd _sigmaSam;
  MatrixXd _impSam;
  MatrixXd _lambdaHistory;
  VectorXd _lambdas;
  double _emConvTol;
  int _nDraws;
  int _nEmIters;
  int _targetIndex;
  int _drawNum;
  int _emIterNum;
  int _optIterCount;
  int _lambdaWindow;
  bool _verboseErrors;
  bool _verboseIters;
  bool _useElasticNet;
  bool _storeGibbsSamples;
};


#endif
