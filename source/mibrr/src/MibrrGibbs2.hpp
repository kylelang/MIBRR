// Title:    Header file for the MibrrGibbs2 Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2016-MAY-14
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

#ifndef MIBRRGIBBS2_H
#define MIBRRGIBBS2_H

#include "MibrrGibbs.hpp"

using namespace std;
using namespace Eigen;

class MibrrGibbs2: public MibrrGibbs {
  
public:
  /////////////////////// CONSTRUCTORS / DESTRUCTOR /////////////////////////////
  
  MibrrGibbs2();
  
  ~MibrrGibbs2();

  ////////////////////////////////// MUTATOR ////////////////////////////////////
  
  void setSimpleIntercept(bool);
  //@param: a value for the simple intercept flag
  
  ///////////////////////// PARAMETER UPDATE FUNCTIONS //////////////////////////

  void updateTaus2(MibrrData2&);
  // @param:  an initialized MibrrData object
  // @effect: update _taus based on current values of other member variables

  void updateBetas2(MibrrData2&);
  // @param:  an initialized MibrrData object
  // @effect: update _betas based on current values of other member variables

  void updateSigma2(MibrrData2&);
  // @param:  an initialized MibrrData object
  // @effect: update _sigma based on current values of other member variables
  
  void updateImputations2(MibrrData2&);
  // @param:  an initialized MibrrData object
  // @effect: update the imputations based on current values of member variables
  
  void doGibbsIteration2(MibrrData2&);
  // @param:  an initialized MibrrData object
  // @effect: run a single iteration of the Gibbs sampler

  //////////////////////// MCEM OPTIMIZATION FUNCTIONS //////////////////////////

  double eNetLambdaObjective2(const std::vector<double>&,
			      std::vector<double>&,
			      void*);
  // @param1: starting values for Lambda
  // @param2: starting values for Lambda's gradient
  // @param3: pointer to additional data
  // @return: the BEN penalty parameters' evaluated EM objective function 
  
  void optimizeMibenLambdas2(bool);
  // @param:  are pre-optimizing (true) or fully optimizing (false) Lambda?
  // @effect: numerically optimize the MIBEN penalty parameters
  
  void updateLambdas2();
  // @effect: update the BEN or LASSO penalty parameters using marginal,
  //          numerical/deterministic optimization.

private:
  bool _simpleIntercept;
};

#endif
