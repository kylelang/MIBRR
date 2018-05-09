// Title:    Header file for MibrrSamplers Class
// Author:   Kyle M. Lang (with some routines adapted from Josef Leydold's and
//           Robert E. Wheeler's C code)
// Created:  2017-NOV-23
// Modified: 2018-MAY-09
// Purpose:  These routines will generate pseudo-random variates for use in the
//           MIBRR routines.
// Note:     Some of these routines were adapted from the C code from other
//           authors (see notes inline). In such cases, all cleverness of the
//           original implementations is due to the original authors, and any
//           bugs induced by the port are my own responsibility. 

//--------------------- COPYRIGHT & LICENSING INFORMATION --------------------//
//  Copyright (C) 2018 Kyle M. Lang <k.m.lang@uvt.nl>                         //
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

#ifndef MIBRRSAMPLERS_H
#define MIBRRSAMPLERS_H

#include "MibrrDefs.h"
#include <random>
#include <Rmath.h>

class MibrrSamplers {

public:
  //////////////////////// CONSTRUCTORS / DESTRUCTOR ///////////////////////////

  MibrrSamplers();
  
  ~MibrrSamplers();

  ///////////////////////////////// MUTATORS ///////////////////////////////////
  
  void seedRng(const unsigned int seed);
  // @param: seed for the pseudo random number generator
  
  ///////////////////////////////// ACCESSORS //////////////////////////////////
  
  unsigned int getSeed() const;
  // @return: the current PRNG seed
  
  ///////////////////////////// SAMPLING FUNCTIONS /////////////////////////////

  double drawNorm();
  // @return: random standard normal variate
  
  VectorXd drawNorm(const int);
  // @param:  number of variates to draw
  // @return: vector of random standard normal variates
  
  double drawNorm(const double, const double);
  // @param1: mean parameter
  // @param2: SD parameter
  // @return: random normal variate

  double drawGamma(const double, const double);
  // @param1: shape parameter
  // @param2: rate parameter
  // @return: random Gamma variate

  double drawInvGamma(const double, const double);
  // @param1: shape parameter
  // @param2: scale parameter
  // @return: random Inverse Gamma variate

  double calcIncGamma(const double, const double, const bool);
  // @param1: shape parameter of the underlying gamma distribution (scale = 1.0)
  // @param2: threshold value cutting off the upper or lower tail (i.e., the
  //          underlying variate whose probability or [1 - probability] is being
  //          returned)
  // @param3: true = lower incomplete gamma, false = upper incomplete gamma
  // @return: value of the incomplete gamma function (i.e., the un-normalized 
  //          area under the upper or lower tail of the Gamma CDF).

  VectorXd drawMvn(const VectorXd&, const MatrixXd&);
  // @param1: mean vector
  // @param2: covariance matrix
  // @return: random multivariate normal variates
  
  //--------------------------------------------------------------------------//
  // The drawInvGauss() function was adapted from source code originally
  // implemented by Robert E. Wheeler (2001-MAR) in the R package SuppDists:
  //--------------------------------------------------------------------------//
  double drawInvGauss(const double, const double);
  // @param1: mean parameter (mu)
  // @param2: shape parameter (lambda)
  // @return: random variate from the inverse Gaussian distribution
  
  //--------------------------------------------------------------------------//
  // The drawGig() function was adapted from source code originally implemented
  // by Josef Leydold in the R package GIGrvg
  //--------------------------------------------------------------------------//
  double drawGig(const double, const double, const double);
  // @param1: the 'lambda' parameter of the GIG distribution
  // @param2: the 'chi' parameter of the GIG distribution
  // @param3: the 'psi' parameter of the GIG distribution
  // @return: a random GIG variate

protected:
  // PRNG objects:
  mt19937_64                        _gen;
  uniform_real_distribution<double> _unif; // defaults to standard uniform
  normal_distribution<double>       _norm; // defaults to standard normal

private:
  double _gigMode() const;
  // @return: analytic mode of the GIG distribution as currently parameterized
  
  double _gigRouNoShift();
  // @return: a GIG variate sampled via the approach that Leydold describes as:
  //          "Ratio-of-uniforms with shift by 'mode', alternative
  //          implementation"
     
  double _gigNewApproach();
  // @return: a GIG variate sampled via the approach that Leydold describes as:
  // "Ratio-of-uniforms without shift"
  
  double _gigRouShiftAlt();
  // @return: a GIG variate sampled via the approach that Leydold describes as:
  //          "New approach, constant hat in log-concave part."
  
  // member variables:
  double       _gigLam;
  double       _gigLam0;
  double       _omega;
  double       _alpha;
  double       _mode;
  unsigned int _seed;
  
};

#endif
