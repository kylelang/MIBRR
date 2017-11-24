// Title:    Header file for MibrrSamplers Class
// Author:   Kyle M. Lang (with some routines adapted from Josef Leydold's and
//           Robert E. Wheeler's C code)
// Created:  2017-NOV-23
// Modified: 2017-NOV-24
// Purpose:  These routines will generate pseudo-random variates for use in the
//           MIBRR routines.
// Note:     Some of these routines were adapted from the C code from other
//           authors (see notes inline). In such cases, all cleverness of the
//           original implementations is due to the original authors, and any
//           bugs induced by the port are my own responsibility. 

//--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------//
//  Copyright (C) 2017 Kyle M. Lang <k.m.lang@uvt.nl>                          //  
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

#ifndef MIBRRSAMPLERS_H
#define MIBRRSAMPLERS_H

#include "/home/kmlang/sg/software/eigen-3.3.4/Eigen/Eigen"

//#include <RcppEigen.h>
#include <Rmath.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <random>

#ifndef ZTOL
#define ZTOL std::numeric_limits<double>::epsilon() * 10
#endif

#ifndef M_LNPI
#define M_LNPI 1.14472988584940017414342735135 // ln(pi)
#endif

using namespace std;
using namespace Eigen;

class MibrrSamplers {

public:
  //////////////////////// CONSTRUCTORS / DESTRUCTOR ////////////////////////////
  
  MibrrSamplers(const unsigned int);
  // @param1: a seed for the pseudo random number generator

  MibrrSamplers();
  
  ~MibrrSamplers();

  ///////////////////////////// SAMPLING FUNCTIONS //////////////////////////////

  double drawInvGamma(const double, const double);
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

  VectorXd drawMvn(const VectorXd&, const MatrixXd&);
  // @param1: mean vector
  // @param2: covariance matrix
  // @return: random multivariate normal variates
  
  //---------------------------------------------------------------------------//
  // The drawInvGauss() function was adapted from source code originally
  // implemented by Robert E. Wheeler (2001-MAR) in the R package SuppDists:
  //---------------------------------------------------------------------------//
  double drawInvGauss(const double, const double);
  // @param1: mean parameter (mu)
  // @param2: shape parameter (lambda)
  // @return: random variate from the inverse Gaussian distribution
  
  //---------------------------------------------------------------------------//
  // The drawGig() function was adapted from source code originally implemented
  // by Josef Leydold in the R package GIGrvg
  //---------------------------------------------------------------------------//
  double drawGig(const double, const double, const double);
  // @param1: the 'lambda' parameter of the GIG distribution
  // @param2: the 'chi' parameter of the GIG distribution
  // @param3: the 'psi' parameter of the GIG distribution
  // @return: a random GIG variate
  
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
  double _gigLam;
  double _gigLam0;
  double _omega;
  double _alpha;
  double _mode;
  
  // PRNG objects:
  mt19937_64                        _gen;
  uniform_real_distribution<double> _unif;
  normal_distribution<double>       _norm;
};

#endif
