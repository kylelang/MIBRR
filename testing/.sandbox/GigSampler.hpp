// Title:    Header file for GigSampler Class
// Author:   Kyle M. Lang (adapted from Josef Leydold's C code)
// Created:  2017-NOV-23
// Modified: 2017-NOV-24
// Purpose:  These routines will generate pseudo-random variates from the
//           Generalized Inverse Gaussian distribution.
// Note:     These routines were adapted from the C code underling Josef
//           Leydold's GIGrvg package. All cleverness of the origninal
//           implementation is due to Josef (and/or ealier authors), and any bugs
//           induced by the port are my own responsibility. 

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

#ifndef GIGSAMPLER_H
#define GIGSAMPLER_H

//#include <R.h>
//#include <RcppEigen.h>
//#include <Rmath.h>
#include <cmath>
#include <vector>
//#include <algorithm>
#include <iostream>
#include <random>

#ifndef ZTOL
#define ZTOL std::numeric_limits<double>::epsilon() * 10
#endif

#ifndef M_LNPI
#define M_LNPI 1.14472988584940017414342735135 // ln(pi)
#endif

using namespace std;

class GigSampler {

public:
  //////////////////////// CONSTRUCTORS / DESTRUCTOR ////////////////////////////
  
  GigSampler(double, double, double);
  // @param1: lambda parameter for distribution
  // @param2: chi parameter for distribution
  // @param3: psi parameter for distribution

  GigSampler();
  
  ~GigSampler();

  //////////////////////////////// ACCESSORS ////////////////////////////////////

  double getLambda() const;
  double getChi() const;
  double getPsi() const;
  double getMode() const;
  vector<double> getOutput() const;
  
  //////////////////////////////// MUTATORS /////////////////////////////////////

  void setLambda(const double);
  void setChi(const double);
  void setPsi(const double);
  void setMode(const double);
  void setOutput(const vector<double>);
  
  void setupGig(const double, const double, const double);
  // @param1: lambda parameter for distribution
  // @param2: chi parameter for distribution
  // @param3: psi parameter for distribution
  // @effect: parameterize the GIG sampler

  ///////////////////////////// SAMPLING FUNCTION ///////////////////////////////
  
  vector<double> drawGig(const int);
  // @param: sample size (positive integer, n)
  // @return: a vector of 'n' random GIG variates
  
private:
  double _gigMode() const;
  // @return: mode of the distribution
  
  double _asymptBesselK(double, double, int, int);
  //@param1: argument for Bessel K_nu()
  //@param2: order of the Bessel function
  //@param3: return logarithm of result when TRUE and raw result when FALSE
  //@param4: return exp(-x) * K_nu(x) when TRUE and K_nu(x) when FALSE
  //@return: the asymptotic value of K_nu(x)
  
  void _rgigRouNoShift();
  // Ratio-of-uniforms with shift by 'mode', alternative implementation
     
  void _rgigNewApproach();
  // Ratio-of-uniforms without shift

  void _rgigRouShiftAlt();
  // New approach, constant hat in log-concave part.
  
  // member variables:
  double _lambda;
  double _lambda0;
  double _chi;
  double _psi;
  double _omega;
  double _alpha;
  double _mode;

  vector<double> _output;
 
  mt19937_64                        _gen;
  uniform_real_distribution<double> _unif;
};

#endif
