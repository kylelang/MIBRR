// Title:    Function definitions for the MibrrSamplers class
// Author:   Kyle M. Lang (with some routines adapted from Josef Leydold's and
//           Robert E. Wheeler's C code)
// Created:  2017-NOV-23
// Modified: 2019-JAN-24
// Purpose:  These routines will generate pseudo-random variates for use in the
//           MIBRR routines.
// Note:     Some of these routines were adapted from the C code from other
//           authors (see notes inline). In such cases, all cleverness of the
//           original implementations is due to the original authors, and any
//           bugs induced by the port are my own responsibility. 

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

#include "MibrrSamplers.h"

///////////////////////// CONSTRUCTORS / DESTRUCTOR ////////////////////////////

MibrrSamplers::MibrrSamplers() {}

MibrrSamplers::~MibrrSamplers() {}

////////////////////////////////// MUTATORS ////////////////////////////////////

void MibrrSamplers::seedRng(const unsigned int seed)
{
  // Store the new seed:
  _seed = seed;
  
  // Re-seed the PRNG:
  _gen.seed(seed);
}

////////////////////////////////// ACCESSORS ///////////////////////////////////

unsigned int MibrrSamplers::getSeed() const { return _seed;                    }

///////////////////////////// SAMPLING FUNCTIONS ///////////////////////////////

double MibrrSamplers::drawNorm(const double mean, const double sd)
{
  return mean + sd * _norm(_gen);
}

//----------------------------------------------------------------------------//

VectorXd MibrrSamplers::drawNorm(const int    n,
				 const double mean,
				 const double sd)
{
  VectorXd out(n);
  for(int i = 0; i < n; i++) out[i] = mean + sd * _norm(_gen);
  return out;
}

//----------------------------------------------------------------------------//

double MibrrSamplers::drawGamma(const double shape, const double rate)
{
  gamma_distribution<double> gam(shape, 1.0 / rate);
  return gam(_gen);
}

//----------------------------------------------------------------------------//

double MibrrSamplers::drawInvGamma(const double shape, const double scale)
{
  gamma_distribution<double> gam(shape, 1.0 / scale);
  return 1.0 / gam(_gen);
}

//----------------------------------------------------------------------------//

double MibrrSamplers::drawScaledInvChiSq(const double df, const double scale)
{
  return drawInvGamma(df / 2.0, (df * scale) / 2.0);
}

//----------------------------------------------------------------------------//

double MibrrSamplers::calcIncGamma(const double shape, 
				   const double cutVal,
				   const bool   lowerTail)
{
  double scale   = 1.0;
  int    lower   = (int)lowerTail;
  int    logTran = 0; // Don't want log transform
  
  return Rf_pgamma(cutVal, shape, scale, lower, logTran) * tgamma(shape);
}

//----------------------------------------------------------------------------//

double MibrrSamplers::drawInvGauss(const double mu, const double lambda)
{ 
  double b      = 0.5 * mu / lambda;
  double a      = mu * b;
  double c      = 4.0 * mu * lambda;
  double d      = pow(mu, 2);
  double outVal = 0.0;
  
  while(outVal <= 0.0) {
    double tmp = _norm(_gen);
    double v   = pow(tmp, 2); // Chi-Squared with df = 1
    
    if (mu <= 0.0)
      throw invalid_argument("The Inverse Gaussian's mean is non-positive.\n");
    
    else if (lambda <= 0.0) 
      throw invalid_argument("The Inverse Gaussian's scale is non-positive.\n");
    
    else {
      double u = _unif(_gen);

      //Find the smallest root:
      double x = mu + a * v - b * sqrt(c * v + d * pow(v, 2));

      // Choose x with probability = mean / (mean + x), else choose d/x:
      outVal = ( u < ( mu / (mu + x) ) ) ? x : d / x; 
    }  
  }// CLOSE while(outVal !> 0.0)
  return outVal;
}// END drawInvGauss()

//----------------------------------------------------------------------------//

VectorXd MibrrSamplers::drawMvn(const VectorXd &meanVec, const MatrixXd &covMat)
{
  int      nVars = meanVec.size();
  MatrixXd covCholesky;
  VectorXd normDraws(nVars);
  
  covCholesky = covMat.llt().matrixL();
  for(int i = 0; i < nVars; i++) normDraws[i] = _norm(_gen);
  
  return meanVec + (covCholesky * normDraws);
}// END drawMvn()

//----------------------------------------------------------------------------//

double MibrrSamplers::drawGig(const double lambda,
			      const double chi,
			      const double psi)
{ 
  // Store lambda parameter:
  _gigLam = lambda;
  
  // Define parameters for the standard distribution:
  _alpha  = sqrt(chi / psi);
  _omega  = sqrt(psi * chi);
  
  // check GIG parameters:
  bool check =
    isinf(_gigLam) || isinf(chi) || isinf(psi) || (chi < 0.0) || (psi < 0.0);
  
  if(check)
    throw invalid_argument("Invalid parameter(s) for GIG distribution.\n");
  
  if(chi < ZTOL) {
    // special cases which are basically Gamma and Inverse Gamma distribution
    if(_gigLam > 0.0) {
      gamma_distribution<double> gam(_gigLam, 2.0 / psi);
      return gam(_gen); 
    }
    else {
      gamma_distribution<double> gam(-_gigLam, 2.0 / psi);
      return (1.0 / gam(_gen)); 
    }
  }
  
  else if(psi < ZTOL) {
    // special cases which are basically Gamma and Inverse Gamma distribution
    if(_gigLam > 0.0) {
      gamma_distribution<double> gam(_gigLam, 2.0 / chi);
      return (1.0 / gam(_gen)); 
    }
    else {
      gamma_distribution<double> gam(-_gigLam, 2.0 / chi);
      return gam(_gen); 
    }
  }
  
  else {
    // Get a positive gigLam for computational purposes:
    _gigLam0 = _gigLam;
    if(_gigLam < 0.0) _gigLam = -_gigLam;

    // Store the mode of the distibution:
    _mode = _gigMode();
    
    // run full generator
    bool check1 = _gigLam > 2.0 || _omega > 3.0;
    bool check2 = _gigLam >= 1.0 - 2.25 * _omega * _omega || _omega > 0.2;
    
    if(check1)      return _gigRouShiftAlt();
    else if(check2) return _gigRouNoShift();
    else            return _gigNewApproach();
  }
} // END drawGig()


////////////////////////////// Private Functions ///////////////////////////////


double MibrrSamplers::_gigMode() const
{
  if(_gigLam >= 1.0) // return mode of fgig(x)
    return (sqrt((_gigLam - 1.0) * (_gigLam - 1.0) + _omega * _omega) +
	    (_gigLam - 1.0)) / _omega;
  
  else // 0 <= _gigLam < 1: use mode of f(1/x)
    return _omega / (sqrt((1.0 - _gigLam) * (1.0 - _gigLam) + _omega * _omega) +
		     (1.0 - _gigLam));
} // END _gigMode()

//----------------------------------------------------------------------------//

double MibrrSamplers::_gigRouNoShift()
//----------------------------------------------------------------------------//
// Type 1: Ratio-of-uniforms without shift.                                   //
//   Dagpunar (1988, Sect.4.6.2); Lehner (1989)                               //
//----------------------------------------------------------------------------//
{
  double nc;      // normalization constant
  double ym, um;  // location of maximum of x * sqrt(f(x)); umax of MBR
  double s, t;    // auxiliary variables
  double U, V, X; // random variables
 
  //--Setup-------------------------------------------------------------------//
  
  // shortcuts
  t = 0.5 * (_gigLam - 1.0);
  s = 0.25 * _omega;
  
  // normalization constant: c = log(sqrt(f(_mode)))
  nc = t * log(_mode) - s * (_mode + 1.0 / _mode);

  // location of maximum of x*sqrt(f(x)):
  // we need the positive root of _omega/2*y^2 - (_gigLam+1)*y - _omega/2 = 0
  ym = ((_gigLam + 1.0) +
	sqrt((_gigLam + 1.0) * (_gigLam + 1.0) + _omega * _omega)) / _omega;

  // boundaries of minimal bounding rectangle:
  // we us the "normalized" density f(x) / f(_mode).
  // hence
  //   upper boundary:      vmax = 1.                                 
  //   left hand boundary:  umin = 0.                             
  //   right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(_mode))
  um = exp(0.5 * (_gigLam + 1.0) * log(ym) - s * (ym + 1.0 / ym) - nc);

  //--Generate sample---------------------------------------------------------//
  
  // Run rejection algorithm:
  bool reject = true;
  while(reject) {
    U = um * _unif(_gen); // U(0, umax)
    V = _unif(_gen);      // U(0, vmax)
    X = U / V;
    
    reject = log(V) > (t * log(X) - s * (X + 1.0 / X) - nc); 
  };
    
  // return the accepted variate
  return (_gigLam0 < 0.0) ? (_alpha / X) : (_alpha * X);
} // END _gigRouNoShift()

//----------------------------------------------------------------------------//

double MibrrSamplers::_gigNewApproach()
//----------------------------------------------------------------------------//
// Type 4: New approach, constant hat in log-concave part.                    //
//   Case: 0 < _gigLam < 1, 0 < _omega < 1                                    //
//----------------------------------------------------------------------------//
{
  // parameters for hat function
  double A[3], Atot; // area below hat
  double k0;         // maximum of PDF
  double k1, k2;     // multiplicative constant

  double x0; // splitting point T-concave / T-convex
  double a;  // auxiliary variable

  double U, V, X; // random numbers
  double hx;      // hat at X
  
  //--Check arguments---------------------------------------------------------//

  if(_gigLam >= 1.0 || _omega > 1.0)
    throw invalid_argument("Invalid parameters for 'Type 4' GIG simulation.\n)");
  
  //--Setup-------------------------------------------------------------------//

  // splitting point
  x0 = _omega / (1.0 - _gigLam);

  // domain [0, x_0]
  k0   = exp((_gigLam - 1.0) * log(_mode) -
	     0.5 * _omega * (_mode + 1.0 / _mode)); // = f(_mode)
  A[0] = k0 * x0;

  // domain [x_0, Infinity]
  if(x0 >= 2.0 / _omega) {
    k1   = 0.0;
    A[1] = 0.0;
    k2   = pow(x0, _gigLam - 1.0);
    A[2] = k2 * 2.0 * exp(-_omega * x0 / 2.0) / _omega;
  }
  else {
    // domain [x_0, 2/_omega]
    k1   = exp(-_omega);
    A[1] = (_gigLam == 0.0) 
      ? k1 * log(2.0 / (_omega * _omega))
      : k1 / _gigLam * (pow(2.0 / _omega, _gigLam) - pow(x0, _gigLam));

    // domain [2/_omega, Infinity]
    k2   = pow(2.0 / _omega, _gigLam - 1.0);
    A[2] = k2 * 2.0 * exp(-1.0) / _omega;
  }

  // total area
  Atot = A[0] + A[1] + A[2];

  //--Generate sample---------------------------------------------------------//
  
  bool accept = false; // flag for rejection algorithm
  while(!accept) {
    // get uniform random number
    V = Atot * _unif(_gen);
    
    while(true) {
      // domain [0, x_0]
      if(V <= A[0]) {
	X  = x0 * V / A[0];
	hx = k0;
	break;
      }
      
      // domain [x_0, 2/_omega]
      V -= A[0];
      if(V <= A[1]) {
	if(_gigLam == 0.0) {
	  X  = _omega * exp(exp(_omega) * V);
	  hx = k1 / X;
	}
	else {
	  X  = pow(pow(x0, _gigLam) + (_gigLam / k1 * V), 1.0 / _gigLam);
	  hx = k1 * pow(X, _gigLam - 1.0);
	}
	break;
      }
      
      // domain [max(x0,2/_omega), Infinity]
      V  -= A[1];
      a  = (x0 > 2.0 / _omega) ? x0 : 2.0 / _omega;
      X  = -2.0 / _omega *
	log(exp(-_omega / 2.0 * a) - _omega / (2.0 * k2) * V);
      hx = k2 * exp(-_omega / 2.0 * X);
      break;
    }; // CLOSE while(true)
    
    // accept or reject?
    U      = _unif(_gen) * hx;
    accept = log(U) <= (_gigLam - 1.0) * log(X) - _omega / 2.0 * (X + 1.0 / X); 

    // Return the accepted variate
    if(accept) return (_gigLam0 < 0.0) ? (_alpha / X) : (_alpha * X);
  };
  throw runtime_error("Something has broken while sampling from the GIG; \
I've returned from a rejection sampling loop without a valid result.");
} // END _gigNewApproach()

//----------------------------------------------------------------------------//

double MibrrSamplers::_gigRouShiftAlt()
//----------------------------------------------------------------------------//
// Type 8:                                                                    //
// Ratio-of-uniforms with shift by 'mode', alternative implementation.        //
//   Dagpunar (1989); Lehner (1989)                                           //
//----------------------------------------------------------------------------//
{
  double nc;      // c = log(f(_mode)) normalization constant
  double s, t;    // auxiliary variables
  double U, V, X; // random variables

  double a, b, c; // coefficent of cubic
  double p, q;    // coefficents of depressed cubic
  double fi, fak; // auxiliary results for Cardano's rule

  double y1, y2; // roots of (1 / x) * sqrt(f((1 / x) + m))

  double uplus, uminus; // maximum and minimum of x * sqrt(f(x + m))

  //--Setup-------------------------------------------------------------------//

  // shortcuts
  t = 0.5 * (_gigLam - 1.);
  s = 0.25 * _omega;
  
  // normalization constant: c = log(sqrt(f(_mode)))
  nc = t * log(_mode) - s * (_mode + 1. / _mode);

  // location of minimum and maximum of (1 / x) * sqrt(f(1 / x + m)):

  // compute coeffients of cubic equation y^3 + a * y^2 + b * y + c = 0
  a = -(2.0 * (_gigLam + 1.0) / _omega + _mode);     // a < 0
  b = (2.0 * (_gigLam - 1.) * _mode / _omega - 1.0);
  c = _mode;

  // we need the roots in (0, _mode) and (_mode, inf)

  // substitute y = z - a / 3 for depressed cubic equation z^3 + p * z + q = 0
  p = b - a * a / 3.0;
  q = (2.0 * a * a * a) / 27.0 - (a * b) / 3.0 + c;

  // use Cardano's rule
  fi  = acos(-q / (2.0 * sqrt(-(p * p * p) / 27.0)));
  fak = 2.0 * sqrt(-p / 3.0);
  y1  = fak * cos(fi / 3.0) - a / 3.0;
  y2  = fak * cos(fi / 3.0 + 4.0 / 3.0 * M_PI) - a / 3.0;

  // boundaries of minimal bounding rectangle:
  // we us the "normalized" density f(x) / f(_mode).
  // hence:
  //   upper bound:      vmax = 1.
  //   left hand bound:  uminus = (y2-_mode) * sqrt(f(y2)) / sqrt(f(_mode))
  //   right hand bound: uplus = (y1-_mode) * sqrt(f(y1)) / sqrt(f(_mode))
  uplus  = (y1 - _mode) * exp(t * log(y1) - s * (y1 + 1.0 / y1) - nc);
  uminus = (y2 - _mode) * exp(t * log(y2) - s * (y2 + 1.0 / y2) - nc);

  //--Generate sample---------------------------------------------------------//

  // Run rejection sampling:
  bool reject = true;
  while(reject) {
    U = uminus + _unif(_gen) * (uplus - uminus); // U(u-, u+) 
    V = _unif(_gen);                             // U(0, vmax)
    X = U / V + _mode;
    
    reject = (X <= 0.0) || (log(V) > (t * log(X) - s * (X + 1.0 / X) - nc)); 
  };
  
  // Return the accepted variate
  return (_gigLam0 < 0.0) ? (_alpha / X) : (_alpha * X);
} // END _gigRouShiftAlt()
