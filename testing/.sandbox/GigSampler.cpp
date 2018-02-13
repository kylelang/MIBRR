// Title:    Function definitions for the GigSampler class
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

#include "GigSampler.hpp"

///////////////////////// CONSTRUCTORS / DESTRUCTOR /////////////////////////////

GigSampler::GigSampler(const double lambda, const double chi, const double psi)
{
  _lambda = lambda;
  _chi    = chi;
  _psi    = psi;

  // Define parameters for the standard distribution:
  _alpha  = sqrt(_chi / _psi);
  _omega  = sqrt(_psi * _chi);
   
  // Initialize the PRNG:
  mt19937_64 _gen();

  // Initialize a standard uniform random number generator:
  uniform_real_distribution<double> _unif(0.0, 1.0);
}

GigSampler::GigSampler() {}

GigSampler::~GigSampler() {}

///////////////////////////////// ACCESSORS /////////////////////////////////////

double         GigSampler::getLambda() const { return _lambda;                  }
double         GigSampler::getChi   () const { return _chi;                     }
double         GigSampler::getPsi   () const { return _psi;                     }
double         GigSampler::getMode  () const { return _mode;                    }
vector<double> GigSampler::getOutput() const { return _output;                  }

///////////////////////////////// MUTATORS //////////////////////////////////////

void GigSampler::setLambda(const double lambda)         { _lambda = lambda;     }
void GigSampler::setChi   (const double chi)            { _chi = chi;           }
void GigSampler::setPsi   (const double psi)            { _psi = psi;           }
void GigSampler::setMode  (const double mode)           { _mode = mode;         }
void GigSampler::setOutput(const vector<double> output) { _output = output;     }

void GigSampler::setupGig(const double lambda,
			  const double chi,
			  const double psi)
{
  _lambda = lambda;
  _chi    = chi;
  _psi    = psi;

  // Define parameters for the standard distribution:
  _alpha  = sqrt(_chi / _psi);
  _omega  = sqrt(_psi * _chi);
}

///////////////////////////// SAMPLING FUNCTION /////////////////////////////////

vector<double> GigSampler::drawGig(const int n)
{  
  // initialize the output vector:
  _output = vector<double>(n);
  
  // check sample size
  if(n <= 0) cerr << "Sample size 'n' must be positive integer." << endl;
  
  // check GIG parameters:
  bool check =
    isinf(_lambda) || isinf(_chi) || isinf(_psi) || (_chi < 0.0) || (_psi < 0.0);
  
  if(check) cerr << "Invalid parameter(s) for GIG distribution." << endl;
  
  if(_chi < ZTOL) {
    // special cases which are basically Gamma and Inverse Gamma distribution
    if(_lambda > 0.0) {
      gamma_distribution<double> _gamma(_lambda, 2.0 / _psi);
      for(int i = 0; i < n; i++) _output[i] = _gamma(_gen); 
    }
    else {
      gamma_distribution<double> _gamma(-_lambda, 2.0 / _psi);
      for(int i = 0; i < n; i++) _output[i] = 1.0 / _gamma(_gen); 
    }
  }
  
  else if(_psi < ZTOL) {
    // special cases which are basically Gamma and Inverse Gamma distribution
    if(_lambda > 0.0) {
      gamma_distribution<double> _gamma(_lambda, 2.0 / _chi);
      for(int i = 0; i < n; i++) _output[i] = 1.0 / _gamma(_gen); 
    }
    else {
      gamma_distribution<double> _gamma(-_lambda, 2.0 / _chi);
      for(int i = 0; i < n; i++) _output[i] = _gamma(_gen); 
    }
  }
  
  else {
    // Get a positive lambda for computational purposes:
    _lambda0 = _lambda;
    if(_lambda < 0.0) _lambda = -_lambda;

    // Store the mode of the distibution:
    _mode = _gigMode();
    
    // run full generator
    bool check1 = _lambda > 2.0 || _omega > 3.0;
    bool check2 = _lambda >= 1.0 - 2.25 * _omega * _omega || _omega > 0.2;
    
    if(check1)      _rgigRouShiftAlt();
    else if(check2) _rgigRouNoShift();
    else            _rgigNewApproach();
  }
  
  return _output;
} // END do_rgig()

////////////////////////////// Private Functions ////////////////////////////////

double GigSampler::_gigMode() const
{
  if(_lambda >= 1.0) // return mode of fgig(x)
    return (sqrt((_lambda - 1.0) * (_lambda - 1.0) + _omega * _omega) +
	    (_lambda - 1.0)) / _omega;
  
  else // 0 <= _lambda < 1: use mode of f(1/x)
    return _omega / (sqrt((1.0 - _lambda) * (1.0 - _lambda) + _omega * _omega) +
		     (1.0 - _lambda));
} // END _gigMode()


void GigSampler::_rgigRouNoShift()
//-----------------------------------------------------------------------------//
// Type 1: Ratio-of-uniforms without shift.                                    //
//   Dagpunar (1988, Sect.4.6.2); Lehner (1989)                                //
//-----------------------------------------------------------------------------//
{
  double nc;      // normalization constant
  double ym, um;  // location of maximum of x * sqrt(f(x)); umax of MBR
  double s, t;    // auxiliary variables
  double U, V, X; // random variables
 
  //--Setup--------------------------------------------------------------------//

  
  // shortcuts
  t = 0.5 * (_lambda - 1.0);
  s = 0.25 * _omega;
  
  // normalization constant: c = log(sqrt(f(_mode)))
  nc = t * log(_mode) - s * (_mode + 1.0 / _mode);

  // location of maximum of x*sqrt(f(x)):
  // we need the positive root of _omega/2*y^2 - (_lambda+1)*y - _omega/2 = 0
  ym = ((_lambda + 1.0) +
	sqrt((_lambda + 1.0) * (_lambda + 1.0) + _omega * _omega)) / _omega;

  // boundaries of minimal bounding rectangle:
  // we us the "normalized" density f(x) / f(_mode).
  // hence
  //   upper boundary:      vmax = 1.                                 
  //   left hand boundary:  umin = 0.                             
  //   right hand boundary: umax = ym * sqrt(f(ym)) / sqrt(f(_mode))
  um = exp(0.5 * (_lambda + 1.0) * log(ym) - s * (ym + 1.0 / ym) - nc);

  //--Generate sample----------------------------------------------------------//

  for(int i = 0; i < _output.size(); i++) {
    // Run rejection algorithm:
    bool reject = true;
    while(reject) {
      U = um * _unif(_gen); // U(0, umax)
      V = _unif(_gen);      // U(0, vmax)
      X = U / V;
      
      reject = log(V) > (t * log(X) - s * (X + 1.0 / X) - nc); 
    };
    
    // store the accepted variate
    _output[i] = (_lambda0 < 0.0) ? (_alpha / X) : (_alpha * X);
  }
  return;
} // END _rgigRouNoShift()


void GigSampler::_rgigNewApproach()
//-----------------------------------------------------------------------------//
// Type 4: New approach, constant hat in log-concave part.                     //
//   Case: 0 < _lambda < 1, 0 < _omega < 1                                     //
//-----------------------------------------------------------------------------//
{
  // parameters for hat function
  double A[3], Atot; // area below hat
  double k0;         // maximum of PDF
  double k1, k2;     // multiplicative constant

  double x0; // splitting point T-concave / T-convex
  double a;  // auxiliary variable

  double U, V, X; // random numbers
  double hx;      // hat at X
  
  //--Check arguments----------------------------------------------------------//

  if(_lambda >= 1.0 || _omega > 1.0) cerr << "invalid parameters" << endl;

  //--Setup--------------------------------------------------------------------//

  // splitting point
  x0 = _omega / (1.0 - _lambda);

  // domain [0, x_0]
  k0   = exp((_lambda - 1.0) * log(_mode) -
	     0.5 * _omega * (_mode + 1.0 / _mode)); // = f(_mode)
  A[0] = k0 * x0;

  // domain [x_0, Infinity]
  if(x0 >= 2.0 / _omega) {
    k1   = 0.0;
    A[1] = 0.0;
    k2   = pow(x0, _lambda - 1.0);
    A[2] = k2 * 2.0 * exp(-_omega * x0 / 2.0) / _omega;
  }
  else {
    // domain [x_0, 2/_omega]
    k1   = exp(-_omega);
    A[1] = (_lambda == 0.0) 
      ? k1 * log(2.0 / (_omega * _omega))
      : k1 / _lambda * (pow(2.0 / _omega, _lambda) - pow(x0, _lambda));

    // domain [2/_omega, Infinity]
    k2   = pow(2.0 / _omega, _lambda - 1.0);
    A[2] = k2 * 2.0 * exp(-1.0) / _omega;
  }

  // total area
  Atot = A[0] + A[1] + A[2];

  //--Generate sample----------------------------------------------------------//
  
  for(int i = 0; i < _output.size(); i++) {
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
	  if(_lambda == 0.0) {
	    X  = _omega * exp(exp(_omega) * V);
	    hx = k1 / X;
	  }
	  else {
	    X  = pow(pow(x0, _lambda) + (_lambda / k1 * V), 1.0 / _lambda);
	    hx = k1 * pow(X, _lambda - 1.0);
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
      accept = log(U) <= (_lambda - 1.0) * log(X) - _omega / 2.0 * (X + 1.0 / X); 
      
      if(accept)
	_output[i] = (_lambda0 < 0.0) ? (_alpha / X) : (_alpha * X);
    }; // CLOSE while(!accept)
  }
  return;
} // END _rgigNewApproach()


void GigSampler::_rgigRouShiftAlt()
//-----------------------------------------------------------------------------//
// Type 8: Ratio-of-uniforms with shift by 'mode', alternative implementation. //
//   Dagpunar (1989); Lehner (1989)                                            //
//-----------------------------------------------------------------------------//
{
  double nc;      // c = log(f(_mode)) normalization constant
  double s, t;    // auxiliary variables
  double U, V, X; // random variables

  double a, b, c; // coefficent of cubic
  double p, q;    // coefficents of depressed cubic
  double fi, fak; // auxiliary results for Cardano's rule

  double y1, y2; // roots of (1 / x) * sqrt(f((1 / x) + m))

  double uplus, uminus; // maximum and minimum of x * sqrt(f(x + m))

  //--Setup--------------------------------------------------------------------//

  // shortcuts
  t = 0.5 * (_lambda - 1.);
  s = 0.25 * _omega;
  
  // normalization constant: c = log(sqrt(f(_mode)))
  nc = t * log(_mode) - s * (_mode + 1. / _mode);

  // location of minimum and maximum of (1 / x) * sqrt(f(1 / x + m)):

  // compute coeffients of cubic equation y^3 + a * y^2 + b * y + c = 0
  a = -(2.0 * (_lambda + 1.0) / _omega + _mode);     // a < 0
  b = (2.0 * (_lambda - 1.) * _mode / _omega - 1.0);
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

  //--Generate sample----------------------------------------------------------//

  for(int i = 0; i < _output.size(); i++) {
    // Run rejection sampling:
    bool reject = true;
    while(reject) {
      U = uminus + _unif(_gen) * (uplus - uminus); // U(u-, u+) 
      V = _unif(_gen);                             // U(0, vmax)
      X = U / V + _mode;
      
      reject = (X <= 0.0) || (log(V) > (t * log(X) - s * (X + 1.0 / X) - nc)); 
    };
    
    // store accepted variate
    _output[i] = (_lambda0 < 0.0) ? (_alpha / X) : (_alpha * X);
  }
  return;
} // END _rgigRouShiftAlt()


double GigSampler::_asymptBesselK(double x,
				  double nu,
				  int    isLog,
				  int    expScale)
//-----------------------------------------------------------------------------//
// Asymptotic expansion of Bessel K_nu(x) function for large nu and x          //
//                                                                             //
//   Reference: Abramowitz & Stegun, p.378, __9.7.8.__                         //
//                                                                             //
//     K_nu(nu * z) ~ sqrt(pi/(2*nu)) * exp(-nu*eta)/(1+z^2)^(1/4)             //
//                       * {1 - u_1(t)/nu + u_2(t)/nu^2 - ... }                //
//                                                                             //
//     where   t = 1 / sqrt(1 + z^2),                                          //
//           eta = sqrt(1 + z^2) + log(z / (1 + sqrt(1+z^2)))                  //
//                                                                             //
//     and u_k(t)  from  p.366  __ 9.3.9 __                                    //
//                                                                             //
//     u0(t) = 1                                                               //
//     u1(t) = (3*t - 5*t^3)/24                                                //
//     u2(t) = (81*t^2 - 462*t^4 + 385*t^6)/1152                               //
//     ...                                                                     //
//                                                                             //
//     with recursion  9.3.10    for  k = 0, 1, .... :                         //
//                                                                             //
//     u_{k+1}(t) = t^2/2 * (1 - t^2) * u'_k(t) +                              //
//                1/8  \int_0^t (1 - 5*s^2)* u_k(s) ds                         //
//-----------------------------------------------------------------------------//
//                                                                             //
// Original implementation in R code (R package "Bessel" v. 0.5-3) by          //
//   Martin Maechler, Date: 23 Nov 2009, 13:39                                 //
//                                                                             //
// Translated into C code by Kemal Dingic, Oct. 2011.                          //
//                                                                             //
// Modified by Josef Leydold on Tue Nov  1 13:22:09 CET 2011                   //
//                                                                             //
// Translated to C++ and modified, for inclusion in MIBRR, by Kyle M. Lang     //
//   on 2017 NOV 24                                                            //
//-----------------------------------------------------------------------------//
{
  double z;                     // rescaled argument for K_nu()
  double sz, t, t2, eta;        // auxiliary variables
  double d, u1t, u2t, u3t, u4t; // (auxiliary) results for Debye polynomials
  double res;                   // value of log(K_nu(x)) [= result]
  
  // rescale: we comute K_nu(z * nu)
  z = x / nu;

  // auxiliary variables:
  sz = hypot(1, z); // = sqrt(1 + z^2)
  t  = 1.0 / sz;
  t2 = t * t;

  eta = (expScale) ? (1.0 / (z + sz)) : sz;
  eta += log(z) - log1p(sz); // = log(z / (1 + sz))

  // evaluate Debye polynomials u_j(t):
  u1t = (t * (3.0 - 5.0 * t2)) / 24.0;
  u2t = t2 * (81.0 + t2 * (-462.0 + t2 * 385.0)) / 1152.0;
  u3t = t * t2 * (30375.0
		  + t2 * (-369603.0
			  + t2 * (765765.0
				  - t2 * 425425.0))) / 414720.0;
  u4t = t2 * t2 * (4465125.0 
		   + t2 * (-94121676.0
			   + t2 * (349922430.0 
				   + t2 * (-446185740.0 
					   + t2 * 185910725.0)))) / 39813120.0;
  d = (-u1t + (u2t + (-u3t + u4t / nu) / nu) / nu) / nu;

  // log(K_nu(x)):
  res = log(1.0 + d) - nu * eta - 0.5 * (log(2.0 * nu * sz) - M_LNPI);

  return (isLog ? res : exp(res));
} // END _asymptBesselK()
