// Title:    Header file for GIG Sampling Routines
// Author:   Kyle M. Lang (slighly modified from source originally authored by
//           Josef Leydold)
// Created:  2017-NOV-23
// Modified: 2017-NOV-23
// Purpose:  These routines will generate pseudo-random variates from the
//           Generalized Inverse Gaussian distribution.
// Note:     These routines were adapted/copied from the C code underling Josef
//           Leydold's GIGrvg package. All cleverness of the origninal
//           implementation is due to Josef, and any bugs induced by the port
//           are my own responsibility. 

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


// define macros for GCC attributes:

#ifdef __GNUC__
#  define ATTRIBUTE__UNUSED  __attribute__ ((unused))
#else
#  define ATTRIBUTE__UNUSED
#endif

///////////////////////////// Function Prototypes ///////////////////////////////

//--Public Fuction-------------------------------------------------------------//

double* drawGig(int, double, double, double);
// @param1: sample size (positive integer, n)
// @param2: lambda parameter for distribution
// @param3: chi parameter for distribution
// @param4: psi parameter for distribution
// @return: random sample of size 'n'

//--Private Functions----------------------------------------------------------//

static double _gig_mode(double, double);
// @param1: lambda parameter for distribution
// @param2: omega parameter for distribution (i.e., sqrt(chi / psi))
// @return: mode of the distribution

static void _rgig_ROU_noshift(double *res,
			      int    n,
			      double lambda,
			      double lambda_old,
			      double omega,
			      double alpha);
// @param1: pointer to results vector
// @param2: sample size (positive integer)
// @param3: lambda parameter of the distribution
// @param4: previous lambda parameter
// @param5: omega parameter for the distribution (i.e., sqrt(chi / psi))
// @param6: alpha parameter for the distribution (i.e., sqrt(psi * chi))

static void _rgig_newapproach1(double *res,
			       int    n,
			       double lambda,
			       double lambda_old,
			       double omega,
			       double alpha);
// @param1: pointer to results vector
// @param2: sample size (positive integer)
// @param3: lambda parameter of the distribution
// @param4: previous lambda parameter
// @param5: omega parameter for the distribution (i.e., sqrt(chi / psi))
// @param6: alpha parameter for the distribution (i.e., sqrt(psi * chi))

static void _rgig_ROU_shift_alt(double *res,
				int    n,
				double lambda,
				double lambda_old,
				double omega,
				double alpha);
// @param1: pointer to results vector
// @param2: sample size (positive integer)
// @param3: lambda parameter of the distribution
// @param4: previous lambda parameter
// @param5: omega parameter for the distribution (i.e., sqrt(chi / psi))
// @param6: alpha parameter for the distribution (i.e., sqrt(psi * chi))

static double _unur_bessel_k_nuasympt(double x,
				      double nu,
				      int    islog,
				      int    expon_scaled);
//@param1: argument for Bessel K_nu()
//@param2: order of the Bessel function
//@param3: return logarithm of result when TRUE and raw result when FALSE
//@param4: return exp(-x) * K_nu(x) when TRUE and K_nu(x) when FALSE
