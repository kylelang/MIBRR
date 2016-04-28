// Title:    Header file to hold my global parameter definitions
// Author:   Kyle M. Lang
// Created:  2014-AUG-25
// Modified: 2015-FEB-28

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

#ifndef MYPARAMS_H
#define MYPARAMS_H

#include <RcppEigen.h>
#include <iostream>
#include <cmath>
#include <string>
//#include <eigen3/Eigen/Dense>
//#include <boost/random.hpp>
//#include <boost/math/distributions.hpp>
//#include <R.h>
#include<Rmath.h>
#include <nlopt.hpp>

#ifndef MACHINE_PRECISION
#define MACHINE_PRECISION std::numeric_limits<double>::epsilon()
#endif

//#ifndef MACHINE_NAN
//#define MACHINE_NAN std::numeric_limits<double>::quiet_NaN()
//#endif

//typedef boost::random::gamma_distribution <double> br_gamma_distribution;
//typedef boost::math::gamma_distribution <double> bm_gamma_distribution;
//typedef boost::random::mt19937 MyGenerator;

typedef Eigen::Array <bool, Eigen::Dynamic, 1> ArrayXb;
typedef Eigen::Array <bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXXb;
typedef Eigen::Matrix <double, 5, 1> Vector5d;
//typedef Eigen::Map <const Eigen::MatrixXd> DataMap;


#endif
