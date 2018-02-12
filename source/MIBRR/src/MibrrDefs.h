// Title:    Header file to hold my global parameter definitions
// Author:   Kyle M. Lang
// Created:  2014-AUG-25
// Modified: 2018-FEB-12

//--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------//
//  Copyright (C) 2018 Kyle M. Lang <k.m.lang@uvt.nl>                          //  
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

#ifndef MIBRRDEFS_H
#define MIBRRDEFS_H

#include <RcppEigen.h>
#include <iostream>
#include <cmath>
#include <string>

#ifndef MACHINE_PRECISION
#define MACHINE_PRECISION std::numeric_limits<double>::epsilon()
#endif

#ifndef ZTOL
#define ZTOL MACHINE_PRECISION * 10
#endif

#ifndef M_LNPI
#define M_LNPI 1.14472988584940017414342735135 // ln(pi)
#endif

typedef Eigen::Array <bool, Eigen::Dynamic, 1             > ArrayXb;
typedef Eigen::Array <bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXXb;

typedef Eigen::Matrix <double,       5,              1> Vector5d;
//typedef Eigen::Matrix <unsigned int, Eigen::Dynamic, 1> VectorXui;

typedef Rcpp::List RList;
//typedef Eigen::Map <const Eigen::MatrixXd> DataMap;

using namespace Eigen;
using namespace std;

#endif
