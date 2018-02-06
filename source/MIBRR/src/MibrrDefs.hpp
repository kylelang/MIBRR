// Title:    Header file to hold my global parameter definitions
// Author:   Kyle M. Lang
// Created:  2014-AUG-25
// Modified: 2017-NOV-17

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

#ifndef MIBRRDEFS_H
#define MIBRRDEFS_H

#include <RcppEigen.h>
#include <iostream>
#include <cmath>
#include <string>
#include <Rmath.h>
#include <algorithm>

#ifndef MACHINE_PRECISION
#define MACHINE_PRECISION std::numeric_limits<double>::epsilon()
#endif

typedef Eigen::Array <bool, Eigen::Dynamic, 1>              ArrayXb;
typedef Eigen::Array <bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXXb;
typedef Eigen::Matrix <double, 5, 1>                        Vector5d;
typedef Rcpp::List                                          RList;
//typedef Eigen::Map <const Eigen::MatrixXd> DataMap;

#endif
