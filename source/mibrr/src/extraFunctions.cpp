// Title:    Additional C++ Function to Export in MIBRR
// Author:   Kyle M. Lang
// Created:  2014-AUG-20
// Modified: 2017-OCT-27

//--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------//
//  Copyright (C) 2017 Kyle M. Lang <k.m.lang@uvt.nl>                          //  
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

#include <RcppEigen.h>
#include "MibrrData.hpp"
#include "MibrrGibbs.hpp"

// [[Rcpp::export]]
Eigen::VectorXd drawInvGamma(int n, double shape, double scale)
{
  MibrrGibbs mibrrGibbs;
  VectorXd   out(n);
  for(int i = 0; i < n; i++) out[i] = mibrrGibbs.drawInvGamma(shape, scale);
  return out;
}

// [[Rcpp::export]]
Eigen::MatrixXd drawMVN(int n, Eigen::VectorXd meanVec, Eigen::MatrixXd covMat)
{
  MibrrData mibrrData;
  int       v = meanVec.size();
  MatrixXd  out(n, v);
  for(int i = 0; i < n; i++) out.row(i) = mibrrData.drawMVN(meanVec, covMat);
  return out;
}

// [[Rcpp::export]]
double calcIncGamma(double shape, double cutVal, bool lowerTail)
{
  MibrrGibbs mibrrGibbs;
  return mibrrGibbs.calcIncGamma(shape, cutVal, lowerTail);
}

// [[Rcpp::export]]
Eigen::VectorXd drawInvGauss(int n, double mu, double lambda)
{
  MibrrGibbs mibrrGibbs;
  VectorXd   out(n);
  for(int i = 0; i < n; i++) out[i] = mibrrGibbs.drawInvGauss(mu, lambda);
  return out;
}

// [[Rcpp::export]]
std::vector<int>
printObsIndices(Eigen::MatrixXd                 data,
		Eigen::VectorXd                 scales,
		std::vector< std::vector<int> > missIndices,
		Eigen::VectorXi                 respCounts,
		bool                            noMiss,
		int                             targetIndex)
{
  MibrrData mibrrData(data, scales, missIndices, respCounts, noMiss);
  return(mibrrData.getObsRows(targetIndex));
}

// [[Rcpp::export]]
std::vector<int>
printMissIndices(Eigen::MatrixXd                 data,
		 Eigen::VectorXd                 scales,
		 std::vector< std::vector<int> > missIndices,
		 Eigen::VectorXi                 respCounts,
		 bool                            noMiss,
		 int                             targetIndex)
{
  MibrrData mibrrData(data, scales, missIndices, respCounts, noMiss);
  return(mibrrData.getMissIndices(targetIndex));
}
