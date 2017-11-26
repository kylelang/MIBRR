// Title:    Test the MIBRR Sampling Routines
// Author:   Kyle M. Lang
// Created:  2017-NOV-23
// Modified: 2017-NOV-25
// Purpose:  Create a simple program to test the MIBRR samplers.

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

#include <string>
#include <iostream>
#include "MibrrSamplers.hpp"

using namespace std;

int main() {
  int    n = 10;
  double shape, cut, out;
  
  MibrrSamplers mibrrSamplers(235711);

  cout << "Inverse Gamma:" << endl;
  for(int i = 0; i < n; i++)
    cout << mibrrSamplers.drawInvGamma(1.0, 1.0) <<endl;

  cout << "\nIncomplete Gamma:" << endl;

  cout << "\nShape: ";
  cin >> shape;

  cout << "\nCut: ";
  cin >> cut;

  out = mibrrSamplers.calcIncGamma(shape, cut, false);
  cout << 100 * out << endl;

  VectorXd mu  = VectorXd::Zero(4);
  MatrixXd sig = MatrixXd::Identity(4, 4);

  cout << "\nMVN:" << endl;

  for(int i = 0; i < n; i++)
    cout << mibrrSamplers.drawMvn(mu, sig).transpose() <<endl;
  
  cout << "\nInverse Gaussian:" << endl;
  
  for(int i = 0; i < n; i++)
    cout << mibrrSamplers.drawInvGauss(1.0, 1.5) <<endl;

  cout << "\nGIG:" << endl;
  
  for(int i = 0; i < n; i++)
    cout << mibrrSamplers.drawGig(1.0, 1.0, 1.0) <<endl;
 
  return 0;
}
