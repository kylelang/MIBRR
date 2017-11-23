// Title:    Test the GIG Sampling Routines
// Author:   Kyle M. Lang
// Created:  2017-NOV-23
// Modified: 2017-NOV-23
// Purpose:  Create a simple program to test the GIG sampler.

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
#include "gigSampler.hpp"

using namespace std;

int main() {
  int    n;
  double lambda, chi, psi;
  std::vector<double> output(n);
  
  cout << "How many variates would you like? ";
  cin >> n;
  
  cout << "What should 'lambda' be? ";
  cin >> lambda;
  
  cout << "What should 'chi' be? ";
  cin >> chi;
  
  cout << "What should 'psi' be? ";
  cin >> psi;
  
  output = drawGig(n, lambda, chi, psi);
  
  for(int i = 0; i < n; i++) cout << output[i] << endl;
  
  return 0;
}
