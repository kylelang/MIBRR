// Title:    Function definitions for the MyErrors Class
// Author:   Kyle M. Lang
// Created:  2016-APR-29
// Modified: 2016-APR-29
// Purpose:  This class contains the C++ exception handling functions for the
//           mibrr package.

//--------------------- COPYRIGHT & LICENSING INFORMATION ---------------------//
//  Copyright (C) 2016 Kyle M. Lang <kyle.lang@ttu.edu>                        //  
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

#include "MyErrors.hpp"

MyErrors::MyErrors() 
{
  _verboseErrors = true;
  _useElasticNet = true;
}

MyErrors::~MyErrors() 
{
}


///// ACCESSORS /////

bool MyErrors::getVerbosity() const 
{ 
  return _verboseErrors; 
}

bool MyErrors::getElasticNetFlag() const 
{ 
  return _useElasticNet; 
}


///// MUTATORS /////

void MyErrors::setVerbosity(bool verboseErrors) 
{ 
  _verboseErrors = verboseErrors; 
}

void MyErrors::setElasticNetFlag(bool useElasticNet) 
{ 
  _useElasticNet = useElasticNet; 
}


///// ERROR HANDLING FUNCTIONS /////

void MyErrors::tauError(int errorCode) const
{
  if(errorCode == 1) {
    Rcpp::Rcout << "\n";
    Rcpp::stop("Ouch! My tau is broken :(\nSomething terrible has occured \
while updating Tau,\nand one of its mean values is non-positive.\n");
  }
  else if (errorCode == 2) {
    Rcpp::Rcout << "\n";
    Rcpp::stop("Ouch! My tau is broken :(\nSomething terrible has occured \
while updating Tau,\nand one of its mean values is non-positive.\n");
  }
}


void MyErrors::betaError(exception &e) const
{
  Rcpp::Rcerr << e.what() << endl;
  Rcpp::stop("Something terrible has occured while updating Beta.\nAbove this \
message, I've printed the that exception I caught.\nBeta luck next time ;)");
}


void MyErrors::lambdaError(int algCount,
			   string outPrefix,
			   string algName) const
{
  if(algCount > 3) {
    if(_verboseErrors) {
      Rcpp::Rcerr << "Lambda " << outPrefix << "optimization failed with "; 
      Rcpp::Rcerr << algName << ".\nNo further optimization algorithms ";
      Rcpp::Rcerr << "are available." << endl;
    }
    Rcpp::stop("Lambda cannot be optimized.");
  }
  else {
    if(_verboseErrors) {
      Rcpp::Rcout << "Lambda " << outPrefix << "optimization failed with ";
      Rcpp::Rcout << algName << "\nRetrying with a different algorithm" << endl;
    }
  }
}


void MyErrors::lambdaError(exception &e,
			   int algCount,
			   string outPrefix,
			   string algName) const
{
  if(algCount > 3) {
    if(_verboseErrors) {
      Rcpp::Rcerr << e.what();
      Rcpp::Rcerr << "Lambda " << outPrefix << "optimization failed with "; 
      Rcpp::Rcerr << algName << ", and returned the preceding exception.\nNo";
      Rcpp::Rcerr << "no further optimization algorithms are available." << endl;
    }
    Rcpp::stop("Lambda cannot be optimized.");
  }
  else {
    if(_verboseErrors) {
      cerr << e.what();
      Rcpp::Rcout << "Lambda " << outPrefix << "optimization failed with ";
      Rcpp::Rcout << algName << ", and returned the preceding exception.\n";
      Rcpp::Rcout << "Retrying with a different algorithm" << endl;
    }
  }
}
