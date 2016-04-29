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

//char MyErrors::getParent() const 
//{ 
//  return _parent; 
//}

bool MyErrors::getVerbosity() const 
{ 
  return _verboseErrors; 
}

bool MyErrors::getElasticNetFlag() const 
{ 
  return _useElasticNet; 
}


///// MUTATORS /////

//void MyErrors::setParent(char parent) 
//{ 
//  _parent = parent; 
//}

void MyErrors::setVerbosity(bool verboseErrors) 
{ 
  _verboseErrors = verboseErrors; 
}

void MyErrors::setElasticNetFlag(bool useElasticNet) 
{ 
  _useElasticNet = useElasticNet; 
}

//void MyErrors::doMiben()
//{
//  _useElasticNet = true;
//}

//void MyErrors::doMibl()
//{
//  _useElasticNet = false;
//}


///// ERROR HANDLING FUNCTIONS /////

void MyErrors::tauError(int errorCode) const
{
  if(errorCode == 1) {
    Rcpp::Rcerr << "\nSomething terrible has occured while updating Tau,\n";
    Rcpp::Rcerr << "and one of its mean values is non-positive.\n";
    Rcpp::Rcerr << "This program will now crash...have a nice day :)\n" << endl;
  }
  else if (errorCode == 2) {
    Rcpp::Rcerr << "\nSomething terrible has occured while updating Tau,\n";
    Rcpp::Rcerr << "and its scale value is non-positive.\n";
    Rcpp::Rcerr << "This program will now crash...have a nice day :)\n" << endl;
  }
  //if(_verboseErrors) {
  //  Rcpp::Rcerr << "P.S. I have written some useful information ";
  //  Rcpp::Rcerr << "in 'tauErrorLog.txt'." << endl;
  //}
  Rcpp::stop("Ouch! My tau is broken :(");
  
  // TODO: Put in code to write extra error information to text file.
}


void MyErrors::betaError(exception &e) const
{
  Rcpp::Rcerr << "\nSomething terrible has occured while updating Beta.\n";
  Rcpp::Rcerr << "While attempting to compute the moments of Beta's ";
  Rcpp::Rcerr << "conditional posterior distribtution,\nI seem to have broken ";
  Rcpp::Rcerr << "mathematics.\n\n Here is the exception I caught:\n";
  Rcpp::Rcerr << e.what() << "\n\n" << endl;
  Rcpp::Rcerr << "This program will now crash...have a nice day :)\n" << endl;

  //if(_verboseErrors) {
  //  Rcpp::Rcerr << "P.S. I have written some useful information ";
  //  Rcpp::Rcerr << "in 'betaErrorLog.txt'." << endl;
  //}
 Rcpp:stop("Oh no...beta luck next time ;)");
}


void MyErrors::lambdaError(int algCount,
			   string outPrefix,
			   string algName) const
{
  if(algCount > 3) {
    if(_verboseErrors) {
      Rcpp::Rcerr << "Lambda " << outPrefix << "optimization failed with "; 
      Rcpp::Rcerr << algName << ".\nNo further optimization algorithms ";
      Rcpp::Rcerr << "are available\nThis program will now crash" << endl;
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
      Rcpp::Rcerr << algName << ".\nNo further optimization algorithms ";
      Rcpp::Rcerr << "are available\nThis program will now crash" << endl;
    }
    Rcpp::stop("Lambda cannot be optimized.");
  }
  else {
    if(_verboseErrors) {
      cerr << e.what();
      Rcpp::Rcout << "Lambda " << outPrefix << "optimization failed with ";
      Rcpp::Rcout << algName << "\nRetrying with a different algorithm" << endl;
    }
  }
}
