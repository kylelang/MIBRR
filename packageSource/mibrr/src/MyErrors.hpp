// Title:    Header file for the MyErrors Class
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

#ifndef MYERRORS_H
#define MYERRORS_H

#include "MyParams.hpp"
//#include "MyData.hpp"
//#include <nlopt.hpp>

using namespace std;

class MyErrors {

 public:

  ///// CONSTRUCTORS / DESTRUCTOR /////

  MyErrors();
  
  ~MyErrors();


  ///// ACCESSORS /////

  bool getVerbosity() const;
  // @returnL the current verbosity level for errors

  bool getElasticNetFlag() const;
  // @return: the current value of the _useElasticNet switch

  
  ///// MUTATORS /////

  void setVerbosity(bool);
  // @param: a new verbosity level flag

  void setElasticNetFlag(bool);
  // @param: a new value of the _useElasticNet flag

  
  ///// EXCEPTION HANDLING FUNCTIONS /////
  
  void tauError(int) const;
  // @param:  the orignial error code
  // @effect: send the appropriate error message to stderr
  
  void betaError(exception&) const;
  // @param:  the orignial exception
  // @effect: send the appropriate error message to stderr
  
  void lambdaError(int, string, string) const;
  // @param1: the ordinality of the current optimization method
  // @param2: the pre-optimize/optimize prefix
  // @param3: the name of the current optimization method
  // @effect: send the appropriate error message to stderr

  void lambdaError(exception&, int, string, string) const;
  // @param1: the exception thrown by nlopt
  // @param2: the ordinality of the current optimization method
  // @param3: the pre-optimize/optimize prefix
  // @param4: the name of the current optimization method
  // @effect: send the appropriate error message to stderr
  
 protected:
  bool _verboseErrors;
  bool _useElasticNet;
};

#endif
