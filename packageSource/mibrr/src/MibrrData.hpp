// Title:    Header file for MibrrData Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2016-MAY-09
// Purpose:  This class contains data- and sampling-related functions
//           used by the MIBRR Gibbs sampler.

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

#ifndef MIBRRDATA_H
#define MIBRRDATA_H

#include "MibrrDefs.hpp"

using namespace std;
using namespace Eigen;

class MibrrData {
  
public:
  //////////////////////// CONSTRUCTORS / DESTRUCTOR ////////////////////////////
  
  MibrrData(MatrixXd&);
  
  MibrrData(MatrixXd&,
	    //VectorXd&,
	    VectorXd&,
	    double);
  
  ~MibrrData();
  
  //////////////////////////////// ACCESSORS ////////////////////////////////////
  
  MatrixXd getIVs(int);
  // @param:  the column-index of the current target variable
  // @return: the IVs of the imputation model with rows corresponding to
  //          missing DV observations deleted
  
  MatrixXd getFullIVs(int);
  // @param:  the column-index of the current target variable
  // @return: the full IVs matrix for the imputation model

  VectorXd getDV(int);
  // @param:  the column-index of the current target variable
  // @return: the (listwise deleted) DV of the imputation model
  
  VectorXd getFullDV(int);
  // @param:  the column-index of the current target variable
  // @return: the full DV vector for the imputation model

  MatrixXd getData();
  // @return: the data matrix

  ArrayXb getNonresponseVector(int);
  // @param:  the column-index of the current target variable  
  // @return: the appropriate nonresponse indicator vector

  ArrayXXb getNonresponseFilter();
  // @return: the nonresponse filter matrix

  //VectorXd getDataMeans();
  // @return: the column-wise means of the data
  
  //double getDataMeans(int);
  // @param:  a column index
  // @return: the column mean for the specified column index
  
  VectorXd getDataScales();
  // @return: the column-wise scales of the data
  
  double getDataScales(int);
  // @param:  a column index
  // @return: the  scale of the data in the specified column

  //////////////////////////////// MUTATORS /////////////////////////////////////

  void setData(MatrixXd&);
  // @param: a new data matrix
  
  void setDV(VectorXd&, int);
  //@param1: a new DV
  //@param2: index of the DV to replace

  void setElement(double newElement,
		  int rowIndex,
	          int colIndex);
  // @param1: a new data element
  // @param2: the row index of the element to replace
  // @param3: the column index of the element to replace

  void setMissingDataCode(double);
  // @param: a new missing data code

  void computeNonresponseFilter();
  // @effect: construct the nonresponse filter matrix

  void computeDataScales();
  // @effect: update the value of _dataScales based on imputed data
  
  void fillMissing(int);
  // @param:  the number of variables to be imputed
  // @effect: initially fill the missing values with (poor) imputaitons
  
  void fillMissing(MatrixXd&);
  // @param:  a matrix of new values to fill in missing target data
  // @effect: fill the missing values with the contents of the supplied matrix
  
  void fillMissing(VectorXd&, int);
  // @param1: a vector containing new imputations
  // @param2: the column index of the variable to fill
  // @effect: fill the missing values in the specified column with the values
  //          in the provided vector
  
  ////////////////////////// DESCRIPTIVE FUNCTIONS //////////////////////////////
  
  int nObs() const;
  // @return: the number of observations

  int nResponses(int);
  // @param:  the index for the target variable
  // @return: the number of responses (i.e., non-missing values)
  // for the target variable

  int nPreds() const;
  // @return: the number of independent variables in the imputation model

  ////////////////////////// RANDOM VARIATE SAMPLERS ////////////////////////////
  
  VectorXd drawMVN(VectorXd&,
		   MatrixXd&);
  // @param1: mean vector
  // @param2: covariance matrix
  // @return: random multivariate normal variates
  
private:
  MatrixXd _data;
  ArrayXXb _nonresponseFilter;
  VectorXd _responseCounts;
  //VectorXd _dataMeans;
  VectorXd _dataScales;
  double   _missingDataCode;
};

#endif
