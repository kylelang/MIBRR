// Title:    Header file for MyData Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2015-FEB-24
// Purpose:  This class contains data- and sampling-related functions
//           used by the MIBRR Gibbs sampler.

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

#ifndef MYDATA_H
#define MYDATA_H


#include "MyParams.hpp"


using namespace std;
using namespace Eigen;


class MyData {

 public:

  ///// CONSTRUCTORS / DESTRUCTOR /////
  
  MyData(MatrixXd&);

  MyData(MatrixXd&,
	 VectorXd&,
	 VectorXd&,
	 double);
  
  ~MyData();


  ///// ACCESSORS /////
  
  MatrixXd getIVs(int);
  // @param the column-index of the current target variable
  // @return the (listwise deleted) IVs of the imputation model

  MatrixXd getFullIVs(int);
  // @param the column-index of the current target variable
  // @return the full IVs matrix for the imputation model

  VectorXd getDV(int);
  // @param the column-index of the current target variable
  // @return the (listwise deleted) DV of the imputation model

  VectorXd getFullDV(int);
  // @param the column-index of the current target variable
  // @return the full DV vector for the imputation model

  MatrixXd getMyData();
  // @return the data matrix

  ArrayXb getNonresponseVector(int);
  // @param the column-index of the current target variable  
  // @return the appropriate nonresponse indicator vector

  ArrayXXb getNonresponseFilter();
  // @return the nonresponse filter matrix

  VectorXd getDataMeans();

  double getDataMeans(int);

  VectorXd getDataScales();

  double getDataScales(int);

  ///// MUTATORS /////

  void setMyData(MatrixXd&);
  // @param a new data matrix
  
  void setDV(VectorXd&, int);
  //@param1 a new DV
  //@param2 index of the DV to replace

  void setElement(double newElement,
		  int rowIndex,
	          int colIndex);
  // @param1 a new data element
  // @param2 the row index of the element to replace
  // @param3 the column index of the element to replace

  void setMissingDataCode(double);
  // @param a new missing data code

  void computeNonresponseFilter();

  void fillMissing(int);
  // @param the number of variables to be imputed
  
  void fillMissing(MatrixXd&);
  // @param a matrix of new values to fill in missing target data

  void fillMissing(VectorXd&, int);
  // @param1 a vector containing new imputations
  // @param2 the column index of the variable to fill


  ///// DESCRIPTIVE FUNCTIONS /////

  int nObs() const;
  // @return the number of observations

  int nResponses(int);
  // @param the index for the target variable
  // @return the number of responses (i.e., non-missing values)
  // for the target variable

  int nPreds() const;
  // @return the number of independent variables 
  // in the imputation model


 ///// RANDOM VARIATE SAMPLERS /////

  double drawInvGamma(double,
		      double);
  // Draw a random Inver Gamma variate
  // @param1 shape parameter
  // @param2 scale parameter

  VectorXd drawMVN(VectorXd&,
		   MatrixXd&);
  // Draw a vector of random multivariate normal variates
  // @param1 mean vector
  // @param2 covariance matrix
  
  double calcIncGamma(double,
		      double,
		      bool);
  // Calculate the incomplete gamma function (i.e., the un-normalized 
  // area under the upper or lower tail of the Gamma CDF).
  // @param1 the shape parameter of the underlying gamma distribution
  // @param2 the threshold value cutting off the upper or lower tail
  // (i.e., the underlying variate whose probability or [1 - probability] 
  // is being returned)
  // @param3 a logic switch: true gives lower incomplete gamma, false
  // gives upper incomplete gamma
  
  double drawInvGauss(double,
		      double);
  // Draw a random variate from the inverse Gaussian distribution
  // @param1 the mean parameter (mu)
  // @param2 the shape parameter (lambda)


 private:
  
  MatrixXd _myData;
  ArrayXXb _nonresponseFilter;
  VectorXd _responseCounts;
  VectorXd _dataMeans;
  VectorXd _dataScales;
  double _missingDataCode;
};


#endif
