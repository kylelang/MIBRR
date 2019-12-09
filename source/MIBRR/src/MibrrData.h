// Title:    Header file for MibrrData Class
// Author:   Kyle M. Lang
// Created:  2014-08-24
// Modified: 2019-01-21
// Purpose:  This class contains data- and sampling-related functions used by
//           the MIBRR Gibbs sampler.

//--------------------- COPYRIGHT & LICENSING INFORMATION --------------------//
//  Copyright (C) 2019 Kyle M. Lang <k.m.lang@uvt.nl>                         //  
//                                                                            //
//  This file is part of MIBRR.                                               //
//                                                                            //
//  This program is free software: you can redistribute it and/or modify it   //
//  under the terms of the GNU General Public License as published by the     //
//  Free Software Foundation, either version 3 of the License, or (at you     //
//  option) any later version.                                                //
//                                                                            //
//  This program is distributed in the hope that it will be useful, but       //
//  WITHOUT ANY WARRANTY; without even the implied warranty of                //
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General  //
//  Public License for more details.                                          //
//                                                                            //
//  You should have received a copy of the GNU General Public License along   //
//  with this program. If not, see <http://www.gnu.org/licenses/>.            //
//----------------------------------------------------------------------------//

#ifndef MIBRRDATA_H
#define MIBRRDATA_H

#include "MibrrDefs.h"
#include <algorithm>

class MibrrData {
  
public:
  //////////////////////// CONSTRUCTORS / DESTRUCTOR ///////////////////////////
    
  MibrrData(const MatrixXd&,
	    vector< vector<int> >,
	    const VectorXi&,
	    const bool);
  // @param1: data matrix
  // @param2: list of indices for missing rows
  // @param3: vector or response counts
  // @param4: logical flag denoting completely observed data
    
  ~MibrrData();

  
  //////////////////////////////// ACCESSORS ///////////////////////////////////

  
  vector<int> getObsRows(int) const;
  // @param:  the column-index of the current target variable
  // @return: the row indices for the observed rows of the target variable
  
  MatrixXd getIVs(int, bool, bool);
  // @param1: the column-index of the current target variable
  // @param2: return the rows corresponding to observed entries in Y or return
  //          the rows corresponding to the missing entries in Y?
  // @param3: standardize the predictor matrix before returning or not?
  // @return: the IVs of the imputation model with rows corresponding to
  //          missing DV observations deleted

  VectorXd getDV(int) const;
  // @param:  the column-index of the current target variable
  // @return: the (listwise deleted) DV of the imputation model
  
  vector<int> getMissIndices(int targetIndex) const;
  // @param:  a column index
  // @return: the row indices of missing values in the specified column
  
  RowVectorXd getMeans(int) const;
  // @param:  column index of the target variable
  // @return: column-wise means of the training data for the target variable
  
  RowVectorXd getScales(int) const;
  // @param:  column index of the target variable
  // @return: column-wise scales of the training data for the target variable

  MatrixXd getData() const;
  // @return: the data matrix

  
  //////////////////////////////// MUTATORS ////////////////////////////////////

  
  void setData(const MatrixXd&);
  // @param: a new data matrix
  
  void setDV(const VectorXd&, const int);
  //@param1: a new DV
  //@param2: index of the DV to replace

  void setElement(const double, const int, const int);
  // @param1: a new data element
  // @param2: the row index of the element to replace
  // @param3: the column index of the element to replace
  
  void fillMissing(const MatrixXd&);
  // @param:  a matrix of new values to fill in missing target data
  // @effect: fill the missing values with the contents of the supplied matrix
  
  void fillMissing(const VectorXd&, const int);
  // @param1: a vector containing new imputations
  // @param2: the column index of the variable to fill
  // @effect: fill the missing values in the specified column with the values
  //          in the provided vector

  void updateMoments(const int);
  // @param:  column index of the target variable
  // @effect: updated the centers and scales for the target variables training
  //          set predictors

  
  ////////////////////////// DESCRIPTIVE FUNCTIONS /////////////////////////////

  
  int nObs() const;
  // @return: the number of observations

  int nPreds() const;
  // @return: the number of independent variables in the imputation model

  int nResp(int) const;
  // @param:  the index for the target variable
  // @return: the number of responses (i.e., non-missing values) for the target
  //          variable

  int nMiss(int) const;
  // @param:  the index for the target variable
  // @return: the number of missing values for the target variable
  
private:
  bool                  _noMiss;
  MatrixXd              _data;
  MatrixXd              _means;
  MatrixXd              _scales;
  VectorXi              _respCounts;
  vector< vector<int> > _missIndices;
};

#endif
