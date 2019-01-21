// Title:    Function definitions for the MibrrData Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2019-JAN-21
// Purpose:  This class contains the data-related functions used by the MIBRR
//           Gibbs sampler.

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

#include "MibrrData.h"

///////////////////////// CONSTRUCTORS / DESTRUCTOR ////////////////////////////

MibrrData::MibrrData(const MatrixXd        &data,
		     vector< vector<int> > missIndices,
		     const VectorXi        &respCounts,
		     const bool            noMiss) 
{
  _data        = data;
  _missIndices = missIndices;
  _respCounts  = respCounts;
  _noMiss      = noMiss;

  MatrixXd tmp = MatrixXd::Zero(respCounts.size(), nPreds());
  _means  = tmp;
  _scales = tmp;
}

MibrrData::~MibrrData() {}


///////////////////////////////// ACCESSORS ////////////////////////////////////


vector<int> MibrrData::getObsRows(int targetIndex) const
{
  int nObs = _data.rows();
    
  // Get an integer sequence of row indices:
  vector<int> allRows(nObs);
  std::iota(allRows.begin(), allRows.end(), 0);
    
  // Get the indices for rows with missing values on targetIndex:
  vector<int> missRows = _missIndices[targetIndex];
  
  // use std::set_difference to find the indices for rows with observed data on
  // targetIndex:
  vector<int>::iterator it;
  vector<int> obsRows(max(missRows.size(), allRows.size()));
  
  it = set_difference(allRows.begin(),
		      allRows.end(),
		      missRows.begin(),
		      missRows.end(),
		      obsRows.begin());
  
  obsRows.resize(it - obsRows.begin());
  
  return obsRows;
}



MatrixXd MibrrData::getIVs(int targetIndex, bool obsRows, bool standardize)
{
  int         nVars    = _data.cols();
  int         rowIndex = 0;
  int         nObs     = _data.rows();
  MatrixXd    outMat   = MatrixXd::Zero(nObs, nVars - 1);
  vector<int> useRows;
  
  if(_noMiss) {// Return full predictor matrix:
    outMat.leftCols(targetIndex) = _data.leftCols(targetIndex);
    
    outMat.rightCols(nVars - (targetIndex + 1)) =
      _data.rightCols(nVars - (targetIndex + 1));
  }
  else {
    if(obsRows) {// Return predictor rows corresponding to observed DVs:
      nObs    = _respCounts[targetIndex];
      useRows = getObsRows(targetIndex);
    }
    else {       // Return predictor rows corresponding to missing DVs:
      nObs    = _data.rows() - _respCounts[targetIndex];
      useRows = _missIndices[targetIndex];
    }
  
    outMat = MatrixXd::Zero(nObs, nVars - 1);
    
    for(int i : useRows) {
      outMat.block(rowIndex, 0, 1, targetIndex) =
	_data.block(i, 0, 1, targetIndex);
      
      outMat.block(rowIndex, targetIndex, 1, nVars - (targetIndex + 1)) = 
	_data.block(i, targetIndex + 1, 1, nVars - (targetIndex + 1));
      
      rowIndex++;
    }
  }
  
  if(standardize) {
    // Standardize the predictor matrix:
    outMat.rowwise() -= _means.row(targetIndex);
    outMat           *= _scales.row(targetIndex).asDiagonal().inverse();
  }
  return outMat;
}


VectorXd MibrrData::getDV(int targetIndex) const
{ 
  if(_noMiss) {
    return _data.col(targetIndex);
  }
  else {
    int         rowIndex = 0;
    VectorXd    outVec   = VectorXd::Zero(_respCounts[targetIndex]);
    vector<int> obsRows  = getObsRows(targetIndex);
    
    for(int i : obsRows) {
      outVec[rowIndex] = _data(i, targetIndex);
      rowIndex++;
    }
    return outVec;
  }
}


vector<int> MibrrData::getMissIndices(int targetIndex) const
{
  return _missIndices[targetIndex];
}


RowVectorXd MibrrData::getMeans(int targetIndex) const
{
  return _means.row(targetIndex);
}


RowVectorXd MibrrData::getScales(int targetIndex) const
{
  return _scales.row(targetIndex);
}


MatrixXd MibrrData::getData() const { return _data;                            }


///////////////////////////////// MUTATORS /////////////////////////////////////


void MibrrData::setData(const MatrixXd &newData) { _data = newData;            }


void MibrrData::setDV(const VectorXd &newDV, const int targetIndex)
{
  _data.col(targetIndex) = newDV;
}


void MibrrData::setElement(const double element, const int row, const int col)
{
  _data(row, col) = element;
}


void MibrrData::fillMissing(const MatrixXd &newTargets)
{
  int nTargets = newTargets.cols();
  
  for(int j = 0; j < nTargets; j++) {
    int         rowIndex = 0;
    vector<int> missRows = _missIndices[j];
    for(int i : missRows) {
      _data(i, j) = newTargets(rowIndex, j);
      rowIndex++;
    }
  }
}// END fillMissing()


void MibrrData::fillMissing(const VectorXd &newTarget, const int targetIndex)
{
  int         rowIndex = 0;
  vector<int> missRows = _missIndices[targetIndex];
  for(int i : missRows) {
    _data(i, targetIndex) = newTarget[rowIndex];
    rowIndex++;
  }
}// END fillMissing()


void MibrrData::updateMoments(const int targetIndex)
{
  _means.row(targetIndex)  = getIVs(targetIndex, true, false).colwise().mean();
  _scales.row(targetIndex) =
    (
     (
      getIVs(targetIndex, true, false).rowwise() - _means.row(targetIndex)
      ).colwise().squaredNorm() * double(1.0 / (nResp(targetIndex) - 1.0))
     ).cwiseSqrt();
}


//////////////////////////// DESCRIPTIVE FUNCTIONS /////////////////////////////


int MibrrData::nObs  ()                const { return _data.rows();            }
int MibrrData::nPreds()                const { return _data.cols() - 1;        }

int MibrrData::nResp (int targetIndex) const
{
  return _respCounts[targetIndex];
}

int MibrrData::nMiss (int targetIndex) const
{
  return _data.rows() - _respCounts[targetIndex];
}
