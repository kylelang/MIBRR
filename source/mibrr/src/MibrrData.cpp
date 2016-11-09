// Title:    Function definitions for the MibrrData Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2016-NOV-08
// Purpose:  This class contains the data-related functions used by the MIBRR
//           Gibbs sampler.

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

#include "MibrrData.hpp"

///////////////////////// CONSTRUCTORS / DESTRUCTOR /////////////////////////////

MibrrData::MibrrData(const MatrixXd &data,
		     const VectorXd &dataScales,
		     vector<vector<int>>    missIndices,
		     const VectorXi &respCounts) 
{
  _data        = data;
  _dataScales  = dataScales;
  _missIndices = missIndices;
  _respCounts  = respCounts;
}

MibrrData::MibrrData() {}

MibrrData::~MibrrData() {}


///////////////////////////////// ACCESSORS /////////////////////////////////////


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



MatrixXd MibrrData::getObsIVs(int targetIndex) const
{ 
  int          nObs     = _data.rows();
  int          nVars    = _data.cols();
  int          rowIndex = 0;
  MatrixXd     tmpMat   = _data * _dataScales.asDiagonal().inverse(); 
  MatrixXd     outMat   = MatrixXd::Zero(_respCounts[targetIndex], nVars - 1);
  vector<int>  obsRows  = getObsRows(targetIndex);

  for(int i : obsRows) {
    outMat.block(rowIndex, 0, 1, targetIndex) =
      tmpMat.block(i, 0, 1, targetIndex);
    
    outMat.block(rowIndex, targetIndex, 1, nVars - (targetIndex + 1)) = 
      tmpMat.block(i, targetIndex + 1, 1, nVars - (targetIndex + 1));
    
    rowIndex++;
  }
  return outMat;
}



MatrixXd MibrrData::getMissIVs(int targetIndex) const
{ 
  int         nMiss    = _data.rows() - _respCounts[targetIndex];
  int         nVars    = _data.cols();
  int         rowIndex = 0;
  MatrixXd    tmpMat   = _data * _dataScales.asDiagonal().inverse(); 
  MatrixXd    outMat   = MatrixXd::Zero(nMiss, nVars - 1);
  vector<int> missRows = _missIndices[targetIndex];
  
  for(int i : missRows) {
    outMat.block(rowIndex, 0, 1, targetIndex) =
      tmpMat.block(i, 0, 1, targetIndex);
    
    outMat.block(rowIndex, targetIndex, 1, nVars - (targetIndex + 1)) = 
      tmpMat.block(i, targetIndex + 1, 1, nVars - (targetIndex + 1));
    
    rowIndex++;
  }
  return outMat;
}



MatrixXd MibrrData::getFullIVs(int targetIndex) const
{ 
  int      nObs   = _data.rows();
  int      nVars  = _data.cols();
  MatrixXd tmpMat = _data * _dataScales.asDiagonal().inverse(); 
  MatrixXd outMat = MatrixXd::Zero(nObs, nVars - 1);
  
  outMat.leftCols(targetIndex) = tmpMat.leftCols(targetIndex);
      
  outMat.rightCols(nVars - (targetIndex + 1)) =
    tmpMat.rightCols(nVars - (targetIndex + 1));
  
  return outMat;
}



VectorXd MibrrData::getObsDV(int targetIndex) const
{ 
  int          nObs     = _data.rows();
  int          rowIndex = 0;
  VectorXd     outVec   = VectorXd::Zero(_respCounts[targetIndex]);
  vector<int>  obsRows  = getObsRows(targetIndex);
  
  for(int i : obsRows) {
    outVec[rowIndex] = _data(i, targetIndex);
    rowIndex++;
  }
  return outVec; 
}



VectorXd MibrrData::getFullDV(int targetIndex) const
{
  return _data.col(targetIndex);
}



double MibrrData::getDataScales(int targetIndex) const
{
  return _dataScales[targetIndex];
}


MatrixXd MibrrData::getData()              const { return _data;                }
VectorXd MibrrData::getDataScales()        const { return _dataScales;          }


///////////////////////////////// MUTATORS //////////////////////////////////////


void MibrrData::setData(const MatrixXd &newData) { _data = newData;             }


void MibrrData::setDV(const VectorXd &newDV, const int targetIndex)
{
  _data.col(targetIndex) = newDV;
}


void MibrrData::setElement(const double element, const int row, const int col)
{
  _data(row, col) = element;
}


void MibrrData::computeDataScales()
{
  int nObs  = _data.rows();
  int nVars = _data.cols();
  
  for(int v = 0; v < nVars; v++) {
    double tmpMean = _data.col(v).mean();
    double tmpVar  =
      (_data.col(v).array() - tmpMean).square().sum() / double(nObs - 1);
    _dataScales[v] = sqrt(tmpVar);
  }
}


void MibrrData::fillMissing(const MatrixXd &newTargets)
{
  int nObs     = _data.rows();
  int nVars    = _data.cols();
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


//////////////////////////// DESCRIPTIVE FUNCTIONS //////////////////////////////


int MibrrData::nObs  ()                const { return _data.rows();             }
int MibrrData::nPreds()                const { return _data.cols() - 1;         }
int MibrrData::nResp (int targetIndex) const {return _respCounts[targetIndex];  }

int MibrrData::nMiss (int targetIndex) const
{
  return _data.rows() - _respCounts[targetIndex];
}


/////////////////////////// RANDOM VARIATE SAMPLERS /////////////////////////////


VectorXd MibrrData::drawMVN(VectorXd &meanVec, MatrixXd &covMat) const
{
  int      nVars = meanVec.size();
  MatrixXd covCholesky;
  VectorXd normDraws(nVars);
  
  covCholesky = covMat.llt().matrixL();
  for(int i = 0; i < nVars; i++) normDraws[i] = norm_rand();
  VectorXd testVec = covCholesky * normDraws;
  
  return meanVec + (covCholesky * normDraws);
}// END drawMVN()
