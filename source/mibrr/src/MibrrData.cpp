// Title:    Function definitions for the MibrrData Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2016-NOV-05
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

MibrrData::MibrrData(const MatrixXd &newData)
{
  _data       = newData;
  _dataScales = RowVectorXd( _data.cols() );
}


MibrrData::MibrrData(const MatrixXd &newData,
		     const VectorXd &newDataScales,
		     const double   missingDataCode) 
{
  _missingDataCode   = missingDataCode;
  _data              = newData;
  _dataScales        = newDataScales;
  _nonresponseFilter = ArrayXXb::Constant(_data.rows(), _data.cols(), false);
  computeNonresponseFilter();
}

MibrrData::MibrrData() {}

MibrrData::~MibrrData() {}


///////////////////////////////// ACCESSORS /////////////////////////////////////


MatrixXd MibrrData::getIVs(int targetIndex) const
{ 
  int      nObs     = _data.rows();
  int      nVars    = _data.cols();
  int      rowIndex = 0;
  MatrixXd tmpMat   = _data * _dataScales.asDiagonal().inverse(); 
  MatrixXd outMat   = MatrixXd::Zero(_respCounts[targetIndex], nVars - 1);
  
  for(int i = 0; i < nObs; i++) {
    if(_nonresponseFilter(i, targetIndex) == false) {
      outMat.block(rowIndex, 0, 1, targetIndex) =
	tmpMat.block(i, 0, 1, targetIndex);
      
      outMat.block(rowIndex, targetIndex, 1, nVars - (targetIndex + 1)) = 
	tmpMat.block(i, targetIndex + 1, 1, nVars - (targetIndex + 1));
      
      rowIndex++;
    }
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


VectorXd MibrrData::getDV(int targetIndex) const
{ 
  int      nObs     = _data.rows();
  int      rowIndex = 0;
  VectorXd outVec   = VectorXd::Zero(_respCounts[targetIndex]);
  
  for(int i = 0; i < nObs; i++) {
    if(_nonresponseFilter(i, targetIndex) == false) {
      outVec[rowIndex] = _data(i, targetIndex);
      rowIndex++;
    }
  }
  return outVec; 
}


VectorXd MibrrData::getFullDV(int targetIndex) const
{
  return _data.col(targetIndex);
}


ArrayXb MibrrData::getNonresponseVector(int targetIndex) const 
{ 
  return _nonresponseFilter.col(targetIndex); 
}


double MibrrData::getDataScales(int targetIndex) const
{
  return _dataScales[targetIndex];
}


ArrayXXb MibrrData::getNonresponseFilter() const { return _nonresponseFilter;   }
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


void MibrrData::setMissingDataCode(const double code)
{
  _missingDataCode = code;
}


void MibrrData::computeNonresponseFilter() 
{
  for(int i = 0; i < _data.rows(); i++) {
    for(int j = 0; j < _data.cols(); j++) {
      bool isMissing = _data(i, j) == _missingDataCode;
      if(isMissing) _nonresponseFilter(i, j) = true;
    }
  }
  
  VectorXd nObsVec  = VectorXd::Constant(_data.cols(), _data.rows()).array();
  VectorXd nMissVec = _nonresponseFilter.colwise().count().cast<double>();
  
  _respCounts = nObsVec - nMissVec;
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


void MibrrData::fillMissing(const int nTargets)
{
  int      nObs  = _data.rows();
  int      nVars = _data.cols();
  MatrixXd mvnDraws(nObs, nTargets);

  // Targets are also zero imputed, initially, to compute dvSigma
  for(int i = 0; i < nObs; i++) {
    for(int j = 0; j < nVars; j++) {
      if(_nonresponseFilter(i, j)) _data(i, j) = 0.0;
    }
  }
  
  if(nTargets > 1) {
    // Compute the mean of non-missing target values:
    VectorXd dvMeans = _data.leftCols(nTargets).colwise().sum();
    for(int v = 0; v < nTargets; v++)
      dvMeans[v] = dvMeans[v] / _respCounts[v];

    // Compute the target covariances with pairwise deletion
    MatrixXd dvSigma = _data.leftCols(nTargets).transpose() *
      _data.leftCols(nTargets) / double(nObs - 1);
    
    for(int i = 0; i < nObs; i++) {
      mvnDraws.row(i) = drawMVN(dvMeans, dvSigma);
    }
    
    for(int i = 0; i < nObs; i++) {
      for(int j = 0; j < nTargets; j++) {
	if( _nonresponseFilter(i, j) ) _data(i, j) = mvnDraws(i, j);
      }
    }  
  }
  else {
    double dvMean = _data.col(0).sum() / _respCounts[0];
    for(int i = 0; i < nObs; i++) {
      if(_nonresponseFilter(i, 0))
	_data(i, 0) = R::rnorm(dvMean, _dataScales[0]);
    }
  }
}// END fillMissing()


void MibrrData::fillMissing(const MatrixXd &newTargets)
{
  int nObs     = _data.rows();
  int nVars    = _data.cols();
  int nTargets = newTargets.cols();

  for(int i = 0; i < nObs; i++) {
    for(int j = 0; j < nTargets; j++) {
      if(_nonresponseFilter(i, j)) _data(i, j) = newTargets(i, j);
    }
  }
  
  for(int i = 0; i < nObs; i++) {
    for(int j = nTargets; j < nVars; j++) {
      if(_nonresponseFilter(i, j)) _data(i, j) = 0.0;
    }
  }
}// END fillMissing()


void MibrrData::fillMissing(const VectorXd &newTarget, const int targetIndex)
{
  int nObs = _data.rows();
  for(int i = 0; i < nObs; i++) {
    if(_nonresponseFilter(i, targetIndex)) _data(i, targetIndex) = newTarget[i];
  }
}// END fillMissing()


//////////////////////////// DESCRIPTIVE FUNCTIONS //////////////////////////////


int MibrrData::nObs()   const { return _data.rows();                            }
int MibrrData::nPreds() const { return _data.cols() - 1;                        }

int MibrrData::nResponses(int targetIndex) const
{
  return _respCounts[targetIndex];
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
