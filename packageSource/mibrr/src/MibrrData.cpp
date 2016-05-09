// Title:    Function definitions for the MibrrData Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2016-MAY-09
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
  
MibrrData::MibrrData(MatrixXd &newData)
{
  _data     = newData;
  //_dataMeans  = RowVectorXd( _data.cols() );
  _dataScales = RowVectorXd( _data.cols() );
}

MibrrData::MibrrData(MatrixXd &newData,
		     //VectorXd &newDataMeans,
		     VectorXd &newDataScales,
		     double missingDataCode) 
{
  _missingDataCode   = missingDataCode;
  _data              = newData;
  //_dataMeans         = newDataMeans;
  _dataScales        = newDataScales;
  _nonresponseFilter = ArrayXXb::Constant(_data.rows(), _data.cols(), false);
  computeNonresponseFilter();
}

MibrrData::~MibrrData() 
{
}

///////////////////////////////// ACCESSORS /////////////////////////////////////

MatrixXd MibrrData::getIVs(int targetIndex)
{ 
  int nObs = _data.rows();
  int nVars = _data.cols();
  int newRowIndex = 0;
  MatrixXd outMat = MatrixXd::Zero(_responseCounts[targetIndex], nVars - 1);
  MatrixXd tmpMat = _data * _dataScales.asDiagonal().inverse(); 
  
  for(int i = 0; i < nObs; i++) {
    if(_nonresponseFilter(i, targetIndex) == false) {
      outMat.block(newRowIndex, 0, 1, targetIndex) =
	tmpMat.block(i, 0, 1, targetIndex);
      
      outMat.block(newRowIndex, targetIndex, 1, nVars - (targetIndex + 1)) = 
	tmpMat.block(i, targetIndex + 1, 1, nVars - (targetIndex + 1));
      
      newRowIndex++;
    }
  }
  return outMat;
}


MatrixXd MibrrData::getFullIVs(int targetIndex)
{ 
  int nObs = _data.rows();
  int nVars = _data.cols();
  MatrixXd outMat = MatrixXd::Zero(nObs, nVars - 1);
  MatrixXd tmpMat = _data * _dataScales.asDiagonal().inverse(); 
  
  outMat.leftCols(targetIndex) = tmpMat.leftCols(targetIndex);
      
  outMat.rightCols(nVars - (targetIndex + 1)) =
    tmpMat.rightCols(nVars - (targetIndex + 1));
  
  return outMat;
}


VectorXd MibrrData::getDV(int targetIndex)
{ 
  int nObs = _data.rows();
  int newRowIndex = 0;
  VectorXd outVec = VectorXd::Zero(_responseCounts[targetIndex]);
  
  for(int i = 0; i < nObs; i++) {
    if(_nonresponseFilter(i, targetIndex) == false) {
      outVec[newRowIndex] = _data(i, targetIndex);
      newRowIndex++;
    }
  }
  return outVec; 
}


VectorXd MibrrData::getFullDV(int targetIndex)
{
  return _data.col(targetIndex);
}


ArrayXb MibrrData::getNonresponseVector(int targetIndex) 
{ 
  return _nonresponseFilter.col(targetIndex); 
}


ArrayXXb MibrrData::getNonresponseFilter()
{
  return _nonresponseFilter;
}


MatrixXd MibrrData::getData()
{
  return _data;
}


//VectorXd MibrrData::getDataMeans()
//{
//  return _dataMeans;
//}


//double MibrrData::getDataMeans(int targetIndex)
//{
//  return _dataMeans[targetIndex];
//}


VectorXd MibrrData::getDataScales()
{
  return _dataScales;
}


double MibrrData::getDataScales(int targetIndex)
{
  return _dataScales[targetIndex];
}

///////////////////////////////// MUTATORS //////////////////////////////////////

void MibrrData::setData(MatrixXd &newData) 
{
  _data = newData;
}


void MibrrData::setDV(VectorXd &newDV,
		      int targetIndex)
{
  _data.col(targetIndex) = newDV;
}


void MibrrData::setElement(double newElement,
			   int rowIndex,
			   int colIndex)
{
  _data(rowIndex, colIndex) = newElement;
}


void MibrrData::setMissingDataCode(double missingDataCode)
{
  _missingDataCode = missingDataCode;
}


void MibrrData::computeNonresponseFilter() 
{
  for(int i = 0; i < _data.rows(); i++) {
    for(int j = 0; j < _data.cols(); j++) {
      bool isMissing = _data(i, j) == _missingDataCode;
      if(isMissing) _nonresponseFilter(i, j) = true;
    }
  }
  
  VectorXd nObsVec = VectorXd::Constant(_data.cols(), _data.rows()).array();
  VectorXd nMissingVec = _nonresponseFilter.colwise().count().cast<double>();
  
  _responseCounts = nObsVec - nMissingVec;
}


void MibrrData::computeDataScales()
{
  for(int i = 0; i < _data.cols(); i++) {
    double tmpMean = _data.col(i).mean();
    double tmpVar = (_data.col(i).array() - tmpMean).square().mean();
    _dataScales[i] = sqrt(tmpVar);
  }
}


void MibrrData::fillMissing(int nTargets)
{
  int nObs = _data.rows();
  int nVars = _data.cols();
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
      dvMeans[v] = dvMeans[v] / _responseCounts[v];

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
    double dvMean = _data.col(0).sum() / _responseCounts[0];
    for(int i = 0; i < nObs; i++) {
      if(_nonresponseFilter(i, 0))
	_data(i, 0) = R::rnorm(dvMean, _dataScales[0]);
    }
  }
}// END fillMissing()


void MibrrData::fillMissing(MatrixXd &newTargets)
{
  int nObs = _data.rows();
  int nVars = _data.cols();
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


void MibrrData::fillMissing(VectorXd &newTarget,
			    int targetIndex)
{
  int nObs = _data.rows();
  for(int i = 0; i < nObs; i++) {
    if(_nonresponseFilter(i, targetIndex))
      _data(i, targetIndex) = newTarget[i];
  }
}// END fillMissing()

//////////////////////////// DESCRIPTIVE FUNCTIONS //////////////////////////////

int MibrrData::nObs() const 
{ 
  return _data.rows(); 
}


int MibrrData::nResponses(int targetIndex)
{
  return _responseCounts[targetIndex];
}


int MibrrData::nPreds() const 
{ 
  return _data.cols() - 1; 
}

/////////////////////////// RANDOM VARIATE SAMPLERS /////////////////////////////
  
VectorXd MibrrData::drawMVN(VectorXd &meanVec, 
			    MatrixXd &covMat)
{
  int nVars = meanVec.size();
  MatrixXd covCholesky;
  VectorXd normDraws(nVars);
  
  covCholesky = covMat.llt().matrixL();
  for(int i = 0; i < nVars; i++) normDraws[i] = norm_rand();
  VectorXd testVec = covCholesky * normDraws;
  
  return meanVec + (covCholesky * normDraws);
}// END drawMVN()
