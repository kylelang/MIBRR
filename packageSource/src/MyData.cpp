// Title:    Function definitions for the MyData Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2015-FEB-23
// Purpose:  This class contains the data- and sampling-related functions
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

#include "MyData.hpp"
#include "MyParams.hpp"

MyData::MyData(MatrixXd &newData)
{
  _myData = newData;
  _dataMeans = RowVectorXd( _myData.cols() );
  _dataScales = RowVectorXd( _myData.cols() );
}

MyData::MyData(MatrixXd &newData,
	       VectorXd &newDataMeans,
	       VectorXd &newDataScales,
	       double missingDataCode) 
{
  _missingDataCode = missingDataCode;
  _myData = newData;
  _dataMeans = newDataMeans;
  _dataScales = newDataScales;
  _nonresponseFilter = 
    ArrayXXb::Constant(_myData.rows(),
		       _myData.cols(),
		       false);
  computeNonresponseFilter();
}

MyData::~MyData() 
{
}

MatrixXd MyData::getIVs(int targetIndex)
{ 
  int nObs = _myData.rows();
  int nVars = _myData.cols();
  int newRowIndex = 0;

  MatrixXd outMat = 
    MatrixXd::Zero( _responseCounts[targetIndex], 
		    nVars - 1 );
  
  MatrixXd tmpMat = _myData * 
    _dataScales.asDiagonal().inverse(); 
  
  for(int i = 0; i < nObs; i++) {
    if(_nonresponseFilter(i, targetIndex) == false) {
      outMat.block(newRowIndex, 0, 1, targetIndex) = 
	tmpMat.block(i, 0, 1, targetIndex);
      
      outMat.block(newRowIndex, targetIndex,
		   1, nVars - (targetIndex + 1) ) = 
	tmpMat.block(i, targetIndex + 1,
		     1, nVars - (targetIndex + 1) );
      newRowIndex++;
    }
  }
  return outMat;
}

MatrixXd MyData::getFullIVs(int targetIndex)
{ 
  int nObs = _myData.rows();
  int nVars = _myData.cols();

  MatrixXd outMat = 
    MatrixXd::Zero(nObs, nVars - 1);
  
  MatrixXd tmpMat = _myData * 
    _dataScales.asDiagonal().inverse(); 
  
  outMat.leftCols(targetIndex) = 
    tmpMat.leftCols(targetIndex);
      
  outMat.rightCols( nVars - (targetIndex + 1) ) = 
    tmpMat.rightCols( nVars - (targetIndex + 1) );
  
  return outMat;
}

VectorXd MyData::getDV(int targetIndex)
{ 
  int nObs = _myData.rows();
  int newRowIndex = 0;

  VectorXd outVec = 
    VectorXd::Zero(_responseCounts[targetIndex]);
  
  for(int i = 0; i < nObs; i++) {
    if(_nonresponseFilter(i, targetIndex) == false) {
      outVec[newRowIndex] = _myData(i, targetIndex);
      newRowIndex++;
    }
  }
  return outVec; 
}

VectorXd MyData::getFullDV(int targetIndex)
{ 
  return _myData.col(targetIndex); 
}

ArrayXb MyData::getNonresponseVector(int targetIndex) 
{ 
  return _nonresponseFilter.col(targetIndex); 
}

ArrayXXb MyData::getNonresponseFilter()
{
  return _nonresponseFilter;
}

MatrixXd MyData::getMyData()
{
  return _myData;
}

VectorXd MyData::getDataMeans()
{
  return _dataMeans;
}

double MyData::getDataMeans(int targetIndex)
{
  return _dataMeans[targetIndex];
}

VectorXd MyData::getDataScales()
{
  return _dataScales;
}

double MyData::getDataScales(int targetIndex)
{
  return _dataScales[targetIndex];
}

void MyData::setMyData(MatrixXd &newData) 
{
  _myData = newData;
}

void MyData::setDV(VectorXd &newDV,
		   int targetIndex)
{
  _myData.col(targetIndex) = newDV;
}

void MyData::setElement(double newElement,
			int rowIndex,
			int colIndex)
{
  _myData(rowIndex, colIndex) = newElement;
}

void MyData::setMissingDataCode(double missingDataCode)
{
  _missingDataCode = missingDataCode;
}

void MyData::computeNonresponseFilter() 
{

  for(int i = 0; i < _myData.rows(); i++) {
    for(int j = 0; j < _myData.cols(); j++) {
      bool isMissing = _myData(i, j) == _missingDataCode;
      if(isMissing) 
	_nonresponseFilter(i, j) = true;
    }
  }

  VectorXd nObsVec = 
    VectorXd::Constant( _myData.cols(), _myData.rows() ).array();

  VectorXd nMissingVec = 
    _nonresponseFilter.colwise().count().cast<double>();
  
  _responseCounts = nObsVec - nMissingVec;

}

void MyData::fillMissing(int nTargets)
{
  int nObs = _myData.rows();
  int nVars = _myData.cols();
  MatrixXd mvnDraws(nObs, nTargets);

  for(int i = 0; i < nObs; i++) {
    for(int j = 0; j < nVars; j++) {
      if( _nonresponseFilter(i, j) )
	_myData(i, j) = 0.0;
    }
  }

  // Targets are also zero imputed, initially, so that 
  // dvSigma is computed appropriately
  if(nTargets > 1) {
    VectorXd dvMeans = _dataMeans.head(nTargets);
    MatrixXd dvSigma = _myData.leftCols(nTargets).transpose() * 
      _myData.leftCols(nTargets) / double(nObs - 1);
    
    for(int i = 0; i < nObs; i++) {
      mvnDraws.row(i) = drawMVN(dvMeans, dvSigma);
    }
    
    for(int i = 0; i < nObs; i++) {
      for(int j = 0; j < nTargets; j++) {
	if( _nonresponseFilter(i, j) )
	  _myData(i, j) = mvnDraws(i, j);
      }
    }  
  }
  else {
    for(int i = 0; i < nObs; i++) {
      if( _nonresponseFilter(i, 0) )
	_myData(i, 0) = R::rnorm(_dataMeans[0], _dataScales[0]);
    }
  }

}// END fillMissing()

void MyData::fillMissing(MatrixXd &newTargets)
{
  int nObs = _myData.rows();
  int nVars = _myData.cols();
  int nTargets = newTargets.cols();

  for(int i = 0; i < nObs; i++) {
    for(int j = 0; j < nTargets; j++) {
      if( _nonresponseFilter(i, j) )
	_myData(i, j) = newTargets(i, j);
    }
  }
  
  for(int i = 0; i < nObs; i++) {
    for(int j = nTargets; j < nVars; j++) {
      if( _nonresponseFilter(i, j) )
	_myData(i, j) = 0.0;
    }
  }
}// END fillMissing()

void MyData::fillMissing(VectorXd &newTarget,
			 int targetIndex)
{
  int nObs = _myData.rows();
  for(int i = 0; i < nObs; i++) {
    if( _nonresponseFilter(i, targetIndex) )
      _myData(i, targetIndex) = newTarget[i];
  }
}// END fillMissing()

int MyData::nObs() const 
{ 
  return _myData.rows(); 
}

int MyData::nResponses(int targetIndex)
{
  return _responseCounts[targetIndex];
}

int MyData::nPreds() const 
{ 
  return _myData.cols() - 1; 
}


///// SAMPLING FUNCTIONS /////

double MyData::drawInvGamma(double shape, 
			    double scale)
{
  return 1.0 / R::rgamma(shape, 1.0 / scale);
}//END drawInvGamma()

VectorXd MyData::drawMVN(VectorXd &meanVec, 
			 MatrixXd &covMat)
{
  int nVars = meanVec.size();
  MatrixXd covCholesky;
  VectorXd normDraws(nVars);
  
  covCholesky = covMat.llt().matrixL();
  
  for(int i = 0; i < nVars; i++) 
    normDraws[i] = norm_rand();
  
  VectorXd testVec = covCholesky * normDraws;
  
  return meanVec + (covCholesky * normDraws);
}// END drawMVN()


double MyData::calcIncGamma(double shape, 
			    double cutVal,
			    bool lowerTail)
{
  double scale = 1.0;
  int lower = (int)lowerTail;
  int logTran = 0; // Don't want log transform
  
  return R::pgamma(cutVal, shape, scale, lower, logTran) * 
    tgamma(shape);
}// END calcIncGamma()


// The following function was adapted from source code originally
// implemented by Robert E. Wheeler (2001-MAR) in the R package SuppDists:

// Note from the original author:
//// Random inverse Gaussian values
//// Follows Mitchael,J.R., Schucany, W.R. and Haas, R.W. (1976). Generating
//// random roots from variates using transformations with multiple roots.
//// American Statistician. 30-2. 88-91.

double MyData::drawInvGauss(double mu, 
			    double lambda)
{ 
  double b = 0.5 * mu / lambda;
  double a = mu * b;
  double c = 4.0 * mu * lambda;
  double d = pow(mu, 2);
  double outVal = 0.0;

  while(outVal <= 0.0) {
    double tmpDraw = norm_rand();
    double v = pow(tmpDraw, 2); // Chi-Squared with df = 1
    if (mu <= 0.0) {
      throw 
	invalid_argument("The Inverse Gaussian's mean is non-positive.\n");  
    }
    else if (lambda <= 0.0) {
      throw 
	invalid_argument("The Inverse Gaussian's scale is non-positive.\n");
    }
    else {
      double u = unif_rand();
      //Find the smallest root:
      double x = mu + a * v - b * sqrt( c * v + d * pow(v, 2) );
      // Choose x with probability = mean / (mean + x), else choose d/x:
      outVal = ( u < ( mu / (mu + x) ) ) ? x : d / x; 
    }  
  }// CLOSE while(outVal !> 0.0)
  return outVal;
}// END drawInvGauss()
