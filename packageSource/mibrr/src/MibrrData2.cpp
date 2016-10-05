// Title:    Function definitions for the MibrrData2 Class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2016-MAY-15
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

#include "MibrrData2.hpp"

///////////////////////// CONSTRUCTORS / DESTRUCTOR /////////////////////////////
  
MibrrData2::MibrrData2(MatrixXd &newData)
{
  //MibrrData::MibrrData(newData);
  _data     = newData;
  //_dataMeans  = RowVectorXd( _data.cols() );
  _dataScales = RowVectorXd( _data.cols() );
}

MibrrData2::MibrrData2(MatrixXd &newData,
		       //VectorXd &newDataMeans,
		       VectorXd &newDataScales,
		       double missingDataCode) 
{
  //MibrrData::MibrrData(newData,
  //		       newDataMeans,
  //		       newDataScales,
  //		       missingDataCode);
  _missingDataCode   = missingDataCode;
  _data              = newData;
  //_dataMeans         = newDataMeans;
  _dataScales        = newDataScales;
  _nonresponseFilter = ArrayXXb::Constant(_data.rows(), _data.cols(), false);
  computeNonresponseFilter();
}

MibrrData2::~MibrrData2() 
{
}

///////////////////////////////// ACCESSORS /////////////////////////////////////

MatrixXd MibrrData2::getIVs2(int targetIndex)
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


MatrixXd MibrrData2::getFullIVs2(int targetIndex)
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

int MibrrData2::nPreds2() const
{
  return _data.cols() - 1;
}
