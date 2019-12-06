// Title:    Function definitions for the MibrrGibbs class
// Author:   Kyle M. Lang
// Created:  2014-08-24
// Modified: 2019-12-06
// Purpose:  This class contains the Gibbs sampling-related functions for the
//           MIBRR package.

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

#include "MibrrGibbs.h"


//////////////////////// CONSTRUCTORS / DESTRUCTOR /////////////////////////////


MibrrGibbs::MibrrGibbs() 
{
  // _betas, _taus, _nDraws, _l1Parms, and _l2Parms need user-supplied starting
  // values
  _sigma             = 0.0;
  _ridge             = 0.0;
  _lambdas           = VectorXd(2);
  _drawNum           = 0;
  _storeGibbsSamples = false;
  _verbose           = true;
  _penType           = 2;
  _fullBayes         = false;
  _savePpSams        = false;
}

//----------------------------------------------------------------------------//

MibrrGibbs::~MibrrGibbs() {}


//////////////////////////////// ACCESSORS /////////////////////////////////////


VectorXd MibrrGibbs::getBetas()           const { return _betas;               }
ArrayXd  MibrrGibbs::getTaus()            const { return _taus;                }
double   MibrrGibbs::getSigma()           const { return _sigma;               }
MatrixXd MibrrGibbs::getBetaSam()         const { return _betaSam;             }
ArrayXXd MibrrGibbs::getTauSam()          const { return _tauSam;              }
VectorXd MibrrGibbs::getSigmaSam()        const { return _sigmaSam;            }
MatrixXd MibrrGibbs::getImpSam()          const { return _impSam;              }
MatrixXd MibrrGibbs::getPpSam()           const { return _ppSam;               }
MatrixXd MibrrGibbs::getLambdaSam()       const { return _lambdaSam;           }
int      MibrrGibbs::getNDraws()          const { return _nDraws;              }
int      MibrrGibbs::getPenType()         const { return _penType;             }
double   MibrrGibbs::getRidge()           const { return _ridge;               }
bool     MibrrGibbs::getVerbosity()       const { return _verbose;             }
bool     MibrrGibbs::getDoImp()           const { return _doImp;               }

//----------------------------------------------------------------------------//

VectorXd MibrrGibbs::getLambdas() const
{ 
  VectorXd outLam;
  if(_penType == 2) outLam = _lambdas;
  else              outLam = VectorXd::Constant(1, _lambdas[0]);
  return outLam; 
}

//----------------------------------------------------------------------------//

double MibrrGibbs::getLambdas(int index) const 
{ 
  return _lambdas[index - 1]; 
}


//////////////////////////////// MUTATORS //////////////////////////////////////


void MibrrGibbs::setBetas      (VectorXd &betas)   { _betas       = betas;     }
void MibrrGibbs::setTaus       (ArrayXd &taus)     { _taus        = taus;      }
void MibrrGibbs::setSigma      (double sigma)      { _sigma       = sigma;     }
void MibrrGibbs::setTargetIndex(int index)         { _targetIndex = index;     }
void MibrrGibbs::setNDraws     (int nDraws)        { _nDraws      = nDraws;    }
void MibrrGibbs::setDoImp      (bool doImp)        { _doImp       = doImp;     }
void MibrrGibbs::beQuiet       ()                  { _verbose     = false;     }
void MibrrGibbs::doFullBayes   ()                  { _fullBayes   = true;      }
void MibrrGibbs::savePpSams    ()                  { _savePpSams  = true;      }
void MibrrGibbs::setPenType    (int penType)       { _penType     = penType;   }
void MibrrGibbs::setRidge      (double ridge)      { _ridge       = ridge;     }
void MibrrGibbs::setLam1Parms  (VectorXd& l1Parms) { _l1Parms     = l1Parms;   }
void MibrrGibbs::setLam2Parms  (VectorXd& l2Parms) { _l2Parms     = l2Parms;   }
void MibrrGibbs::setFinalRep   (bool finalRep)     { _finalRep    = finalRep;  }
void MibrrGibbs::setLambdas    (VectorXd& lambdas) { _lambdas     = lambdas;   }

//----------------------------------------------------------------------------//

void MibrrGibbs::setLambdas(double lambda1, double lambda2) 
{
  _lambdas[0] = lambda1;
  _lambdas[1] = lambda2;
}

//----------------------------------------------------------------------------//

void MibrrGibbs::setLambdas(double lambda) 
{
  _lambdas[0] = lambda;
  _lambdas[1] = 0.0;
}

//----------------------------------------------------------------------------//

void MibrrGibbs::setLambdaParms(VectorXd& lambdaParms)
{
  _l1Parms = lambdaParms.head(2);
  _l2Parms = lambdaParms.tail(2);
}

//----------------------------------------------------------------------------//

void MibrrGibbs::startParameters(VectorXd &betaStarts,
				 ArrayXd  &tauStarts,
				 double   sigmaStart,
				 double   lambda1Start,
				 double   lambda2Start,
				 bool     useBetaMeans)
{
  _betas      = betaStarts;
  _taus       = tauStarts;
  _sigma      = sigmaStart;
  _lambdas[0] = lambda1Start;
  _lambdas[1] = lambda2Start;

  _useBetaMeans = useBetaMeans;
  _betaMeans    = VectorXd::Zero(_betas.size());
}

//----------------------------------------------------------------------------//

void MibrrGibbs::startParameters(VectorXd &betaStarts,
				 ArrayXd  &tauStarts,
				 double   sigmaStart,
				 VectorXd &lambdaStartVec,
				 bool     useBetaMeans)
{
  _betas   = betaStarts;
  _taus    = tauStarts;
  _sigma   = sigmaStart;
  _lambdas = lambdaStartVec;

  _useBetaMeans = useBetaMeans;
  _betaMeans    = VectorXd::Zero(_betas.size());
}

//----------------------------------------------------------------------------//

void MibrrGibbs::startGibbsSampling(const MibrrData &mibrrData)
{
  if(_doImp) _impSam = MatrixXd(_nDraws, mibrrData.nMiss(_targetIndex));
  else       _impSam = MatrixXd::Zero(1, 1);
  
  if(_savePpSams) _ppSam = MatrixXd(_nDraws, mibrrData.nResp(_targetIndex));
  else            _ppSam = MatrixXd::Zero(1, 1);
  
  if(_fullBayes) _lambdaSam = MatrixXd(_nDraws, 2);
  else           _lambdaSam = MatrixXd::Zero(1, 1);
  
  _betaSam  = MatrixXd::Zero(_nDraws, _betas.size()); /////////////////////////////////////////////
  _tauSam   = MatrixXd::Zero(_nDraws, _taus.size()); ///////////////////////////////////////////////
  _sigmaSam = VectorXd::Zero(_nDraws); ////////////////////////////////////////////////////////////
  
  _storeGibbsSamples = true;
}


////////////////////////// PARAMETER UPDATE FUNCTIONS //////////////////////////


void MibrrGibbs::updateLambdas(const MibrrData &mibrrData)
{
  int    nPreds = mibrrData.nPreds();
  double lam1   = _lambdas[0];
  double tauSum = _taus.sum();
  
  if(_penType == 1) {// MIBL Version
    double shape = _l1Parms[0] + nPreds;
    double rate  = _l1Parms[1] + tauSum / 2.0;
    
    _lambdas[0] = sqrt(drawGamma(shape, rate));
  }
  else {             // MIBEN Version
    double lam2 = _lambdas[1];
    
    // Sample new value of lambda1:
    double shape = _l1Parms[0] + nPreds / 2.0;
    double rate  = _l1Parms[1] + tauSum / (8.0 * lam2 * _sigma);
    
    _lambdas[0] = sqrt(drawGamma(shape, rate));
    
    // Sample new value of lambda2:
    double chi = _l2Parms[0] + (pow(lam1, 2) * tauSum) / (4 * _sigma);
    double psi =
      _l2Parms[1] + ((tauSum / (_taus - 1.0).sum()) *
		     _betas.tail(nPreds).array().square().sum()) / _sigma;
    
    _lambdas[1] = drawGig(1.0, chi, psi);
  }

  // Add the updated lambdas to their Gibbs sample:
  if(_storeGibbsSamples) _lambdaSam.row(_drawNum) = _lambdas.transpose();
}// END updateLambdas()

//----------------------------------------------------------------------------//

void MibrrGibbs::updateTaus(const MibrrData &mibrrData)
{
  int     nPreds   = mibrrData.nPreds();
  double  tauScale = -1.0; // -1 to ensure an exception if try() fails
  ArrayXd tauMeans;
  
  try {
    if(_penType == 2) {// MIBEN Version
      ArrayXd tauMeansNumerator = ArrayXd::Constant(nPreds, _lambdas[0]);
      
      tauMeans = tauMeansNumerator.sqrt() /
	(2.0 * _lambdas[1] * _betas.tail(nPreds).array().abs());
      
      tauScale = _lambdas[0] / (4.0 * _lambdas[1] * _sigma);
    }
    else {             // MIBL Version
      double tauMeansNumerator = pow(_lambdas[0], 2) * _sigma;
      
      tauMeans =
	(tauMeansNumerator / _betas.tail(nPreds).array().square()).sqrt();
      
      tauScale = pow(_lambdas[0], 2);
    }
    if((tauMeans <= 0.0).any()) throw 1;
    if(tauScale <= 0.0)         throw 2;
  }
  catch(int e) { tauError(e); }
  
  // Draw new values of the auxiliary penalty parameters:
  ArrayXd tmpDraws(nPreds);
  for(int i = 0; i < nPreds; i++)
    tmpDraws[i] = drawInvGauss(tauMeans[i], tauScale);
  
  if(_penType == 2) _taus = (tmpDraws + 1.0) / tmpDraws; // MIBEN Version   
  else              _taus = 1.0 / tmpDraws;              // MIBL Version
  
  // Add the updated Taus to their Gibbs sample:
  if(_storeGibbsSamples) _tauSam.row(_drawNum) = _taus.transpose();
}// END updateTaus()

//----------------------------------------------------------------------------//

void MibrrGibbs::updateBetas(MibrrData &mibrrData)
{
  int             nPreds = mibrrData.nPreds();
  int             nObs   = mibrrData.nResp(_targetIndex);
  LDLT <MatrixXd> aMatChol;
    
  // Compute the IV's crossproducts matrix:
  MatrixXd aMat = mibrrData.getIVs(_targetIndex, true, true).transpose() *
    mibrrData.getIVs(_targetIndex, true, true);
  
  // Compute an appropriate penalty term:
  MatrixXd penalty;
  if(_penType == 2) {     // MIBEN Version
    VectorXd transTaus = _taus / (_taus - 1.0);
    penalty = _lambdas[1] * transTaus.asDiagonal();    
  }
  else if(_penType == 1) {// MIBL Version
    VectorXd transTaus = 1.0 / _taus;
    penalty = transTaus.asDiagonal();
  }
  else                    // Basic ridge penalty
    penalty = (aMat.diagonal() * _ridge).asDiagonal();
  
  // Penalize the IV's crossproducts matrix:
  aMat += penalty;
  
  VectorXd betaMeans;
  MatrixXd betaCov;
  try {
    // Compute the Cholesky decomposition of the IV's crossproducts, and use it
    // to compute the moments of beta's fully conditional posterior:
    aMatChol.compute(aMat);
    _betaMeans.tail(nPreds) = 
      aMatChol.solve(mibrrData.getIVs(_targetIndex, true, true).transpose() *
		     mibrrData.getDV(_targetIndex)); 
    betaCov = aMatChol.solve(_sigma * MatrixXd::Identity(nPreds, nPreds));
  }
  catch(exception &e) { betaError(e); }
  
  // Draw new values of the regression slope coefficients:
  _betas.tail(nPreds) = drawMvn(_betaMeans.tail(nPreds), betaCov);
  
  // Compute the mean of the intercept's distribution:
  _betaMeans[0] = mibrrData.getDV(_targetIndex).mean();
  
  // Draw a new value of the intercept term:  
  _betas[0] = drawNorm(_betaMeans[0], sqrt(_sigma / double(nObs)));
  
  if(_storeGibbsSamples) {
    VectorXd betas = _betas;

    // Back-transform the standardized betas to their raw metric:
    if(_finalRep) {
      betas.tail(nPreds).array() /= mibrrData.getScales(_targetIndex).array();
      betas[0] -= mibrrData.getMeans(_targetIndex) * betas.tail(nPreds);
    }
    
    _betaSam.row(_drawNum) = betas.transpose();
  }
}// END updateBetas()

//----------------------------------------------------------------------------//

void MibrrGibbs::updateSigma(MibrrData &mibrrData)
{
  double par1, par2; // parameters of sigma's posterior distribution
  int    nObs   = mibrrData.nResp(_targetIndex);
  int    nPreds = mibrrData.nPreds();
  double nBetas = double(nPreds) + 1.0;

  VectorXd betas;
  if((_penType == 0) | _useBetaMeans) betas = _betaMeans;
  else                                betas = _betas;
  
  // Compute the OLS residuals:
  VectorXd resid = mibrrData.getDV(_targetIndex) -
    (betas[0] * VectorXd::Ones(nObs) +
     mibrrData.getIVs(_targetIndex, true, true) * betas.tail(nPreds)
     );
  
  // Compute the residual sum of squares:
  double sse = resid.squaredNorm();
  
  if(_penType == 2) {     // MIBEN Version
    double scaleSum =
      (_taus / (_taus - 1.0) * betas.tail(nPreds).array().square()).sum();
 
    par1 = (double(nObs) / 2.0) + nBetas; 
    par2 = 0.5 * (sse + _lambdas[1] * scaleSum +
		  (pow(_lambdas[0], 2) / (4.0 * _lambdas[1])) * _taus.sum());
 
    bool   validDraw = false;
    double testDraw;
    while(!validDraw) {// Rejection sampling to draw a Sigma variate
      testDraw = drawInvGamma(par1, par2);
      
      double igShape = pow(_lambdas[0], 2) / (8.0 * testDraw * _lambdas[1]);
      double igVal   = calcIncGamma(0.5, igShape, false);
      double thresh  = _unif(_gen);
      
      validDraw =
	log(thresh) <= (nBetas * log(tgamma(0.5))) - (nBetas * log(igVal));
      
      Rcpp::checkUserInterrupt();
    };
    _sigma = testDraw;
  }
  else if(_penType == 1) {// MIBL Version
    VectorXd transTaus = 1.0 / _taus;
    double   penalty   = betas.tail(nPreds).transpose() *
      transTaus.asDiagonal() * betas.tail(nPreds);
    
    par1 = 0.5 * ((double(nObs) - 1.0) + nBetas);
    par2 = 0.5 * (sse + penalty);
    
    // Draw a new sigma variate:
    _sigma = drawInvGamma(par1, par2);
  }
  else {                  // Basic Ridge Version
    par1 = double(nObs) - nBetas;
    par2 = sse / par1;
    
    // Draw a new sigma variate:
    _sigma = drawInvChiSq(par1, par2);
  }
  
  // Add the updated sigma to its Gibbs sample:
  if(_storeGibbsSamples) _sigmaSam[_drawNum] = _sigma;
}// END updateSigma()

//----------------------------------------------------------------------------//

void MibrrGibbs::updateImputations(MibrrData &mibrrData, bool ppSam)
{
  int nPreds = mibrrData.nPreds();
  
  int size;
  if(ppSam) size = mibrrData.nResp(_targetIndex);
  else      size = mibrrData.nMiss(_targetIndex);
  
  // Draw a vector of imputations from the posterior predictive distribution
  // of the missing data:
  VectorXd imps = _betas[0] * VectorXd::Ones(size) +
    mibrrData.getIVs(_targetIndex, ppSam, true) * _betas.tail(nPreds) +
    drawNorm(size, 0.0, sqrt(_sigma));
  
  if(ppSam) {
    // Add the updated posterior predictive sample to its Gibbs sample:
    _ppSam.row(_drawNum) = imps.transpose();
  }
  else {
    // Replace the missing data in the target variable with the imputations:
    mibrrData.fillMissing(imps, _targetIndex);
    
    // Add the updated imputations to their Gibbs sample:
    if(_storeGibbsSamples) _impSam.row(_drawNum) = imps.transpose();
  }
}// END updateImputations()

//----------------------------------------------------------------------------//

void MibrrGibbs::doGibbsIteration(MibrrData &mibrrData)
{
  if(_fullBayes & (_penType != 0)) updateLambdas(mibrrData);
  
  if(_penType != 0) updateTaus(mibrrData);
  
  updateBetas(mibrrData);
  updateSigma(mibrrData);
  
  if(_doImp) updateImputations(mibrrData, false);
  
  if(_storeGibbsSamples) {
    // Save a posterior predictive sample of Y_obs:
    if(_savePpSams) updateImputations(mibrrData, true);
    _drawNum++;
  }
}// END doGibbsIteration()


///////////////////////// EXCEPTION HANDLING FUNCTIONS /////////////////////////


void MibrrGibbs::tauError(int errorCode) const
{
  // Print the most recent values of each parameter:
  dumpParameters();
  
  if(errorCode == 1) {
    Rcpp::Rcout << "\n";
    Rcpp::stop("Ouch! My tau is broken :(\nSomething terrible has occured \
while updating Tau,\nand one of its mean values is non-positive.\n");
  }
  else if (errorCode == 2) {
    Rcpp::Rcout << "\n";
    Rcpp::stop("Ouch! My tau is broken :(\nSomething terrible has occured \
while updating Tau,\nand one of its scale values is non-positive.\n");
  }
}

//----------------------------------------------------------------------------//

void MibrrGibbs::betaError(exception &e) const
{
  Rcpp::Rcerr << e.what() << endl;
  Rcpp::stop("Something terrible has occured while updating Beta.\nAbove this \
message, I've printed the that exception I caught.\nBeta luck next time.");
}

//----------------------------------------------------------------------------//

void MibrrGibbs::dumpParameters() const
{
  dbg(_targetIndex);
  dbg(_penType);
  dbg(_sigma);
  dbg(_betas);
  dbg(_taus);
  dbg(_lambdas);
}
