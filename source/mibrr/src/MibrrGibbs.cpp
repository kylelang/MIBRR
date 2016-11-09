// Title:    Function definitions for the MibrrGibbs class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2016-NOV-08
// Purpose:  This class contains the Gibbs sampling-related functions for the
//           MIBRR package.

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

#include "MibrrGibbs.hpp"

//////////////////////// CONSTRUCTORS / DESTRUCTOR //////////////////////////////
  
MibrrGibbs::MibrrGibbs() 
{
  // _beta, _tau, _nEmIters, and _nDraws need user-supplied starting values
  _sigma             = 0.0;
  _lambdas           = VectorXd(2);
  _emConvTol         = 1.0e-5;
  _lambdaWindow      = 1;
  _drawNum           = 0;
  _emIterNum         = 0;
  _optIterCount      = 0;
  _storeGibbsSamples = false;
  _verbose           = true;
  _useElasticNet     = true;
  _doImputation      = true;
  _twoPhaseOpt       = true;
}


MibrrGibbs::~MibrrGibbs() {}

//////////////////////////////// ACCESSORS //////////////////////////////////////

VectorXd MibrrGibbs::getBetas()           const { return _betas;                }
ArrayXd  MibrrGibbs::getTaus()            const { return _taus;                 }
double   MibrrGibbs::getSigma()           const { return _sigma;                }
MatrixXd MibrrGibbs::getBetaSam()         const { return _betaSam;              }
ArrayXXd MibrrGibbs::getTauSam()          const { return _tauSam;               }
VectorXd MibrrGibbs::getSigmaSam()        const { return _sigmaSam;             }
MatrixXd MibrrGibbs::getImpSam()          const { return _impSam;               }
MatrixXd MibrrGibbs::getLambdaHistory()   const { return _lambdaHistory;        }
int      MibrrGibbs::getNDraws()          const { return _nDraws;               }
int      MibrrGibbs::getNEmIters()        const { return _nEmIters;             }
bool     MibrrGibbs::getVerbosity()       const { return _verbose;              }
bool     MibrrGibbs::getElasticNetFlag()  const { return _useElasticNet;        }
bool     MibrrGibbs::getDoImputation()    const { return _doImputation;         }
bool     MibrrGibbs::getSimpleIntercept() const { return _simpleIntercept;      }

VectorXd MibrrGibbs::getLambdas() const
{ 
  VectorXd outLam;
  if(_useElasticNet) outLam = _lambdas;
  else               outLam = VectorXd::Constant(1, _lambdas[0]);
  return outLam; 
}


double MibrrGibbs::getLambdas(int index) const 
{ 
  return _lambdas[index - 1]; 
}


//////////////////////////////// MUTATORS ///////////////////////////////////////


void MibrrGibbs::setBetas      (VectorXd &betas)   { _betas = betas;            }
void MibrrGibbs::setTaus       (ArrayXd &taus)     { _taus = taus;              }
void MibrrGibbs::setSigma      (double sigma)      { _sigma = sigma;            }
void MibrrGibbs::setNEmIters   (int nEmIters)      { _nEmIters = nEmIters;      }
void MibrrGibbs::setTargetIndex(int index)         { _targetIndex = index;      }
void MibrrGibbs::setNDraws     (int nDraws)        { _nDraws = nDraws;          }
void MibrrGibbs::beQuiet       ()                  { _verbose = false;          }
void MibrrGibbs::doBl          ()                  { _useElasticNet = false;    }
void MibrrGibbs::doPrediction  ()                  { _doImputation = false;     }
void MibrrGibbs::useSimpleInt  ()                  { _simpleIntercept = true;   }
void MibrrGibbs::setLambdas    (VectorXd& lambdas) { _lambdas = lambdas;        }


void MibrrGibbs::setLambdas(double lambda1, double lambda2) 
{
  _lambdas[0] = lambda1;
  _lambdas[1] = lambda2;
}


void MibrrGibbs::setLambdas(double lambda) 
{
  _lambdas[0] = lambda;
  _lambdas[1] = 0.0;
}


void MibrrGibbs::setLambdas() 
{  
  int      startRow      = _emIterNum - _lambdaWindow;
  int      nCols         = _useElasticNet ? 2 : 1;
  VectorXd pooledLambdas = 
    _lambdaHistory.block(startRow, 0, _lambdaWindow, nCols).colwise().mean();

  if(_useElasticNet) _lambdas = pooledLambdas.transpose();
  else               _lambdas[0] = pooledLambdas[0];
}


void MibrrGibbs::startParameters(VectorXd &betaStarts,
				 ArrayXd  &tauStarts,
				 double   sigmaStart,
				 double   lambda1Start,
				 double   lambda2Start)
{
  _betas         = betaStarts;
  _taus          = tauStarts;
  _sigma         = sigmaStart;
  _lambdas[0]    = lambda1Start;
  _lambdas[1]    = lambda2Start;
  int lamNum     = _useElasticNet ? 2 : 1;
  _lambdaHistory = MatrixXd(_nEmIters, lamNum);
}


void MibrrGibbs::startParameters(VectorXd &betaStarts,
				 ArrayXd  &tauStarts,
				 double   sigmaStart,
				 VectorXd &lambdaStartVec)
{
  _betas         = betaStarts;
  _taus          = tauStarts;
  _sigma         = sigmaStart;
  _lambdas       = lambdaStartVec;
  int lamNum     = _useElasticNet ? 2 : 1;
  _lambdaHistory = MatrixXd(_nEmIters, lamNum);
}


void MibrrGibbs::setupOptimizer(int    nEmIters,
				int    lambdaWindow,
				double emConvTol,
				bool   twoPhaseOpt)
{
  _emConvTol    = emConvTol;
  _nEmIters     = nEmIters;
  _lambdaWindow = lambdaWindow;
  _twoPhaseOpt  = twoPhaseOpt;
}


void MibrrGibbs::startGibbsSampling(const MibrrData &mibrrData)
{
  int nMiss          = mibrrData.nMiss(_targetIndex);
  _storeGibbsSamples = true;
  
  _betaSam  = MatrixXd(_nDraws, _betas.size());
  _tauSam   = MatrixXd(_nDraws, _taus.size());
  _sigmaSam = VectorXd(_nDraws);
  _impSam   = MatrixXd(_nDraws, nMiss);
}


void MibrrGibbs::restartParameters(MibrrData &mibrrData)
{
  _betas = _betaSam.bottomRows(_betaSam.rows() - 1).colwise().mean().transpose();
  _taus  = _tauSam.bottomRows(_tauSam.rows() - 1).colwise().mean().transpose();
  _sigma = _sigmaSam.tail(_sigmaSam.size() - 1).mean();
  
  VectorXd meanImps = _impSam.bottomRows(_impSam.rows() - 1).colwise().mean();
  mibrrData.fillMissing(meanImps, _targetIndex);
}


/////////////////////////// RANDOM VARIATE SAMPLERS /////////////////////////////


double MibrrGibbs::drawInvGamma(double shape, double scale) const
{
  return 1.0 / R::rgamma(shape, 1.0 / scale);
}//END drawInvGamma()


double MibrrGibbs::calcIncGamma(const double shape, 
				const double cutVal,
				const bool   lowerTail)
{
  double scale   = 1.0;
  int    lower   = (int)lowerTail;
  int    logTran = 0; // Don't want log transform
  
  return R::pgamma(cutVal, shape, scale, lower, logTran) * tgamma(shape);
}// END calcIncGamma()


double MibrrGibbs::drawInvGauss(const double mu, const double lambda)
{ 
  double b      = 0.5 * mu / lambda;
  double a      = mu * b;
  double c      = 4.0 * mu * lambda;
  double d      = pow(mu, 2);
  double outVal = 0.0;

  while(outVal <= 0.0) {
    double tmpDraw = norm_rand();
    double         v = pow(tmpDraw, 2); // Chi-Squared with df = 1
    if (mu <= 0.0) {
      throw invalid_argument("The Inverse Gaussian's mean is non-positive.\n");  
    }
    else if (lambda <= 0.0) { 
      throw invalid_argument("The Inverse Gaussian's scale is non-positive.\n");
    }
    else {
      double u = unif_rand();
      //Find the smallest root:
      double x = mu + a * v - b * sqrt(c * v + d * pow(v, 2));
      // Choose x with probability = mean / (mean + x), else choose d/x:
      outVal = ( u < ( mu / (mu + x) ) ) ? x : d / x; 
    }  
  }// CLOSE while(outVal !> 0.0)
  return outVal;
}// END drawInvGauss()


////////////////////////// PARAMETER UPDATE FUNCTIONS ///////////////////////////


void MibrrGibbs::updateTaus(const MibrrData &mibrrData)
{
  int     nPreds   = mibrrData.nPreds();
  double  tauScale = -1.0; // -1 to ensure an exception if try() fails
  ArrayXd tauMeans;
  
  try {
    if(_useElasticNet) {// Miben Version
      ArrayXd tauMeansNumerator = ArrayXd::Constant(nPreds, _lambdas[0]);
      
      tauMeans = tauMeansNumerator.sqrt() /
	(2.0 * _lambdas[1] * _betas.tail(nPreds).array().abs());
      
      tauScale = _lambdas[0] / (4.0 * _lambdas[1] * _sigma);
    }
    else {// MIBL Version
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
  
  ArrayXd newTaus;
  if(_useElasticNet) newTaus = (tmpDraws + 1.0) / tmpDraws; // MIBEN Version   
  else               newTaus = 1.0 / tmpDraws;              // MIBL Version
  
  _taus = newTaus; // Store the updated Taus
  
  // Add the updated Taus to their Gibbs sample:
  if(_storeGibbsSamples) _tauSam.row(_drawNum) = newTaus.transpose();
}// END updateTaus()



void MibrrGibbs::updateBetas(const MibrrData &mibrrData)
{
  int             nPreds = mibrrData.nPreds();
  int             nObs   = mibrrData.nResp(_targetIndex);
  MatrixXd        tmpMat;
  LDLT <MatrixXd> aMatrixCholesky;
  
  if(_useElasticNet) {// MIBEN Version
    VectorXd transformedTaus = _taus / (_taus - 1.0);
    tmpMat = _lambdas[1] * transformedTaus.asDiagonal();    
  }
  else {// MIBL Version
    VectorXd transformedTaus = 1.0 / _taus;
    tmpMat = transformedTaus.asDiagonal();
  }
  
  MatrixXd aMatrix;
  aMatrix = mibrrData.getObsIVs(_targetIndex).transpose() *
    mibrrData.getObsIVs(_targetIndex) + tmpMat;
  
  VectorXd betaMeans;
  MatrixXd betaCovariances;
  try {
    // Compute the Cholesky decomposition of the IV's penalized crossproducts,
    // and use it to compute the moments of beta's fully conditional posterior:
    aMatrixCholesky.compute(aMatrix);
    betaMeans =
      aMatrixCholesky.solve(mibrrData.getObsIVs(_targetIndex).transpose() *
			    mibrrData.getObsDV(_targetIndex)); 
    
    tmpMat = _sigma * MatrixXd::Identity(nPreds, nPreds);   
    betaCovariances = aMatrixCholesky.solve(tmpMat);
  }
  catch(exception &e) { betaError(e); }

  // Draw a new value of the intercept term:
  double intSd = sqrt(_sigma / double(nObs));
  double intMean;
  if(_simpleIntercept)
    intMean = mibrrData.getObsDV(_targetIndex).mean();
  else
    intMean = mibrrData.getObsDV(_targetIndex).mean() -
      mibrrData.getObsIVs(_targetIndex).colwise().mean() * _betas.tail(nPreds);
  
  VectorXd newBetas(nPreds + 1);  
  newBetas[0] = R::rnorm(intMean, intSd);
  
  // Draw new values of the regression slope coefficients:
  newBetas.tail(nPreds) = mibrrData.drawMVN(betaMeans, betaCovariances);
  
  _betas = newBetas; // Store the updated Betas
  
  // Add the updated Betas to their Gibbs sample:
  if(_storeGibbsSamples) _betaSam.row(_drawNum) = newBetas.transpose();
}// END updateBetas()
    


void MibrrGibbs::updateSigma(const MibrrData &mibrrData)
{
  double   newSigma;
  int      nPreds        = mibrrData.nPreds();
  int      nObs          = mibrrData.nResp(_targetIndex);
  VectorXd tmpBiasVector = VectorXd::Ones(nObs);
  
  // Compute the regularized residual sum of squares:
  double sse =
    (mibrrData.getObsDV(_targetIndex) - _betas[0] * tmpBiasVector -
     mibrrData.getObsIVs(_targetIndex) * _betas.tail(nPreds)).transpose() *
    (mibrrData.getObsDV(_targetIndex) - _betas[0] * tmpBiasVector -
     mibrrData.getObsIVs(_targetIndex) * _betas.tail(nPreds));
  
  if(_useElasticNet) {// MIBEN Version
    double scaleSum =
      (_taus / (_taus - 1.0) * _betas.tail(nPreds).array().square()).sum();
    
    double sigmaShape = (double(nObs) / 2.0) + double(nPreds);
    double sigmaScale =
      0.5 * (sse + _lambdas[1] * scaleSum +
	     (pow(_lambdas[0], 2) / (4.0 * _lambdas[1])) * _taus.sum());
    
    bool   isDrawValid = false;
    double testDraw;
    while(!isDrawValid) {// Rejection sampling to draw a Sigma variate
      testDraw         = drawInvGamma(sigmaShape, sigmaScale);
      double threshold = unif_rand();
      double igShape   = pow(_lambdas[0], 2) / (8.0 * testDraw * _lambdas[1]);
      double igDraw    = calcIncGamma(0.5, igShape, false);
      isDrawValid      = log(threshold) <=
	(double(nPreds) * log(tgamma(0.5))) - (double(nPreds) * log(igDraw));
      Rcpp::checkUserInterrupt();
    };
    newSigma = testDraw;
  }
  else {// MIBL Version
    VectorXd transformedTaus = 1.0 / _taus;
    MatrixXd tmpMat = transformedTaus.asDiagonal();
    
    double penaltyTerm =
      _betas.tail(nPreds).transpose() * tmpMat * _betas.tail(nPreds);
    double sigmaShape = 0.5 * ((double(nObs) - 1.0) + double(nPreds));
    double sigmaScale = 0.5 * (sse + penaltyTerm);
    
    // Draw a new sigma variate:
    newSigma = drawInvGamma(sigmaShape, sigmaScale);
  }
  
  _sigma = newSigma; // Store the updated Sigma
  
  // Add the updated sigma to its Gibbs sample:
  if(_storeGibbsSamples) _sigmaSam[_drawNum] = newSigma;
}// END updateSigma()



void MibrrGibbs::updateImputations(MibrrData &mibrrData)
{
  int         nMiss         = mibrrData.nMiss(_targetIndex);
  int         nPreds        = mibrrData.nPreds();
  VectorXd    tmpBiasVector = VectorXd::Ones(nMiss);
  VectorXd    errorVector(nMiss);
  
  // Draw the residual error terms for the imputation model:
  for(int i = 0; i < nMiss; i++) errorVector[i] = R::rnorm(0.0, sqrt(_sigma));
  
  // Draw a vector of imputations from the posterior predictive distribution
  // of the missing data:
  VectorXd tmpImps = _betas[0] * tmpBiasVector +
    mibrrData.getMissIVs(_targetIndex) * _betas.tail(nPreds) + errorVector;
  
  // Replace the missing data in the target variable with the imputations:
  mibrrData.fillMissing(tmpImps, _targetIndex);
  
  // Add the updated imputations to their Gibbs sample:
  if(_storeGibbsSamples) _impSam.row(_drawNum) = tmpImps.transpose();
}// END updateImputations()



void MibrrGibbs::doGibbsIteration(MibrrData &mibrrData)
{  
  updateTaus       (mibrrData);
  updateBetas      (mibrrData);
  updateSigma      (mibrrData);
  updateImputations(mibrrData);
  
  if(_storeGibbsSamples) _drawNum++;
}// END doGibbsIteration()



double MibrrGibbs::lambdaObjective(const std::vector<double> &lambdas,
				   std::vector<double>       &grad,
				   void                      *extras)
{
  // Check for finite, well-defined lambda values:
  bool nanLam1 = std::isnan(lambdas[0]) | lambdas[0] != lambdas[0];
  bool nanLam2 = std::isnan(lambdas[1]) | lambdas[1] != lambdas[1];
  
  if(nanLam1 & nanLam2) 
    throw invalid_argument("Both elastic net penalty parameters are NaN.");  
  else if(nanLam1) 
    throw invalid_argument("The LASSO penalty parameter is NaN.");  
  else if(nanLam2) 
    throw invalid_argument("The ridge penalty parameter is NaN.");  
  
  _optIterCount++; // Track the number of function evaluations
  
  int nPreds = _tauSam.cols();
  int nSams  = _tauSam.rows();

  // BEGIN compute LL and gradient terms.
  ArrayXd tmpArray = 2.0 * _sigmaSam;
  double  w1       = nPreds * log(lambdas[0]);
  ArrayXd w2       = lambdas[1] / tmpArray;
  ArrayXd w3       = 1 / tmpArray;
  
  tmpArray   = 8.0 * _sigmaSam * lambdas[1];
  ArrayXd w4 = pow(lambdas[0], 2) / tmpArray;
  double  w5 = nPreds / lambdas[0];
  double  w6 = (nPreds * lambdas[0]) / (4.0 * lambdas[1]);
  
  tmpArray   = 4.0 * _sigmaSam * lambdas[1];
  ArrayXd w7 = lambdas[0] / tmpArray;
  double  w8 = (nPreds * pow(lambdas[0], 2)) / (8.0 * pow(lambdas[1], 2));
  
  tmpArray   = 8.0 * _sigmaSam * pow(lambdas[1], 2);
  ArrayXd w9 = pow(lambdas[0], 2) / tmpArray;
  
  ArrayXd igArray(nSams); 
  for(int i = 0; i < nSams; i++)
    igArray[i] = calcIncGamma(0.5, w4[i], false);
  
  ArrayXd term1 = igArray.log();
  ArrayXd term2 =
    ((_tauSam / (_tauSam - 1)) *
     _betaSam.rightCols(nPreds).array().square()).rowwise().sum();
  // END compute LL and gradient terms.
  
  if(!grad.empty()) {// Calculate the gradient vector:
    ArrayXd term3 = (1 / igArray) * (1 / w4.sqrt()) * (-1.0 * w4).exp() *
      (1 / _sigmaSam.array());
    grad[0] = (w5 + (w6 * term3) - w7 * _tauSam.rowwise().sum()).mean();
    grad[1] =
      (-1.0 * (w8 * term3) - (w3 * term2) + w9 * _tauSam.rowwise().sum()).mean();
  }
  
  // Return the objective value:
  return (w1 - (nPreds * term1) - (w2 * term2) -
	  (w4 * _tauSam.rowwise().sum())).mean();
}// END lambdaObjective()



// As suggested by the nlopt authors, specify a simple  wrapper function for
// lambdaObjective() so that nlopt will run inside the MibrrGibbs class:
double objectiveWrap(const std::vector<double> &lambdas, 
		     std::vector<double>       &grad, 
		     void                      *extras) 
{
  MibrrGibbs *obj = static_cast<MibrrGibbs *>(extras);   
  return obj -> lambdaObjective(lambdas, grad, extras);
}// END objectiveWrap()



void MibrrGibbs::optimizeMibenLambdas(const bool preOptimize) 
{
  string              algName, outPrefix; 
  bool                optimized = false;
  std::vector<double> lambdas(2), grad(2), lamBounds(2);
  nlopt::algorithm    optAlg;
  
  lambdas[0]   = _lambdas[0];
  lambdas[1]   = _lambdas[1];
  grad[0]      = 0.0;
  grad[1]      = 0.0;
  lamBounds[0] = 1.0e-4; // Set low bounds slightly above zero to avoid dividing 
  lamBounds[1] = 1.0e-4; // by zero when no regularization is needed.
  _optMethod   = 0;
  
  while(!optimized) {//Try different algorithms until convergence
    if(preOptimize) {
      _optPrefix = "pre-";
      if(_optMethod == 0) {
	optAlg   = nlopt::LN_BOBYQA;
	_algName = "BOBYQA";
      }
      else if(_optMethod == 1) {
	optAlg   = nlopt::LN_COBYLA;
	_algName = "COBYLA";
      }
      else if(_optMethod == 2) {
	optAlg   = nlopt::LN_SBPLX;
	_algName = "SBPLX";
      }
      else if(_optMethod == 3) {
	optAlg   = nlopt::LN_PRAXIS;
	_algName = "PRAXIS";
      }
    }
    else {
      _optPrefix = "";
      if(_optMethod == 0) {
	optAlg   = nlopt::LD_MMA;
	_algName = "MMA";
      }
      else if(_optMethod == 1) {
	optAlg   = nlopt::LD_VAR2;
	_algName = "VAR2";
      }
      else if(_optMethod == 2) {
	optAlg   = nlopt::LD_VAR1;
	_algName = "VAR1";
      }
      else if(_optMethod == 3) {
	optAlg   = nlopt::LD_LBFGS;
	_algName = "LBFGS";
      }
    }// END if(preOptimize}
    
    nlopt::opt myOptimizer(optAlg, 2);       // Initialize the optimizer object
    myOptimizer.set_max_objective(objectiveWrap, this);
    myOptimizer.set_ftol_rel(_emConvTol);    // Stopping criterion
    myOptimizer.set_lower_bounds(lamBounds); // Force positive Lambdas	
    
    double        maxLL = 0.0;
    nlopt::result myResult;
    try {
      myResult = myOptimizer.optimize(lambdas, maxLL);  
      if(myResult < 0) {          // Catch nlopt failure codes
	_optMethod ++;            // Try the next algorithm
	lambdas[0] = _lambdas[0]; // Reset the lambdas
	lambdas[1] = _lambdas[1];
	lambdaError();
      }
      else if(myResult > 0) {// Successful convergence!
	if(_verbose) {
	  cout << "Lambdas " << _optPrefix << "optimized with " << _algName;
	  cout << " in " << _optIterCount << " iterations" << endl;
	}	
	optimized = true;
      } 
    }
    catch(int &e) {
      throw e;
    }
    catch(exception &e) {       // Catch Lambdas == NaN
      _optMethod ++;            // Try the next algorithm
      lambdas[0] = _lambdas[0]; // Reset the lambdas
      lambdas[1] = _lambdas[1]; 
      lambdaError(e);
    }  
    _optIterCount = 0;
  }// END while(!optimized)      
  // Store the updated penalty parameters:
  _lambdas[0]                    = lambdas[0];
  _lambdas[1]                    = lambdas[1];
  _lambdaHistory.row(_emIterNum) = _lambdas.transpose();
}// END optimizeMibenLambdas()



void MibrrGibbs::updateLambdas()
{
  if(_useElasticNet) {// MIBEN version
    // For MIBEN, optimization can be done in two stages:
    // 1) A rough, gradient-free, pre-optimization is done to move the estimates
    //    into the neighborhood of the MLE.
    // 2) Gradient-based optimization is used to fine-tune the estimates from
    //    Step (1).
    if(_twoPhaseOpt) optimizeMibenLambdas(true);  // Pre-optimization
    optimizeMibenLambdas(false);                  // Optimization    
  }
  else {// MIBL version:
    int nPreds = _tauSam.cols();
    
    // For MIBL, optimization is done via the closed-form update rule given by
    // Park and Casella (2008).
    double lambdaDenom = (_tauSam.colwise().mean()).sum();
    double newLambda   = sqrt((2.0 * double(nPreds)) / lambdaDenom); 
    _lambdas[0]                   = newLambda;
    _lambdaHistory(_emIterNum, 0) = _lambdas[0];
  }
  
  // Do some housekeeping:
  _emIterNum++;
  _storeGibbsSamples = false;
  _drawNum           = 0;
}// END updateLambdas()



///////////////////////// EXCEPTION HANDLING FUNCTIONS //////////////////////////



void MibrrGibbs::tauError(int errorCode) const
{
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



void MibrrGibbs::betaError(exception &e) const
{
  Rcpp::Rcerr << e.what() << endl;
  Rcpp::stop("Something terrible has occured while updating Beta.\nAbove this \
message, I've printed the that exception I caught.\nBeta luck next time.");
}



void MibrrGibbs::lambdaError() const
{
  if(_optMethod > 3) {
    if(_verbose) {
      Rcpp::Rcerr << "Lambda " << _optPrefix << "optimization failed with "; 
      Rcpp::Rcerr << _algName << ".\nNo further optimization algorithms ";
      Rcpp::Rcerr << "are available." << endl;
    }
    Rcpp::stop("Lambda cannot be optimized.");
  }
  else {
    if(_verbose) {
      Rcpp::Rcout << "Lambda " << _optPrefix << "optimization failed with ";
      Rcpp::Rcout << _algName << "\nRetrying with a different algorithm" << endl;
    }
  }
}



void MibrrGibbs::lambdaError(exception &e) const
{
  if(_optMethod > 3) {
    if(_verbose) {
      Rcpp::Rcerr << e.what();
      Rcpp::Rcerr << "Lambda " << _optPrefix << "optimization failed with "; 
      Rcpp::Rcerr << _algName << ", and returned the preceding exception.\nNo";
      Rcpp::Rcerr << "no further optimization algorithms are available." << endl;
    }
    Rcpp::stop("Lambda cannot be optimized.");
  }
  else {
    if(_verbose) {
      cerr << e.what();
      Rcpp::Rcout << "Lambda " << _optPrefix << "optimization failed with ";
      Rcpp::Rcout << _algName << ", and returned the preceding exception.\n";
      Rcpp::Rcout << "Retrying with a different algorithm" << endl;
    }
  }
}
