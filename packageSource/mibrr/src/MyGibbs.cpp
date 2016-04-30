// Title:    Function definitions for the MyGibbs class
// Author:   Kyle M. Lang
// Created:  2014-AUG-24
// Modified: 2016-APR-29
// Purpose:  This class contains the Gibbs sampling-related functions
//           for the MIBRR package.

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

#include "MyGibbs.hpp"

MyGibbs::MyGibbs() 
{
  // _betas, _taus, _nEmIters, and _nDraws  
  // all need user-supplied starting values.
  _sigma             = 0.0;
  _lambdas           = VectorXd(2);
  _storeGibbsSamples = false;
  _drawNum           = 0;
  _emIterNum         = 0;
  _lambdaWindow      = 1;
  _emConvTol         = 1.0e-5;
  _optIterCount      = 0;
}


MyGibbs::~MyGibbs() 
{
}


///// ACCESSORS /////

VectorXd MyGibbs::getBetas() const 
{ 
  return _betas; 
}


ArrayXd MyGibbs::getTaus() const 
{ 
  return _taus; 
}


double MyGibbs::getSigma() const 
{ 
  return _sigma; 
}


MatrixXd MyGibbs::getBetaSam() const 
{ 
  return _betaSam; 
}


ArrayXXd MyGibbs::getTauSam() const 
{ 
  return _tauSam; 
}


VectorXd MyGibbs::getSigmaSam() const 
{ 
  return _sigmaSam; 
}


MatrixXd MyGibbs::getImpSam() const
{
  return _impSam;
}


MatrixXd MyGibbs::getLambdaHistory() const
{
  return _lambdaHistory;
}


VectorXd MyGibbs::getLambdas() const 
{ 
  VectorXd outLam;
  if(_useElasticNet) outLam = _lambdas;
  else               outLam = VectorXd::Constant(1, _lambdas[0]);
  return outLam; 
}


double MyGibbs::getLambdas(int lambdaNumber) const 
{ 
  return _lambdas[lambdaNumber - 1]; 
}


int MyGibbs::getNDraws() const
{
  return _nDraws;
}


///// MUTATORS /////

void MyGibbs::setBetas(VectorXd &betas) 
{
  _betas = betas;
}


void MyGibbs::setTaus(ArrayXd &taus) 
{
  _taus = taus;
}


void MyGibbs::setSigma(double sigma) 
{
  _sigma = sigma;
}


void MyGibbs::setLambdas(VectorXd& lambdas) 
{  
  _lambdas = lambdas;
}


void MyGibbs::setLambdas(double lambda1, 
			 double lambda2) 
{
  _lambdas[0] = lambda1;
  _lambdas[1] = lambda2;
}


void MyGibbs::setLambdas(double lambda) 
{
  _lambdas[0] = lambda;
  _lambdas[1] = 0.0;
}


void MyGibbs::setLambdas() 
{  
  int startRow = _emIterNum - _lambdaWindow;
  int nCols = _useElasticNet ? 2 : 1;
  VectorXd pooledLambdas = 
    _lambdaHistory.block(startRow, 0, _lambdaWindow, nCols).colwise().mean();

  if(_useElasticNet) _lambdas = pooledLambdas.transpose();
  else _lambdas[0] = pooledLambdas[0];
}


void MyGibbs::setNDraws(int nDraws) 
{
  _nDraws = nDraws;
}


void MyGibbs::setNEmIters(int nEmIters) 
{
  _nEmIters = nEmIters;
}


void MyGibbs::setTargetIndex(int targetIndex) 
{
  _targetIndex = targetIndex;
}


void MyGibbs::setVerbosity(bool verboseIters,
			   bool verboseErrors) 
{
  _verboseIters  = verboseIters;
  _verboseErrors = verboseErrors;
}


void MyGibbs::doMiben()
{
  _useElasticNet = true;
}


void MyGibbs::doMibl()
{
  _useElasticNet = false;
}


void MyGibbs::startGibbsSampling(MyData &myData)
{
  int nObs = myData.nObs();
  _storeGibbsSamples = true;
  
  _betaSam  = MatrixXd(_nDraws, _betas.size());
  _tauSam   = MatrixXd(_nDraws, _taus.size());
  _sigmaSam = VectorXd(_nDraws);
  _impSam   = MatrixXd(_nDraws, nObs);
}


void MyGibbs::stopGibbsSampling()
{
  _storeGibbsSamples = false;
}


void MyGibbs::setupOptimizer(int nEmIters,
			     int lambdaWindow,
			     double emConvTol)
{
  _nEmIters = nEmIters;
  _lambdaWindow = lambdaWindow;
  if(_useElasticNet) _emConvTol = emConvTol;
}


//void MyGibbs::setupOptimizer(int nEmIters,
//			     int lambdaWindow)
//{
//  _nEmIters     = nEmIters;
//  _lambdaWindow = lambdaWindow;
//}


void MyGibbs::startParameters(VectorXd &betaStarts,
			      ArrayXd &tauStarts,
			      double sigmaStart,
			      double lambda1Start,
			      double lambda2Start)
{
  _betas         = betaStarts;
  _taus          = tauStarts;
  _sigma         = sigmaStart;
  _lambdas[0]    = lambda1Start;
  _lambdas[1]    = lambda2Start;
  int lamNum     = _useElasticNet ? 2 : 1;
  _lambdaHistory = MatrixXd(_nEmIters, lamNum);
}


void MyGibbs::startParameters(VectorXd &betaStarts,
			      ArrayXd &tauStarts,
			      double sigmaStart,
			      VectorXd &lambdaStartVec)
{
  _betas         = betaStarts;
  _taus          = tauStarts;
  _sigma         = sigmaStart;
  _lambdas       = lambdaStartVec;
  int lamNum     = _useElasticNet ? 2 : 1;
  _lambdaHistory = MatrixXd(_nEmIters, lamNum);
}


void MyGibbs::restartParameters(MyData &myData)
{
  _betas = _betaSam.bottomRows(_betaSam.rows() - 1).colwise().mean().transpose();
  _taus = _tauSam.bottomRows(_tauSam.rows() - 1).colwise().mean().transpose();
  _sigma = _sigmaSam.tail(_sigmaSam.size() - 1).mean();
  
  VectorXd meanImps = _impSam.bottomRows(_impSam.rows() - 1).colwise().mean();
  myData.fillMissing(meanImps, _targetIndex);
}


///// PARAMETER UPDATE FUNCTIONS /////

void MyGibbs::updateTaus(MyData &myData)
{
  double lambda1 = _lambdas[0];
  double lambda2 = _lambdas[1];
  int nPreds = myData.nPreds(); 
  ArrayXd tauMeans;
  double tauScale = -1.0;// -1 to ensure an exception if try() fails
  
  try {
    if(_useElasticNet) {// Miben Version
      ArrayXd tauMeansNumerator = ArrayXd::Constant(nPreds, lambda1);
      tauMeans = tauMeansNumerator.sqrt() /
	(2.0 * lambda2 * _betas.tail(nPreds).array().abs());
      tauScale = lambda1 / (4.0 * lambda2 * _sigma); 
    }
    else {// MIBL Version
      double tauMeansNumerator = pow(lambda1, 2) * _sigma;
      tauMeans =
	(tauMeansNumerator / _betas.tail(nPreds).array().square()).sqrt();
      tauScale = pow(lambda1, 2);
    }
    
    if((tauMeans <= 0.0).any()) throw 1;
    if(tauScale <= 0.0) throw 2;
  }
  catch (int e) { tauError(e); }
  
  // Draw new values of the auxiliary penalty parameters:
  ArrayXd tmpDraws(nPreds);
  for(int i = 0; i < nPreds; i++)
    tmpDraws[i] = myData.drawInvGauss(tauMeans[i], tauScale);
  
  ArrayXd newTaus;
  if(_useElasticNet) newTaus = (tmpDraws + 1.0) / tmpDraws; // MIBEN Version   
  else               newTaus = 1.0 / tmpDraws;              // MIBL Version
  
  _taus = newTaus;// Store the updated Taus
  
  // Add the updated Taus to their Gibbs sample:
  if(_storeGibbsSamples) _tauSam.row(_drawNum) = newTaus.transpose();
}// END updateTaus ()


void MyGibbs::updateBetas(MyData &myData)
{
  int nPreds = myData.nPreds();
  int nObs = myData.nResponses(_targetIndex);
  LDLT <MatrixXd> aMatrixCholesky;
  MatrixXd tmpMat;

  if(_useElasticNet) {// MIBEN Version
    double lambda2 = _lambdas[1];
    VectorXd transformedTaus = _taus / (_taus - 1.0);
    tmpMat = lambda2 * transformedTaus.asDiagonal();    
  }
  else {// MIBL Version
    VectorXd transformedTaus = 1.0 / _taus;
    tmpMat = transformedTaus.asDiagonal();
  }
  
  MatrixXd aMatrix = myData.getIVs(_targetIndex).transpose() *
    myData.getIVs(_targetIndex) + tmpMat;
  
  VectorXd betaMeans;
  MatrixXd betaCovariances;
  try {
    // Compute the Cholesky decomposition of the IV's penalized crossproducts,
    // and use it to compute the moments of beta's fully conditional posterior:
    aMatrixCholesky.compute(aMatrix);
    betaMeans = aMatrixCholesky.solve(myData.getIVs(_targetIndex).transpose() *
				      myData.getDV(_targetIndex)); 
    
    tmpMat = _sigma * MatrixXd::Identity(nPreds, nPreds);
    betaCovariances = aMatrixCholesky.solve(tmpMat);
  }
  catch(exception &e) { betaError(e); }
  
  // Draw a new value of the intercept term:
  VectorXd newBetas(nPreds + 1);
  newBetas[0] = R::rnorm(myData.getDV(_targetIndex).mean(),
			 sqrt(_sigma / double(nObs)));
  
  // Draw new values of the regression slope coefficients:
  newBetas.tail(nPreds) = myData.drawMVN(betaMeans, betaCovariances);
  
  _betas = newBetas;// Store the updated Betas
  
  // Add the updated Betas to their Gibbs sample:
  if(_storeGibbsSamples) _betaSam.row(_drawNum) = newBetas.transpose();
}// END updateBetas ()


void MyGibbs::updateSigma(MyData &myData)
{
  int nPreds = myData.nPreds();
  int nObs = myData.nResponses(_targetIndex);
  double lambda1 = _lambdas[0];
  double lambda2 = _lambdas[1];
  VectorXd tmpBiasVector = VectorXd::Ones (nObs);
  double newSigma;

  // Compute the regularized residual sum of squares:
  double sse =
    (myData.getDV(_targetIndex) - _betas[0] * tmpBiasVector -
     myData.getIVs(_targetIndex) * _betas.tail(nPreds)).transpose() *
    (myData.getDV(_targetIndex) - _betas[0] * tmpBiasVector -
     myData.getIVs(_targetIndex) * _betas.tail(nPreds));
  
  if(_useElasticNet) {// MIBEN Version
    double scaleSum =
      (_taus / (_taus - 1.0) * _betas.tail(nPreds).array().square()).sum();
    
    double sigmaShape = (double(nObs) / 2.0) + double(nPreds);
    double sigmaScale = 0.5 * (sse + lambda2 * scaleSum +
			       (pow(lambda1, 2) / (4.0 * lambda2)) * _taus.sum()
			       );
    
    bool isDrawValid = false;
    double testDraw;
    while(!isDrawValid) {// Rejection sampling to draw a Sigma variate
      testDraw = myData.drawInvGamma(sigmaShape, sigmaScale);
      double thresholdDraw = unif_rand();
      double uiGammaShape = pow(lambda1, 2) / (8.0 * testDraw * lambda2);
      double uiGammaDraw = myData.calcIncGamma(0.5, uiGammaShape, false);
      isDrawValid = log(thresholdDraw) <= (double(nPreds) * log(tgamma(0.5))) -
	(double(nPreds) * log(uiGammaDraw));
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
    newSigma = myData.drawInvGamma(sigmaShape, sigmaScale);
  }
  
  _sigma = newSigma;// Store the updated Sigma
  
  // Add the updated sigma to its Gibbs sample:
  if(_storeGibbsSamples) _sigmaSam[_drawNum] = newSigma;
}// END updateSigma ()


void MyGibbs::updateImputations(MyData &myData)
{
  int nObs = myData.nObs();
  int nPreds = myData.nPreds();
  ArrayXb nonresponseVector = myData.getNonresponseVector(_targetIndex);
  VectorXd tmpBiasVector = VectorXd::Ones(nObs);
  double tmpSd = sqrt(_sigma);
  VectorXd errorVector(nObs);
  
  // Draw the residual error terms for the imputation model:
  for(int i = 0; i < nObs; i++) errorVector[i] = R::rnorm(0.0, tmpSd);
  
  // Draw a vector of imputations from the posterior predictive distribution
  // of the missing data:
  VectorXd tmpImps = _betas[0] * tmpBiasVector +
    myData.getFullIVs(_targetIndex) * _betas.tail(nPreds) + errorVector;
  
  // Replace the missing data in the target variable with the imputations:
  for(int i = 0; i < nObs; i++)
    if(nonresponseVector[i] == true)
      myData.setElement(tmpImps[i], i, _targetIndex);
  
  // Add the updated imputations to their Gibbs sample:
  if(_storeGibbsSamples) _impSam.row(_drawNum) = tmpImps.transpose();
}// END updateImputations ()


void MyGibbs::doGibbsIteration(MyData &myData)
{  
  updateTaus(myData);
  updateBetas(myData);
  updateSigma(myData);
  updateImputations(myData);
  
  if(_storeGibbsSamples) _drawNum++;
}// END doGibbsIteration ()


///// MCEM OPTIMIZATION FUNCTIONS /////

double MyGibbs::eNetLambdaObjective(const std::vector<double> &lambdaVec,
				    std::vector<double> &gradVec,
				    void *extraOptData)
{
  bool nanLambdas = isnan(lambdaVec[0]) | isnan(lambdaVec[1]) |
    lambdaVec[0] != lambdaVec[0] | lambdaVec[1] != lambdaVec[1];
  
  if(nanLambdas) 
    throw invalid_argument("An Elastic Net penalty parameter is NaN.");  
  
  _optIterCount++; // Track the number of function evaluations
  
  int nPreds = _tauSam.cols();
  int nSams = _tauSam.rows();

  // BEGIN compute LL and gradient terms.
  ArrayXd tmpArray = 2.0 * _sigmaSam;
  double w1 = nPreds * log(lambdaVec[0]);
  ArrayXd w2 = lambdaVec[1] / tmpArray;
  ArrayXd w3 = 1 / tmpArray;
  
  tmpArray = 8.0 * _sigmaSam * lambdaVec[1];
  ArrayXd w4 = pow(lambdaVec[0], 2) / tmpArray;
  double w5 = nPreds / lambdaVec[0];
  double w6 = (nPreds * lambdaVec[0]) / (4.0 * lambdaVec[1]);
  
  tmpArray = 4.0 * _sigmaSam * lambdaVec[1];
  ArrayXd w7 = lambdaVec[0] / tmpArray;
  double w8 = (nPreds * pow(lambdaVec[0], 2)) / (8.0 * pow(lambdaVec[1], 2));
  
  tmpArray = 8.0 * _sigmaSam * pow(lambdaVec[1], 2);
  ArrayXd w9 = pow(lambdaVec[0], 2) / tmpArray;
  
  ArrayXd uiGammaArray(nSams); 
  for(int i = 0; i < nSams; i++)
    uiGammaArray[i] = R::pgamma(w4[i], 0.5, 1.0, 0, 0) * tgamma(0.5);
  
  ArrayXd term1 = uiGammaArray.log();
  ArrayXd term2 = ((_tauSam / (_tauSam - 1)) *
		   _betaSam.rightCols(nPreds).array().square()).rowwise().sum();
  // END compute LL and gradient terms.

  if(!gradVec.empty()) {// Calculate the gradient vector:
    ArrayXd term3 = (1 / uiGammaArray) * (1 / w4.sqrt()) * (-1.0 * w4).exp() *
      (1 / _sigmaSam.array());
    gradVec[0] = (w5 + (w6 * term3) - w7 * _tauSam.rowwise().sum()).mean();
    gradVec[1] =
      (-1.0 * (w8 * term3) - (w3 * term2) + w9 * _tauSam.rowwise().sum()).mean();
  }
  
  // Calculate and return the objective value:
  return (w1 - (nPreds * term1) - (w2 * term2) -
	  (w4 * _tauSam.rowwise().sum())).mean();
}// END eNetLambdaObjective()


// As suggested by the nlopt authors, specify a simple  wrapper function for
// eNetLambdaObjective() so that nlopt will run inside the MyGibbs class:
double eNetObjectiveWrap(const std::vector<double> &lambdaVec, 
			 std::vector<double> &gradVec, 
			 void *extraOptData) 
{
  MyGibbs *obj = static_cast<MyGibbs *>(extraOptData);   
  return obj -> eNetLambdaObjective(lambdaVec, gradVec, extraOptData);
}// END eNetObjectiveWrap()


void MyGibbs::optimizeMibenLambdas(bool preOptimize) 
{
  std::vector<double> lambdaVec(2), gradVec(2), lamBounds(2);
  lambdaVec[0] = _lambdas[0];
  lambdaVec[1] = _lambdas[1];
  gradVec[0] = 0.0;
  gradVec[1] = 0.0;
  lamBounds[0] = 1.0e-4; // Set low bounds slightly above zero to avoid dividing 
  lamBounds[1] = 1.0e-4; // by zero when no regularization is needed.
  nlopt::algorithm optAlg;
  int optMethod = 0;
  string algName, outPrefix; 
  bool optimized = false;
  
  while(!optimized) {//Try different algorithms until convergence
    if(preOptimize) {
      outPrefix = "pre-";
      if(optMethod == 0) {
	optAlg = nlopt::LN_BOBYQA;
	algName = "BOBYQA";
      }
      else if(optMethod == 1) {
	optAlg = nlopt::LN_COBYLA;
	algName = "COBYLA";
      }
      else if(optMethod == 2) {
	optAlg = nlopt::LN_SBPLX;
	algName = "SBPLX";
      }
      else if(optMethod == 3) {
	optAlg = nlopt::LN_PRAXIS;
	algName = "PRAXIS";
      }
    }
    else {
      outPrefix = "";
      if(optMethod == 0) {
	optAlg = nlopt::LD_MMA;
	algName = "MMA";
      }
      else if(optMethod == 1) {
	optAlg = nlopt::LD_VAR2;
	algName = "VAR2";
      }
      else if(optMethod == 2) {
	optAlg = nlopt::LD_VAR1;
	algName = "VAR1";
      }
      else if(optMethod == 3) {
	optAlg = nlopt::LD_LBFGS;
	algName = "LBFGS";
      }
    }// END if(preOptimize}
    
    nlopt::opt myOptimizer(optAlg, 2);       // Initialize the optimizer object
    myOptimizer.set_max_objective(eNetObjectiveWrap, this);
    myOptimizer.set_ftol_rel(_emConvTol);    // Stopping criterion
    myOptimizer.set_lower_bounds(lamBounds); // Force positive Lambdas	
    
    nlopt::result myResult;
    double maxLL = 0.0;
    try {
      myResult = myOptimizer.optimize(lambdaVec, maxLL);  
      if(myResult < 0) {            // Catch nlopt failure codes
	optMethod ++;               // Try the next algorithm
	lambdaVec[0] = _lambdas[0]; // Reset the lambdas
	lambdaVec[1] = _lambdas[1];
	lambdaError(optMethod, outPrefix, algName);
	//if(optMethod > 3) {
	//  if(_verboseErrors) {
	//  cerr << "Lambda " << outPrefix << "optimization failed with "; 
	//  cerr << algName << "\nNo more optimization algorithms available";
	//  cerr << "\nThis program will now crash" << endl;
	//}
	//throw -42;// If there are no more algorithms to try, 
	  // an uncaught exception will crash the program
	//}
	//else {
	//  if(_verboseErrors) {
	//    cerr << "Lambda " << outPrefix << "optimization failed with ";
	//    cerr << algName << "\nRetrying with a different algorithm" << endl;
	//  }
	//}
      }
      else if(myResult > 0) {// Successful convergence!
	if(_verboseIters) {
	  cout << "Lambdas " << outPrefix << "optimized with " << algName;
	  cout << " in " << _optIterCount << " iterations" << endl;
	}	
	optimized = true;
      } 
    }
    catch(int &e) {
      throw e;
    }
    catch(exception &e) {//Catch Lambdas == NaN
      optMethod ++;// Try the next algorithm
      lambdaVec[0] = _lambdas[0];
      lambdaVec[1] = _lambdas[1];// Reset lambdas
      lambdaError(e, optMethod, outPrefix, algName);
	//if(optMethod > 3) {
	//if(_verboseErrors) {
	//cerr << e.what(); 
	//cerr << "\nLambda " << outPrefix << "optimization failed with "; 
	//cerr << algName << "\nNo more optimization algorithms available";
	//cerr << "\nThis program will now crash" << endl;
	//}
	//throw -42;// If there are no more algorithms to try, 
	// an uncaught exception will crash the program
	//}
	//else {
	//if(_verboseErrors) {
	//cerr << e.what();
	//cerr << "\nLambda " << outPrefix << "optimization failed with ";
	//cerr << algName << "\nRetrying with a different algorithm" << endl;
	//}
	//}
	}  
    _optIterCount = 0;
  }// END while(!optimized)      
  // Store the updated penalty parameters:
  _lambdas[0] = lambdaVec[0];
  _lambdas[1] = lambdaVec[1];
  _lambdaHistory.row(_emIterNum) = _lambdas.transpose();
}// END optimizeMibenLambdas()



void MyGibbs::updateLambdas()
{
  if(_useElasticNet) {// MIBEN version
    // For MIBEN, optimization is done in two stages:
    // 1) A rough, gradient-free, pre-optimization is done to move the estimates
    //    into the neighborhood of the MLE.
    // 2) Gradient-based optimization is used to fine-tune the estimates from
    //    Step (1).
    optimizeMibenLambdas(true);  // Pre-optimization
    optimizeMibenLambdas(false); // Optimization    
  }
  else {// MIBL version:
    int nPreds = _tauSam.cols();
    
    // For MIBL, optimization is done via the closed-form update rule given by
    // Park and Casella (2008).
    double lambdaDenominator = (_tauSam.colwise().mean()).sum();
    double newLambda = sqrt((2.0 * double(nPreds)) / lambdaDenominator); 
    _lambdas[0] = newLambda;
    _lambdaHistory(_emIterNum, 0) = _lambdas[0];
  }
  
  // Do some housekeeping:
  _emIterNum++;
  _storeGibbsSamples = false;
  _drawNum = 0;
}// END updateLambdas()
