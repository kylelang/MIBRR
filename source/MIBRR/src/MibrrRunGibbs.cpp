//
// Created by Amir Masoud Abdol on 20/03/2020.
//

#include "MibrrData.h"
#include "MibrrGibbs.h"

json runGibbs(Eigen::MatrixXd           data,
              int                       nTargets,
              std::vector<std::vector<int>> missList,
              Eigen::VectorXi           respCounts,
              Eigen::VectorXd           lambda1,
              Eigen::VectorXd           lambda2,
              Eigen::VectorXd           l1Parms,
              Eigen::VectorXd           l2Parms,
              Eigen::VectorXd           sigmaStarts,
              Eigen::MatrixXd           tauStarts,
              Eigen::MatrixXd           betaStarts,
              int                       burnSams,
              int                       totalSams,
              int                       penType,
              double                    ridge,
              bool                      verbose,
              bool                      fullBayes,
              bool                      noMiss,
              bool                      savePpSams,
              bool                      useBetaMeans,
              bool                      finalRep,
              std::vector<unsigned int> seeds,
              int                       chain,
              bool                      intercept)
{
  // Disable multithreading for Eigen ops:
  Eigen::setNbThreads(1);

  // Unpack the list of missing row indices:
  std::vector< std::vector<int> > missIndices;
  for(int v = 0; v < nTargets; v++) missIndices.push_back(missList[v]);

  // Initialize the various classes needed below:
  MibrrData  mibrrData(data, missIndices, respCounts, noMiss);
  MibrrGibbs *mibrrGibbs = new MibrrGibbs[nTargets];

  // Initialize all parameters and setup the Gibbs sampler:
  for(int j = 0; j < nTargets; j++) {
    Eigen::VectorXd betaStartVec  = betaStarts.col(j);
    Eigen::ArrayXd  tauStartArray = tauStarts.col(j).array();

    // Define the type of regularization:
    mibrrGibbs[j].setPenType(penType);

    if(penType == 0)// Doing basic ridge?
      mibrrGibbs[j].setRidge(ridge);

    if(fullBayes) {// Fully Bayesian estimation?
      mibrrGibbs[j].doFullBayes();
      mibrrGibbs[j].setLam1Parms(l1Parms);

      if(penType == 2)// Doing MIBEN?
        mibrrGibbs[j].setLam2Parms(l2Parms);
    }

    if(savePpSams) mibrrGibbs[j].savePpSams();
    if(!intercept) mibrrGibbs[j].noIntercept();

    mibrrGibbs[j].seedRng(seeds[j]);

    mibrrGibbs[j].startParameters(betaStartVec,
                                  tauStartArray,
                                  sigmaStarts[j],
                                  lambda1[j],
                                  lambda2[j],
                                  useBetaMeans);

    mibrrGibbs[j].setTargetIndex(j);
    mibrrGibbs[j].setDoImp(!noMiss);
    mibrrGibbs[j].setNDraws(totalSams - burnSams);
    mibrrGibbs[j].setFinalRep(finalRep);

    if(!verbose) mibrrGibbs[j].beQuiet();
  }// CLOSE for(in j ==0; j < nTargets; j++)

  for(int i = 0; i < totalSams; i++) {// LOOP over Gibbs iterations
    // Print a nice progress message:
    if(verbose) {
      int marg, max;
      bool check0;
      if(i < burnSams) {
        marg   = burnSams % 5;
        max    = burnSams - marg;
        check0 = (i % (max / 5) == 0) & ((burnSams - i) > marg);
        if(check0) {
          std::cout << "Chain " << chain <<
                      ": Doing Gibbs burn-in iteration " << (i + 1) << " of " <<
                      burnSams << endl;
        }
      }
      else {
        marg   = (totalSams - burnSams) % 5;
        max    = (totalSams - burnSams) - marg;
        check0 = ((i - burnSams) % (max / 5) == 0) & ((totalSams - i) > marg);
        if(check0) {
          std::cout << "Chain " << chain <<
                      ": Doing Gibbs sampling iteration " << (i + 1) - burnSams <<
                      " of " << totalSams - burnSams << endl;
        }
      }
    }

    // Improve the output's aesthetics:
    bool check1 = verbose & ((i == burnSams - 1) || (i == totalSams - 1));
    if(check1) std::cout << "\n";

    for(int j = 0; j < nTargets; j++) {// LOOP over target variables
      // Compute new centers and scales for the jth target's predictors:
      mibrrData.updateMoments(j);

      // Update the Gibbs samples:
      mibrrGibbs[j].doGibbsIteration(mibrrData);

      // Start saving iterations after burn-in:
      if((i + 1) == burnSams)
        mibrrGibbs[j].startGibbsSampling(mibrrData);
    }
  }// CLOSE for (int i = 0; i < totalSams; i++)

//  RList outList(nTargets);
//  for(int j = 0; j < nTargets; j++)
//    outList[j] =
//        RList::create(Rcpp::Named("imps"  ) = mibrrGibbs[j].getImpSam(),
//                      Rcpp::Named("ppSams") = mibrrGibbs[j].getPpSam(),
//                      Rcpp::Named("beta"  ) = mibrrGibbs[j].getBetaSam(),
//                      Rcpp::Named("tau"   ) = mibrrGibbs[j].getTauSam(),
//                      Rcpp::Named("sigma" ) = mibrrGibbs[j].getSigmaSam(),
//                      Rcpp::Named("lambda") = mibrrGibbs[j].getLambdaSam()
//        );

//  return outList;

  return json{};
}// END runGibbs()
