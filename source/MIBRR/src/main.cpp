//
// Created by Amir Masoud Abdol on 13/03/2020.
//

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "MibrrDefs.h"
#include "MibrrGibbs.h"

#include "utils/ioutils.h"

using json = nlohmann::json;

int main() {

  std::ifstream configfile("../config.json");

  json configs;
  configfile >> configs;

  // Reading the CSV file directory to an Eigen Matrix
  auto csvdata = csvToEigenMatrix<Eigen::MatrixXd, Eigen::VectorXd, double>(
      configs["datafile"], true, true);

  // Reading the data from a JSON _2d_ array (object)
  Eigen::MatrixXd data =
      toEigenMatrix<Eigen::MatrixXd, Eigen::VectorXd, double>(configs["data"]);

  int nTargets = configs["nTargets"].get<int>();

  auto missList = configs["missLists"].get<std::vector<std::vector<int>>>();

  // Reading a 1D array into a vector
  // I was hopping that my template makes it easier to read different types of
  // data but for some reason Eigen doesn't like when I construct an integer
  // vector the way I do. I'm going to fix this at some point.
  Eigen::VectorXd respCounts_d =
      toEigenVector<Eigen::VectorXd, double>(configs["respCounts"]);
  Eigen::VectorXi respCounts = respCounts_d.cast<int>();

  Eigen::VectorXd lambda1 =
      toEigenVector<Eigen::VectorXd, double>(configs["lambda1"]);

  Eigen::VectorXd lambda2 =
      toEigenVector<Eigen::VectorXd, double>(configs["lambda2"]);

  Eigen::VectorXd l1Parms =
      toEigenVector<Eigen::VectorXd, double>(configs["l1Parms"]);

  Eigen::VectorXd l2Parms =
      toEigenVector<Eigen::VectorXd, double>(configs["l2Parms"]);

  Eigen::VectorXd sigmaStarts =
      toEigenVector<Eigen::VectorXd, double>(configs["sigmaStarts"]);

  Eigen::MatrixXd tauStarts =
      toEigenMatrix<Eigen::MatrixXd, Eigen::VectorXd, double>(
          configs["tauStarts"]);

  Eigen::MatrixXd betaStarts =
      toEigenMatrix<Eigen::MatrixXd, Eigen::VectorXd, double>(
          configs["betaStarts"]);

  int burnSams = configs["burnSams"].get<int>();
  int totalSams = configs["totalSams"].get<int>();
  int penType = configs["penType"].get<int>();
  double ridge = configs["ridge"].get<double>();
  bool verbose = configs["verbose"].get<bool>();
  bool fullBayes = configs["fullBayes"].get<bool>();
  bool noMiss = configs["noMiss"].get<bool>();
  bool savePpSams = configs["savePpSams"].get<bool>();
  bool useBetaMeans = configs["useBetaMeans"].get<bool>();
  bool finalRep = configs["finalRep"].get<bool>();
  std::vector<unsigned int> seeds =
      configs["seeds"].get<std::vector<unsigned int>>();
  int chain = configs["chain"].get<int>();
  bool intercept = configs["intercept"].get<bool>();

  //  auto j_out =
  //      runGibbs(data, nTargets, missList, respCounts, lambda1, lambda2,
  //      l1Parms,
  //               l2Parms, sigmaStarts, tauStarts, betaStarts, burnSams,
  //               totalSams, penType, ridge, verbose, fullBayes, noMiss,
  //               savePpSams, useBetaMeans, finalRep, seeds, chain, intercept);

  return 0;
}
