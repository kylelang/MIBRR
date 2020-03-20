//
// Created by Amir Masoud Abdol on 13/03/2020.
//

#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>

using json = nlohmann::json;

#include "MibrrDefs.h"
#include "MibrrGibbs.h"

template <typename VectorType, typename T>
VectorType to_eigen_vector(json &j_vector) {
  return Eigen::Map<VectorType, Eigen::Unaligned>(
      j_vector.get<std::vector<T>>().data(), j_vector.size());
}

template <typename MatrixType, typename VectorType, typename T>
MatrixType to_eigen_matrix(json &j_matrix) {
  int n_cols, n_rows;
  n_rows = j_matrix.size();
  n_cols = j_matrix[0].size();
  Eigen::MatrixXd matrix(n_rows, n_cols);
  for (int i = 0; i < n_rows; ++i) {
    matrix.row(i) = to_eigen_vector<VectorType, T>(j_matrix[i]);
  }

  return matrix;
}

int main() {

  std::ifstream configfile("../config.json");

  json configs;
  configfile >> configs;

  Eigen::MatrixXd data =
      to_eigen_matrix<Eigen::MatrixXd, Eigen::VectorXd, double>(
          configs["data"]);

  int nTargets = configs["nTargets"].get<int>();

  auto missList = configs["missList"].get<std::vector<std::vector<int>>>();

  Eigen::VectorXd respCounts_d =
      to_eigen_vector<Eigen::VectorXd, double>(configs["respCounts"]);
  // to_eigen_vector<Eigen::VectorXi, int>(configs["respCounts"]);
  Eigen::VectorXi respCounts = respCounts_d.cast<int>();

  Eigen::VectorXd lambda1 =
      to_eigen_vector<Eigen::VectorXd, double>(configs["lambda1"]);

  Eigen::VectorXd lambda2 =
      to_eigen_vector<Eigen::VectorXd, double>(configs["lambda2"]);

  Eigen::VectorXd l1Parms =
      to_eigen_vector<Eigen::VectorXd, double>(configs["l1Parms"]);

  Eigen::VectorXd l2Parms =
      to_eigen_vector<Eigen::VectorXd, double>(configs["l2Parms"]);

  Eigen::VectorXd sigmaStarts =
      to_eigen_vector<Eigen::VectorXd, double>(configs["sigmaStarts"]);

  Eigen::MatrixXd tauStarts =
      to_eigen_matrix<Eigen::MatrixXd, Eigen::VectorXd, double>(
          configs["tauStarts"]);

  Eigen::MatrixXd betaStarts =
      to_eigen_matrix<Eigen::MatrixXd, Eigen::VectorXd, double>(
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

  auto j_out =
      runGibbs(data, nTargets, missList, respCounts, lambda1, lambda2, l1Parms,
               l2Parms, sigmaStarts, tauStarts, betaStarts, burnSams, totalSams,
               penType, ridge, verbose, fullBayes, noMiss, savePpSams,
               useBetaMeans, finalRep, seeds, chain, intercept);

  return 0;
}
