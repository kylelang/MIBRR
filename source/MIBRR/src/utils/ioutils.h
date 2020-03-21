//
// Created by Amir Masoud Abdol on 21/03/2020.
//

#ifndef MIBRR_SERIALIZERS_H
#define MIBRR_SERIALIZERS_H

#include "parser.hpp"

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

template <typename VectorType, typename T>
VectorType to_eigen_vector(std::vector<T> row) {
  return Eigen::Map<VectorType, Eigen::Unaligned>(row.data(), row.size());
}

template <typename MatrixType, typename VectorType, typename T>
MatrixType to_eigen_matrix(const std::vector<std::vector<T>> &data) {
  int n_cols, n_rows;
  n_rows = data.size();
  n_cols = data[0].size();
  Eigen::MatrixXd matrix(n_rows, n_cols);
  for (int i = 0; i < n_rows; ++i) {
    int n = data[i].size();
    assert(n_cols == n && "data is not representing a matrix");
    matrix.row(i) = to_eigen_vector<VectorType, T>(data[i]);
  }
  return matrix;
}

template <typename MatrixType, typename VectorType, typename T = double>
MatrixType csv_to_eigen(const std::string filename, bool has_index = true,
                        bool has_header = true) {
  std::vector<std::vector<T>> vdata;

  std::ifstream datafile(filename);
  aria::csv::CsvParser csvparser(datafile);

  bool is_index = has_index;

  for (auto &row : csvparser) {
    if (has_header) {
      has_header = false;
      continue;
    }

    vdata.push_back(std::vector<double>());

    bool is_index = has_index;
    for (auto field : row) {
      if (is_index) {
        is_index = false;
        continue;
      }

      if (field == "NA")
        field = "NAN";

      vdata.back().push_back(stod(field));
    }
  }

  return to_eigen_matrix<MatrixType, VectorType, T>(vdata);
}

#endif // MIBRR_SERIALIZERS_H
