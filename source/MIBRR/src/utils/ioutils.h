//
// Created by Amir Masoud Abdol on 21/03/2020.
//

#ifndef MIBRR_SERIALIZERS_H
#define MIBRR_SERIALIZERS_H

#include <cmath>

#include "rapidcsv.h"

namespace rapidcsv {
template <>
void Converter<double>::ToVal(const std::string &pStr, double &pVal) const {
  if (pStr == "NA") {
    pVal = NAN;
  } else {
    pVal = stod(pStr);
  }
}
} // namespace rapidcsv

template <typename VectorType, typename T>
VectorType toEigenVector(json &j_vector) {
  return Eigen::Map<VectorType, Eigen::Unaligned>(
      j_vector.get<std::vector<T>>().data(), j_vector.size());
}

template <typename MatrixType, typename VectorType, typename T>
MatrixType toEigenMatrix(json &j_matrix) {
  int n_cols, n_rows;
  n_rows = j_matrix.size();
  n_cols = j_matrix[0].size();
  Eigen::MatrixXd matrix(n_rows, n_cols);

  for (int i = 0; i < n_rows; ++i) {
    matrix.row(i) = toEigenVector<VectorType, T>(j_matrix[i]);
  }

  return matrix;
}

template <typename VectorType, typename T>
VectorType toEigenVector(std::vector<T> row) {
  return Eigen::Map<VectorType, Eigen::Unaligned>(row.data(), row.size());
}

template <typename MatrixType, typename VectorType, typename T>
MatrixType toEigenMatrix(const std::vector<std::vector<T>> &data) {
  int n_cols, n_rows;
  n_rows = data.size();
  n_cols = data[0].size();
  Eigen::MatrixXd matrix(n_rows, n_cols);
  for (int i = 0; i < n_rows; ++i) {
    int n = data[i].size();
    assert(n_cols == n && "data is not representing a matrix");
    matrix.row(i) = toEigenVector<VectorType, T>(data[i]);
  }
  return matrix;
}

template <typename MatrixType, typename VectorType, typename T = double>
MatrixType csvToEigenMatrix(const std::string filename, bool has_index = true,
                            bool has_header = true) {

  std::vector<std::vector<T>> vdata;

  rapidcsv::Document doc(filename, rapidcsv::LabelParams());

  for (int i = 0; i < doc.GetRowCount(); ++i) {
    vdata.push_back(doc.GetRow<double>(i));
  }

  return toEigenMatrix<MatrixType, VectorType, T>(vdata);
}

auto findMisses(Eigen::MatrixXd &data) {
  std::vector<std::vector<int>> misses;

  for (int i = 0; i < data.cols(); ++i) {
    misses.push_back(std::vector<int>());
    for (int j = 0; j < data.rows(); ++j) {
      if (isnan(data.col(i)[j])) {
        misses.back().push_back(j);
      }
    }
  }

  return misses;
}

#endif // MIBRR_SERIALIZERS_H
