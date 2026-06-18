#pragma once
#include <vector>
// Import Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
using mat_type = Eigen::MatrixXcd;

Matrix::t get_sparse_from_eigen(Eigen::MatrixXcd mat)
{
  std::vector<double> values;
  std::vector<int> rows;
  std::vector<int> cols;

  for (int i = 0; i < mat.rows(); i++)
  {

    for (int j = 0; j < mat.cols(); j++)
    {

      if (std::abs(mat.coeff(i, j)) > 1e-9)
      {
        // std::cout<< "her"<<std::endl;
        if (std::abs(mat.coeff(i, j).real()) > 1e-9)
        {
          values.push_back(mat.coeff(i, j).real());
          rows.push_back(i); // row index
          cols.push_back(j); // col index (here it is equal to k)
          values.push_back(mat.coeff(i, j).real());
          rows.push_back(i + mat.cols()); // row index
          cols.push_back(j + mat.cols()); // col index (here it is equal to k)
                                          // it.index(); // inner index, here it is equal to it.row()
        }

        if (std::abs(mat.coeff(i, j).imag()) > 1e-9)
        {
          values.push_back(mat.coeff(i, j).imag());
          rows.push_back(i + mat.cols()); // row index
          cols.push_back(j);              // col index (here it is equal to k) // correct with minus?
          values.push_back(-1. * mat.coeff(i, j).imag());
          rows.push_back(j + mat.cols()); // row index
          cols.push_back(i + mat.cols());
        }
      }
    }
  }

  Matrix::t Alpha = Matrix::sparse(2 * mat.rows(), 2 * mat.cols(), nint(rows), nint(cols), ndou(values));
  return Alpha;
}
