#pragma once
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>

#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "definitions.hpp"
// defenitions
using namespace mosek::fusion;
using namespace monty;
const double pi = std::acos(-1.0);
//
std::shared_ptr<ndarray<int, 1>> nint(const std::vector<int> &X) { return new_array_ptr<int>(X); }
std::shared_ptr<ndarray<double, 1>> ndou(const std::vector<double> &X) { return new_array_ptr<double>(X); }

template <class T>
int getIndex(std::vector<T> v, T K)
{
  auto it = find(v.begin(), v.end(), K);

  // If element was found
  if (it != v.end())
  {

    // calculating the index
    // of K
    int index = it - v.begin();
    return index;
  }
  else
  {
    // If the element is not
    // present in the vector
    return -1;
  }
}
struct matrix_organizer
{
  std::vector<int_pair> matrix_positions;
  std::vector<double> matrix_values;
  std::vector<double> b;
  int variable_index{0};
  bool has_elements_ = false;

  void add_values(int_pair position, double value)
  {
    has_elements_ = true;

    auto index = getIndex(matrix_positions, position);

    // Check if the target value was found
    if (index < 0)
    {
      matrix_positions.push_back(position);
      matrix_values.push_back(value);
    }
    else
    {
      matrix_values[index] += value;
    }
  }
  Matrix::t make_matrix(int dim1, int dim2)
  {
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> T;
    int i = 0;
    for (auto &p : matrix_positions)
    {
      if (std::abs(matrix_values[i]) > 1e-9)
      {
        rows.push_back(p.first);
        cols.push_back(p.second);
        T.push_back(matrix_values[i]);
      }
      // std::cout<< "values "<<matrix_values[i]<<std::endl;

      i++;
    }

    return Matrix::sparse(dim1, dim2, nint(rows), nint(cols), ndou(T));
  }
};