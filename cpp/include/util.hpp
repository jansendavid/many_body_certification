#pragma once
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>
#include"spins.hpp"
#include<unordered_map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include"definitions.hpp"
// defenitions
const double pi = std::acos(-1.0);
//
struct matrix_organizer{
  std::vector<int_pair> matrix_positions;
  std::vector<double> matrix_values;
  std::vector<double> b;
  int variable_index{0};
  bool has_elements_=false;
  
  void add_values(int_pair position,double value)
  {
    has_elements_=true;

    auto index=getIndex(matrix_positions, position);
    
    
    // Check if the target value was found
    if(index<0)
      {
	matrix_positions.push_back(position);
	matrix_values.push_back(value);
	
      }
    else
      {
	matrix_values[index]+=value;
      }
    
  }
  Matrix::t make_matrix(int dim1, int dim2)
  {
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> T;
    int i=0;
    for(auto& p:matrix_positions )
      {
        if(std::abs(matrix_values[i])>1e-9)
        {
	rows.push_back(p.first);
	cols.push_back(p.second);
	T.push_back(matrix_values[i]);
        }
  //std::cout<< "values "<<matrix_values[i]<<std::endl;

	i++;
      }
    
  
  return Matrix::sparse(dim1, dim2, nint(rows), nint(cols), ndou(T));
  }
};