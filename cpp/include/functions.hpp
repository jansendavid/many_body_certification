#pragma once
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xnpy.hpp>
using namespace mosek::fusion;
using namespace monty;
void store_matrix_to_numpy_eff(const std::string &filename, size_t dim, Variable::t X){
  // Matrix in fusion
   auto Mat_fusion=Matrix::dense(2*dim, 2*dim, X->level());
  std::vector<size_t> shape = { dim, dim };
  // Matrix in xtensor
  xt::xarray<std::complex<double>, xt::layout_type::dynamic> M(shape, xt::layout_type::row_major);

  for(int i=0; i<dim; i++)
       {
       for(int j=0; j<dim;j++)
    {
      // imag X_3 -X_3^T
      // top left is X_3
      M(i,j)=std::complex(Mat_fusion->get(i,j)+Mat_fusion->get(dim+i,dim+j), Mat_fusion->get(j,dim+i)-Mat_fusion->get(i,dim+j));
     
    }

    }

  xt::dump_npy(filename, M);
  return;}



void store_matrix_to_numpy(const std::string &filename, size_t dim, Variable::t X){
  // Matrix in fusion
   auto Mat_fusion=Matrix::dense(2*dim, 2*dim, X->level());
  std::vector<size_t> shape = { dim, dim };
  // Matrix in xtensor
  xt::xarray<std::complex<double>, xt::layout_type::dynamic> M(shape, xt::layout_type::row_major);

  for(int i=0; i<dim; i++)
       {
       for(int j=0; j<dim;j++)
    {

      M(i,j)=std::complex(Mat_fusion->get(i,j), Mat_fusion->get(i,dim+j));
     
    }

    }

  xt::dump_npy(filename, M);
  return;}





// struct mat_el{
//   mat_el(int ind_1, int ind_2, double coeff):i1_(ind_1), i2_(ind_2), coeff_(coeff){}
//   int i1_;
//   int i2_;
//   double coeff_;
//   mat_el()=default;
// };
// struct mat_el_ref{
//   mat_el_ref(Variable::t var, int i1, int i2,double coeff):var_(var),i1_(i1), i2_(i2), coeff_(coeff){}
//   Variable::t var_;
//   double coeff_;
//    int i1_;
//   int i2_;
//   mat_el_ref()=default;
// };
// std::ostream& operator<<(std::ostream& os, const mat_el_ref& op)
// {
//   os <<"("<<op.i1_ << ","<<op.i2_<< ") coeff: "<< op.coeff_<<std::endl;
//     return os;
// }



// std::ostream& operator<<(std::ostream& os, const mat_el& op)
// {
//   os <<"("<<op.i1_ << ","<<op.i2_<< ") coeff: "<< op.coeff_<<std::endl;
//     return os;
// }
