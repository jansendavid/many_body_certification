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
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include"spins.hpp"
using namespace mosek::fusion;
using namespace monty;

// function to extract all integers
std::vector<int> extractIntegerWords(std::string str)
{
  
  std::stringstream ss;
  
    /* Storing the whole string into string stream */
    ss << str;
 
    /* Running loop till the end of the stream */
    std::string temp;
    int found;
    std::vector<int> numbers;
    while (!ss.eof()) {
 
        /* extracting word by word from stream */
        ss >> temp; 
        /* Checking the given word is integer or not */
        if (std::stringstream(temp) >> found)
	  {

	    numbers.push_back(found);
	  }
 
        /* To save from space at the end of string */
        temp = "";
    }
    return numbers;
}
std::vector<op_vec> load_basis_from_file(std::string filename, int Lx, bool bilayer=false)

{
   std::vector<op_vec> basis;
std::vector<int> offset;
int shift=2;
if(bilayer)
{
  offset={2,Lx,Lx};
  shift=3;
}   
else{
  offset={Lx,Lx};
}
  std::fstream new_file;
  new_file.open(filename, std::ios::in); 
    
    // Checking whether the file is open.
    if (new_file.is_open()) { 
      std::string sa;
        // Read data from the file object and put it into a string.
        while (getline(new_file, sa)) { 
            // Print the data of the string.
	  
	  std::replace( sa.begin(), sa.end(), ',', ' ');
	  std::replace( sa.begin(), sa.end(), '(', ' ');
	  std::replace( sa.begin(), sa.end(), ')', ' ');
	  auto sa_copy=sa;//std::strcpy(sa.data(), sa.data());
	  sa.erase(std::remove(sa.begin(), sa.end(), ' '), sa.end());
	  auto numbs=extractIntegerWords(sa_copy);
	 
	  sa.erase(std::remove_if(std::begin(sa), std::end(sa), [](auto ch) { return std::isdigit(ch); }),sa.end());
	  op_vec state;
	  int i=0; // index thar gives number
	  
	  for(int j=0; j<sa.size(); j++){
      std::vector<int> indices;
      if(bilayer){
        indices={numbs[i],numbs[i+1], numbs[i+3]};
      }
      else{
indices={numbs[i],numbs[i+1]};

      }
		state.push_back(spin_op(sa.substr(j,1), indices, offset));
	
	      i+=2;
	    }
	  basis.push_back(state);
        }
        
        // Close the file object.
        new_file.close(); 
    }
   return basis;

 }

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
