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
using namespace mosek::fusion;
using namespace monty;
using mat_type=Eigen::MatrixXcd;
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
std::vector<op_vec> load_basis_from_file(std::string filename, int Lx)
{
  std::vector<op_vec> basis;
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
		state.push_back(spin_op(sa.substr(j,1), {numbs[i],numbs[i+1]}, Lx));
	
	      i+=2;
	    }
	  basis.push_back(state);
        }
        
        // Close the file object.
        new_file.close(); 
    }
  return basis;

}
std::string find_state(op_vec nf,std::map<std::string, mom_ref> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> mp, int Lx )
{
  // pass in normal form

 
	  auto all_t=generate_all_translations(nf, Lx);
	  bool found=false;
	  
	     	 for(auto op_t: all_t)
	   {
	     
	     
      // 	     // generate all y translations, inc=1
	      auto all_ty=generate_all_translations_y(op_t, Lx,1);
	      std::cout<< "start"<<std::endl;
	      
	     for(auto op_ty: all_ty)
	      {
		
		auto key=print_op(op_ty);
		auto it = mp.find(key);
		std::cout<< print_op(op_ty)<<std::endl;	
		    if (it != mp.end())
		      {
			return key;
		      }
		    else{
		    
		    }
     
       	        }
       	       }
		 std::cout<< "RDM not found "<<print_op(nf)<< std::endl;
		 std::string s="";
		 return s;   

      
}
    //xt::xarray<std::complex<double>, xt::layout_type::dynamic>;
void rdm(int degree,std::map<std::string, mom_ref> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, int Ly, int Lx)
{

  //mat_type sigz({2,2});
  
   Eigen::MatrixXcd sigz=Eigen::MatrixXcd::Zero(2,2);
  sigz(0,0) = std::complex<double>(1.,0);
  sigz(1,1) = std::complex<double>(-1,0);
   Eigen::MatrixXcd sigx=Eigen::MatrixXcd::Zero(2,2);
   //mat_type sigx({2,2});
  sigx(0,1) = std::complex<double>(1.,0);
  sigx(1,0) = std::complex<double>(1,0);
   Eigen::MatrixXcd sigy=Eigen::MatrixXcd::Zero(2,2);
   //  mat_type sigy({2,2});
  sigy(0,1) = std::complex<double>(0,-1.);
  sigy(1,0) = std::complex<double>(0,1);
  // std::cout<< sigy<<std::endl;
   Eigen::MatrixXcd unit=Eigen::MatrixXcd::Zero(2,2);
   //mat_type unit({2,2});
  unit(0,0) = std::complex<double>(1,0);
  unit(1,1) = std::complex<double>(1,0);
  // std::cout<< sigy<<std::endl;
  
  std::vector<std::string> terms;
  std::map<std::string, mat_type> sigma_map;
  sigma_map.insert({"1",unit});
  sigma_map.insert({"x",sigx});
  sigma_map.insert({"y",sigy});
  sigma_map.insert({"z",sigz});
  
  auto dirs=std::vector<std::string>{"1","x", "y", "z"};
   std::vector<std::string> tots;

       for(auto d1: dirs)
	 {

	   if(degree==1)
	     {tots.push_back(d1); continue;}
	          for(auto d2: dirs)
	 {
	   	   if(degree==2)
	     {tots.push_back(d1+d2); continue;}
	   	          for(auto d3: dirs)
	 {
	   	   	   if(degree==3)
	     {tots.push_back(d1+d2+d3); continue;}
			   for(auto d4: dirs)
	 {
	   	   	   if(degree==4)
	     {tots.push_back(d1+d2+d3+d4); continue;}
	 }
	 }
	 }
	 }
       int dim=std::pow(2, degree);
       auto arr_zero = std::make_shared<ndarray<double,2>>(shape(dim,dim));
       for(int l=0; l<dim;l++)
	 {
       for(int n=0; n<dim;n++)
	 {
	   (*arr_zero)(l,n)=0;
	 }
	 }
       auto Mat_fusion_zero=Matrix::dense( arr_zero);
       auto expr_imag=Expr::constTerm(Mat_fusion_zero);
       auto expr_real=Expr::constTerm(Mat_fusion_zero);
       for(auto t: tots)
	 {
	   mat_type mat;
	   op_vec state;
	   
	   for(int i=0; i<t.size(); i++)
	     {
	       std::string key=t.substr(i,i+1);
	       if(key!="1")
		 {
	       state.push_back(spin_op(key, {0,i}, Lx));
		 }
	       if(i==0)
		 {
		   mat=sigma_map[key];
		   //  std::cout<<key<< " and "<< t[i]<< " "<<sigma_map[key]<<std::endl;
		 }
	       else
		 {
		   
		   mat=Eigen::KroneckerProduct(mat,sigma_map[key]).eval();
		   
	     }

	 }
	  
	   auto mat_r=mat.real();
	   auto a = std::make_shared<ndarray<double,2>>(shape(int(mat.rows()),int(mat.rows())));
	   
	   for(int l=0; l<int(mat.rows()); l++)
	     {

for(int n=0; n<int(mat.rows()); n++)
  {
    (*a)(l,n)=mat_r(l,n);
   }
	     }
	   if(state.size()>0)
	     {
      auto [fac, nf] =get_normal_form(state);
      
      auto state_key=find_state(nf,refs,map,Lx );
	         auto [ key_p,coeff_map_p]=map.at(state_key);
	     }
 			   auto Mat_fusion=Matrix::dense( a);
			   //auto expr_real=Expr::constTerm(Mat_fusion);
			   expr_real=Expr::add( expr_real, Mat_fusion);
	   
	 }
  
  return;
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
