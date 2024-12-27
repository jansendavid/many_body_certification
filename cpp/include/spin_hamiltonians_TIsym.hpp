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
#include"symmetries.hpp"
using namespace mosek::fusion;
using namespace monty;
using string_pair=std::pair<std::string, std::string>;

// file with spin hamiltonians for translation invariant system


// HAMILTONIANS WHEN USING TRANSLATION SYMMETRY 
// NOTE: all are of the form vec{S_i}vec{S_j}

// Expression::t define_xxz1d( std::map<std::string, mom_ref> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J, double Delta)
// {

//   std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
//   auto ham=Expr::constTerm(0.);
//   for(auto term:dirs)
//     {
//       op_vec v={spin_op(term.first, {0}),spin_op(term.second, {1})};
//       auto [fac, nf] =get_normal_form(v);
//       auto [ key,coeff_map]=map.at(print_op(nf));  
//       auto el=refs.at(key);
//       double coeff=J;
//       if(term==string_pair("z","z")){coeff=Delta;}
      
//       if(std::abs(fac.imag())>1e-9 or std::abs(coeff_map.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       ham=Expr::add(ham,Expr::mul(coeff*fac.real()*coeff_map.real()/4., el.var_));
//     }


  


//   return ham;
// }

//  Matrix::t define_xxz2d_dual( std::map<std::string, int> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J, double Delta, int Ly, int Lx)
// {
//     std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
//   auto ham=Expr::constTerm(0.);
//   std::vector<int> rows;
//   std::vector<double> vals;

//   for(auto term:dirs)
//     {
   
//       //for(int i=0; i<Ly; i++)
//       //{
// 	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
// 	  //op_vec v_p={spin_op(term.first, {std::min(i,(i+1)%Ly),0}, Lx),spin_op(term.second, {std::max(i,(i+1)%Ly),0}, Lx)};
// 	  op_vec v_p={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {1,0}, Lx)};	  
//       auto [fac_p, nf_p] =get_normal_form(v_p);
     
//       auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
    
//       auto el_p=refs.at(key_p);
 
//       double coeff=J;
//       if(term==string_pair("z","z")){coeff=Delta;}
      
//       if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       {


// 	auto index_1=getIndex(rows,  el_p);
  
// 	if(index_1<0)
// 	  {
// 	  rows.push_back(el_p);
// 	  vals.push_back(coeff*fac_p.real()*coeff_map_p.real()/4.);
// 	  }
// 	else{
// 	  vals[index_1]+=coeff*fac_p.real()*coeff_map_p.real()/4.;
// 	}
//        }
   
//       op_vec v_t={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {0,1}, Lx)};
//       auto [fac_t, nf_t] =get_normal_form(v_t);
//       auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
//       auto el_t=refs.at(key_t);
//       coeff=J;
//       if(term==string_pair("z","z")){coeff=Delta;}
      
//       if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
    
//       {
  

// 		auto index_1=getIndex(rows,  el_t);
// 	if(index_1<0)
// 	  {
//   rows.push_back(el_t);
// 	  vals.push_back(coeff*fac_t.real()*coeff_map_t.real()/4);
// 	  }
// 	else{
// 	  vals[index_1]+=coeff*fac_t.real()*coeff_map_t.real()/4.;
// 	}

	 
	
//       }
    
    
//     }

// std::vector<int> cols(rows.size(),0);
//  Matrix::t Alpha_e_1  = Matrix::sparse(refs.size(),1, nint(rows), nint(cols), ndou(vals)); 
 
//   return Alpha_e_1; 
    

// }

std::vector<double> define_xxz2d_sos( std::map<std::string, int> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J, double Delta, int Ly, int Lx)
{
    std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
  std::vector<double> vals(refs.size(), 0);

  for(auto term:dirs)
    {
   
      //for(int i=0; i<Ly; i++)
      //{
	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
	  //op_vec v_p={spin_op(term.first, {std::min(i,(i+1)%Ly),0}, Lx),spin_op(term.second, {std::max(i,(i+1)%Ly),0}, Lx)};
	  op_vec v_p={spin_op(term.first, {0,0}, {Lx, Ly}),spin_op(term.second, {1,0}, {Lx, Ly})};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
     
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
    
      auto el_p=refs.at(key_p);
 
      double coeff=J;
      if(term==string_pair("z","z")){coeff=Delta;}
      
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      {

  vals[el_p]+=coeff*fac_p.real()*coeff_map_p.real()/4.;
       }
   
      op_vec v_t={spin_op(term.first, {0,0}, {Lx, Ly}),spin_op(term.second, {0,1}, {Lx, Ly})};
      auto [fac_t, nf_t] =get_normal_form(v_t);
      auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
      auto el_t=refs.at(key_t);
      coeff=J;
      if(term==string_pair("z","z")){coeff=Delta;}
      
      if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
    
      {
  
vals[el_t]+=coeff*fac_t.real()*coeff_map_t.real()/4.;
	
      }
    
    
    }

  return vals; 
    

}
std::vector<double> define_heisenberg_bilayer_sos( std::map<std::string, int> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J_perpendicular,double J_parallel, double J_x,int layers, int Ly, int Lx)
{
    std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
  std::vector<double> vals(refs.size(), 0);
// J_parallel terms
  for(auto term:dirs)
    {
   for(int i=0; i<2; i++)
   {
    {
	  op_vec v_p={spin_op(term.first, {i,0,0}, {layers, Lx, Ly}),spin_op(term.second, {i,1,0}, {layers, Lx, Ly})};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
     
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
    
      auto el_p=refs.at(key_p);
 

      
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      {

  vals[el_p]+=J_parallel*fac_p.real()*coeff_map_p.real()/4.;
       }
    }
     {
	  op_vec v_p={spin_op(term.first, {i,0,0}, {layers, Lx, Ly}),spin_op(term.second, {i,0,1}, {layers, Lx, Ly})};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
     
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
    
      auto el_p=refs.at(key_p);
 

      
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      {

  vals[el_p]+=J_parallel*fac_p.real()*coeff_map_p.real()/4.;
       }
    }
   }
 
    
    }
// J x
    for(auto term:dirs)
    {

    {
	  op_vec v_p={spin_op(term.first, {0,0,0}, {layers, Lx, Ly}),spin_op(term.second, {1,0,0}, {layers, Lx, Ly})};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
     
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
    
      auto el_p=refs.at(key_p);
 

      
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      {

  vals[el_p]+=J_perpendicular*fac_p.real()*coeff_map_p.real()/4.;
       }
    }
   
   
 
    
    }
// J_x
        for(auto term:dirs)
    {

    {
	  op_vec v_p={spin_op(term.first, {0,0,0}, {layers, Lx, Ly}),spin_op(term.second, {1,0,1}, {layers, Lx, Ly})};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
     
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
    
      auto el_p=refs.at(key_p);
 

      
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      {

  vals[el_p]+=J_x*fac_p.real()*coeff_map_p.real()/4.;
       }
    }
   
   
    {
	  op_vec v_p={spin_op(term.first, {0,0,0}, {layers, Lx, Ly}),spin_op(term.second, {1,1,0}, {layers, Lx, Ly})};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
     
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
    
      auto el_p=refs.at(key_p);
 

      
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      {

  vals[el_p]+=J_x*fac_p.real()*coeff_map_p.real()/4.;
       }
    }
   
    
    }

  return vals; 
    

}

// Expression::t define_xxz2d( std::map<std::string, mom_ref> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J, double Delta, int Ly, int Lx)
// {
//     std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
//   auto ham=Expr::constTerm(0.);
//   for(auto term:dirs)
//     {
//       //for(int i=0; i<Ly; i++)
//       //{
// 	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
// 	  //op_vec v_p={spin_op(term.first, {std::min(i,(i+1)%Ly),0}, Lx),spin_op(term.second, {std::max(i,(i+1)%Ly),0}, Lx)};
// 	  op_vec v_p={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {1,0}, Lx)};	  
//       auto [fac_p, nf_p] =get_normal_form(v_p);
//       auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
//       auto el_p=refs.at(key_p);
//       double coeff=J;
//       if(term==string_pair("z","z")){coeff=Delta;}
      
//       if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       std::cout<< " adding "<< el_p.i1_<<std::endl;
//       ham=Expr::add(ham,Expr::mul(Ly*coeff*fac_p.real()*coeff_map_p.real()/4., el_p.var_));

//       // transverse part part S_(0,0)S_(0,1),S_(1,0)S_(1,1)... etc
//       //op_vec v_t={spin_op(term.first, {i,0}, Lx),spin_op(term.second, {i,1}, Lx)};
//       op_vec v_t={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {0,1}, Lx)};
//       auto [fac_t, nf_t] =get_normal_form(v_t);
//       auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
//       auto el_t=refs.at(key_t);
//       coeff=J;
//       if(term==string_pair("z","z")){coeff=Delta;}
      
//       if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       std::cout<< " adding "<< el_t.i1_<<std::endl;
//       ham=Expr::add(ham,Expr::mul(Ly*coeff*fac_t.real()*coeff_map_t.real()/4., el_t.var_));
//     }
//   //}

//   return ham;
    

// }
// Matrix::t define_J1J2_2d_dual( std::map<std::string, int> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J1, double J2, int Ly, int Lx)
// {
//     std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
//   auto ham=Expr::constTerm(0.);
//   std::vector<int> rows;
//   std::vector<double> vals;
//   for(auto term:dirs)
//     {
      
// 	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
//       op_vec v_p={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {1,0}, Lx)};
//       auto [fac_p, nf_p] =get_normal_form(v_p);
//       auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
//       auto el_p=refs.at(key_p);
//       if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }

//       {
// 	auto index_1=getIndex(rows,  el_p);
// 	if(index_1<0)
// 	  {
// 	  rows.push_back(el_p);
// 	  vals.push_back(Ly*J1*fac_p.real()*coeff_map_p.real()/4.);
// 	  }
// 	else{
// 	  vals[index_1]+=Ly*J1*fac_p.real()*coeff_map_p.real()/4.;
// 	}
//       }

//       op_vec v_t={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {0,1}, Lx)};
//       auto [fac_t, nf_t] =get_normal_form(v_t);
//       auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
//       auto el_t=refs.at(key_t);

      
      
//       if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       {

// 		auto index_1=getIndex(rows,  el_t);
// 	if(index_1<0)
// 	  {
//   rows.push_back(el_t);
// 	  vals.push_back(Ly*J1*fac_t.real()*coeff_map_t.real()/4.);
// 	  }
// 	else{
// 	  vals[index_1]+=Ly*J1*fac_t.real()*coeff_map_t.real()/4.;
// 	}
//       }
      
//     }
//     for(auto term:dirs)
//     {

// 	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
//       op_vec v_p={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {2,0}, Lx)};
//       auto [fac_p, nf_p] =get_normal_form(v_p);
//       auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
//       auto el_p=refs.at(key_p);
//       if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//             {
// 	auto index_1=getIndex(rows,  el_p);
// 	if(index_1<0)
// 	  {
// 	  rows.push_back(el_p);
// 	  vals.push_back(Ly*J2*fac_p.real()*coeff_map_p.real()/4.);
// 	  }
// 	else{
// 	  vals[index_1]+=Ly*J2*fac_p.real()*coeff_map_p.real()/4.;
// 	}
// 	    }
//       op_vec v_t={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {0,2}, Lx)};
//       auto [fac_t, nf_t] =get_normal_form(v_t);
//       auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
//       auto el_t=refs.at(key_t);      
//       if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       {
//       		auto index_1=getIndex(rows,  el_t);
// 	if(index_1<0)
// 	  {
//   rows.push_back(el_t);
// 	  vals.push_back(Ly*J2*fac_t.real()*coeff_map_t.real()/4.);
// 	  }
// 	else{
// 	  vals[index_1]+=Ly*J2*fac_t.real()*coeff_map_t.real()/4.;
// 	}
//       }

// }
  


// std::vector<int> cols(rows.size(),0);
//  Matrix::t Alpha_e_1  = Matrix::sparse(refs.size(),1, nint(rows), nint(cols), ndou(vals)); 
 
//   return Alpha_e_1; 
    

// }



// Expression::t define_J1J2_2d( std::map<std::string, mom_ref> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J1, double J2, int Ly, int Lx)
// {
//     std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
//   auto ham=Expr::constTerm(0.);
//   for(auto term:dirs)
//     {
//       for(int i=0; i<Ly; i++)
// 	{
// 	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
//       op_vec v_p={spin_op(term.first, {std::min(i,(i+1)%Ly),0}, Lx),spin_op(term.second, {std::max(i,(i+1)%Ly),0}, Lx)};
//       auto [fac_p, nf_p] =get_normal_form(v_p);
//       auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
//       auto el_p=refs.at(key_p);
//       if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       ham=Expr::add(ham,Expr::mul(J1*fac_p.real()*coeff_map_p.real()/4., el_p.var_));
//       op_vec v_t={spin_op(term.first, {i,0}, Lx),spin_op(term.second, {i,1}, Lx)};
//       auto [fac_t, nf_t] =get_normal_form(v_t);
//       auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
//       auto el_t=refs.at(key_t);

      
      
//       if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       ham=Expr::add(ham,Expr::mul(J1*fac_t.real()*coeff_map_t.real()/4., el_t.var_));

//     }
//     }
//     for(auto term:dirs)
//     {
//       for(int i=0; i<Ly; i++)
// 	{
// 	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
//       op_vec v_p={spin_op(term.first, {std::min(i,(i+2)%Ly),0}, Lx),spin_op(term.second, {std::max(i,(i+2)%Ly),0}, Lx)};
//       auto [fac_p, nf_p] =get_normal_form(v_p);
//       auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
//       auto el_p=refs.at(key_p);
//       if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       ham=Expr::add(ham,Expr::mul(J2*fac_p.real()*coeff_map_p.real()/4., el_p.var_));

      
//       op_vec v_t={spin_op(term.first, {i,0}, Lx),spin_op(term.second, {i,2}, Lx)};
//       auto [fac_t, nf_t] =get_normal_form(v_t);
//       auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
//       auto el_t=refs.at(key_t);

      
      
//       if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
// 	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
//       }
//       ham=Expr::add(ham,Expr::mul(J2*fac_t.real()*coeff_map_t.real()/4., el_t.var_));

//     }
//     }

//   return ham;
    

// }

