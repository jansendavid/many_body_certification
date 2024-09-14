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

using namespace mosek::fusion;
using namespace monty;
using string_pair=std::pair<std::string, std::string>;

// file with spin hamiltonians for translation invariant system
std::pair<int,int> get_sec(op_vec op)
{
  // computes the sector of a given vector of operators
  // NOTE: all are of the form vec{S_i}vec{S_j} in the Hamiltonian is needed. For e.g., TFI, small modifications must be made
  int sxy=1;
  int syz=1;
  for(auto a: op)
    {
      if(a.dir_=="x"){sxy*=-1;}
      if(a.dir_=="y"){syz*=-1;
	sxy*=-1;}
      if(a.dir_=="z"){syz*=-1;}
      
    }
  return std::pair<int,int>(sxy, syz);
}
void add_state(basis_structure& states, op_vec op, std::map<std::pair<int,int>, int> map_sec)
{
  // adds a state to a basis
      auto [fac, nf] =get_normal_form(op);
	     auto sign=get_sec( nf);
	     if(nf.size()>0)
	       {
	     states.at(map_sec.at(sign)).push_back(nf);
	       }
	     
	    

  return;}

void add_state_symm(basis_structure& states, op_vec op, std::map<std::pair<int,int>, int> map_sec, int L, std::string permuts)
{
  // adds a state to a basis
      auto [fac, nf] =get_normal_form(op);
      auto sign=get_sec( nf);
      if(nf.size()>0)
	{
	  auto all_t=generate_all_translations(nf, L);
	  bool found=false;
	  
	     	 for(auto op_t: all_t)
	   {
	     
	     
      // 	     // generate all y translations, inc=1
	      auto all_ty=generate_all_translations_y(op_t, L,1);
	     for(auto op_ty: all_ty)
	      {
      // 	     // generate all permutations
	   // 	 std::vector<op_vec> all_p;	 
	   // 	 if(permuts=="xyz")
	   // 	   {
	   //    all_p=generate_all_permutations_xyz(op_ty);
	   // 	   }
	   // 	 else if(permuts=="xy")
	   // 	   {
	   // 	     all_p=generate_all_permutations_xy(op_ty);
	   // 	   }

	   //   for(auto op_p: all_p)
	   // { 
	     
	     int cnt = count(states.at(map_sec.at(sign)).begin(), states.at(map_sec.at(sign)).end(), op_ty);
	     if(cnt>0){
	       found=true;
	       //}
	        }
	       }
       	   }
	     if(not found)
	       {   states.at(map_sec.at(sign)).push_back(nf);}
       	   }
	     
	    

  return;}





// HAMILTONIANS WHEN USING TRANSLATION SYMMETRY 
// NOTE: all are of the form vec{S_i}vec{S_j}

Expression::t define_xxz1d( std::map<std::string, mom_ref> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J, double Delta)
{

  std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
  auto ham=Expr::constTerm(0.);
  for(auto term:dirs)
    {
      op_vec v={spin_op(term.first, {0}),spin_op(term.second, {1})};
      auto [fac, nf] =get_normal_form(v);
      auto [ key,coeff_map]=map.at(print_op(nf));  
      auto el=refs.at(key);
      double coeff=J;
      if(term==string_pair("z","z")){coeff=Delta;}
      
      if(std::abs(fac.imag())>1e-9 or std::abs(coeff_map.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      ham=Expr::add(ham,Expr::mul(coeff*fac.real()*coeff_map.real()/4., el.var_));
    }


  


  return ham;
}

void define_xxz2d_dual( std::map<std::string, int> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J, double Delta, int Ly, int Lx,Expression::t& second_block)
{
    std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
  auto ham=Expr::constTerm(0.);
  std::vector<int> rows;
  std::vector<double> vals;
  for(auto term:dirs)
    {
      //for(int i=0; i<Ly; i++)
      //{
	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
	  //op_vec v_p={spin_op(term.first, {std::min(i,(i+1)%Ly),0}, Lx),spin_op(term.second, {std::max(i,(i+1)%Ly),0}, Lx)};
	  op_vec v_p={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {1,0}, Lx)};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
      auto el_p=refs.at(key_p);
      double coeff=J;
      if(term==string_pair("z","z")){coeff=Delta;}
      
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      {

	//  second_block[el_p]=Expr::add(second_block[el_p],-1.*Ly*coeff*fac_p.real()*coeff_map_p.real()/4.);
	  rows.push_back(el_p);
	  vals.push_back(-1.*Ly*coeff*fac_p.real()*coeff_map_p.real()/4.);
      }
      // block[el_p]=Expr::add(block[el_p],-1.*Ly*coeff*fac_p.real()*coeff_map_p.real()/4. );
      //      ham=Expr::add(ham,Expr::mul(Ly*coeff*fac_p.real()*coeff_map_p.real()/4., el_p.var_));

      // transverse part part S_(0,0)S_(0,1),S_(1,0)S_(1,1)... etc
      //op_vec v_t={spin_op(term.first, {i,0}, Lx),spin_op(term.second, {i,1}, Lx)};
      op_vec v_t={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {0,1}, Lx)};
      auto [fac_t, nf_t] =get_normal_form(v_t);
      auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
      auto el_t=refs.at(key_t);
      coeff=J;
      if(term==string_pair("z","z")){coeff=Delta;}
      
      if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
    
      {
  
	//	  second_block[el_t]=Expr::add(second_block[el_t],-1.*Ly*coeff*fac_t.real()*coeff_map_t.real()/4);
	  rows.push_back(el_t);
	  vals.push_back(-1.*Ly*coeff*fac_t.real()*coeff_map_t.real()/4);
      }
    
      //   block[el_t]=Expr::add(block[el_t], -1.*Ly*coeff*fac_t.real()*coeff_map_t.real()/4);
      //      ham=Expr::add(ham,Expr::mul(Ly*coeff*fac_t.real()*coeff_map_t.real()/4., el_t.var_));
    }
std::vector<int> cols(rows.size(),0);
 Matrix::t Alpha_e_1  = Matrix::sparse(refs.size(),1, nint(rows), nint(cols), ndou(vals)); 
 second_block=Expr::add(second_block,  Alpha_e_1);
  return;
    

}


Expression::t define_xxz2d( std::map<std::string, mom_ref> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J, double Delta, int Ly, int Lx)
{
    std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
  auto ham=Expr::constTerm(0.);
  for(auto term:dirs)
    {
      //for(int i=0; i<Ly; i++)
      //{
	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
	  //op_vec v_p={spin_op(term.first, {std::min(i,(i+1)%Ly),0}, Lx),spin_op(term.second, {std::max(i,(i+1)%Ly),0}, Lx)};
	  op_vec v_p={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {1,0}, Lx)};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
      auto el_p=refs.at(key_p);
      double coeff=J;
      if(term==string_pair("z","z")){coeff=Delta;}
      
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      ham=Expr::add(ham,Expr::mul(Ly*coeff*fac_p.real()*coeff_map_p.real()/4., el_p.var_));

      // transverse part part S_(0,0)S_(0,1),S_(1,0)S_(1,1)... etc
      //op_vec v_t={spin_op(term.first, {i,0}, Lx),spin_op(term.second, {i,1}, Lx)};
      op_vec v_t={spin_op(term.first, {0,0}, Lx),spin_op(term.second, {0,1}, Lx)};
      auto [fac_t, nf_t] =get_normal_form(v_t);
      auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
      auto el_t=refs.at(key_t);
      coeff=J;
      if(term==string_pair("z","z")){coeff=Delta;}
      
      if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      ham=Expr::add(ham,Expr::mul(Ly*coeff*fac_t.real()*coeff_map_t.real()/4., el_t.var_));
    }
  //}

  return ham;
    

}

Expression::t define_J1J2_2d( std::map<std::string, mom_ref> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J1, double J2, int Ly, int Lx)
{
    std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
  auto ham=Expr::constTerm(0.);
  for(auto term:dirs)
    {
      for(int i=0; i<Ly; i++)
	{
	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
      op_vec v_p={spin_op(term.first, {std::min(i,(i+1)%Ly),0}, Lx),spin_op(term.second, {std::max(i,(i+1)%Ly),0}, Lx)};
      auto [fac_p, nf_p] =get_normal_form(v_p);
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
      auto el_p=refs.at(key_p);
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      ham=Expr::add(ham,Expr::mul(J1*fac_p.real()*coeff_map_p.real()/4., el_p.var_));
      op_vec v_t={spin_op(term.first, {i,0}, Lx),spin_op(term.second, {i,1}, Lx)};
      auto [fac_t, nf_t] =get_normal_form(v_t);
      auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
      auto el_t=refs.at(key_t);

      
      
      if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      ham=Expr::add(ham,Expr::mul(J1*fac_t.real()*coeff_map_t.real()/4., el_t.var_));

    }
    }
    for(auto term:dirs)
    {
      for(int i=0; i<Ly; i++)
	{
	  // parallel part S_(0,0)S_(1,0),S_(1,0)S_(2,0)... etc
      op_vec v_p={spin_op(term.first, {std::min(i,(i+2)%Ly),0}, Lx),spin_op(term.second, {std::max(i,(i+2)%Ly),0}, Lx)};
      auto [fac_p, nf_p] =get_normal_form(v_p);
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
      auto el_p=refs.at(key_p);
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      ham=Expr::add(ham,Expr::mul(J2*fac_p.real()*coeff_map_p.real()/4., el_p.var_));

      
      op_vec v_t={spin_op(term.first, {i,0}, Lx),spin_op(term.second, {i,2}, Lx)};
      auto [fac_t, nf_t] =get_normal_form(v_t);
      auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
      auto el_t=refs.at(key_t);

      
      
      if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      ham=Expr::add(ham,Expr::mul(J2*fac_t.real()*coeff_map_t.real()/4., el_t.var_));

    }
    }

  return ham;
    

}
