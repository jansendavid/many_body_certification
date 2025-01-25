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

std::vector<double> define_xxz2d_sos( std::map<std::string, int> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J, double Delta, int Ly, int Lx)
{
    std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
  std::vector<double> vals(refs.size(), 0);

  for(auto term:dirs)
    {
   
     
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
std::vector<double> define_J1J22d_sos( std::map<std::string, int> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J1,double J2, int Ly, int Lx)
{
    std::vector<string_pair> dirs{string_pair("x","x"),string_pair("z","z"),string_pair("y","y")};
  std::vector<double> vals(refs.size(), 0);
// J1 nearest neighbour interaction
  for(auto term:dirs)
    {
   
     
	  op_vec v_p={spin_op(term.first, {0,0}, {Lx, Ly}),spin_op(term.second, {1,0}, {Lx, Ly})};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
     
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
    
      auto el_p=refs.at(key_p);
 
      double coeff=J1;
 
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

      if(std::abs(fac_t.imag())>1e-9 or std::abs(coeff_map_t.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
    
      {
  
vals[el_t]+=coeff*fac_t.real()*coeff_map_t.real()/4.;
	
      }
    
    
    }
// J2 next nearest neighbour interaction
 for(auto term:dirs)
    {
   
     
	  op_vec v_p={spin_op(term.first, {0,0}, {Lx, Ly}),spin_op(term.second, {1,1}, {Lx, Ly})};	  
      auto [fac_p, nf_p] =get_normal_form(v_p);
     
      auto [ key_p,coeff_map_p]=map.at(print_op(nf_p));  
    
      auto el_p=refs.at(key_p);
 
      double coeff=J2;
 
      if(std::abs(fac_p.imag())>1e-9 or std::abs(coeff_map_p.imag())>1e-9){
	std::cout<< "error: Hamiltonian contains complex elements "<<std::endl;
      }
      {

  vals[el_p]+=coeff*fac_p.real()*coeff_map_p.real()/4.;
       }
   
      op_vec v_t={spin_op(term.first, {0,0}, {Lx, Ly}),spin_op(term.second, {1,Lx-1}, {Lx, Ly})};
      auto [fac_t, nf_t] =get_normal_form(v_t);
      auto [ key_t,coeff_map_t]=map.at(print_op(nf_t));  
      auto el_t=refs.at(key_t);

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
// J perpendicular
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
for(int i=0; i<layers; i++)
{
        for(auto term:dirs)
    {

    {
	  op_vec v_p={spin_op(term.first, {i,0,0}, {layers, Lx, Ly}),spin_op(term.second, {(i+1)%layers,0,1}, {layers, Lx, Ly})};	  
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
	  op_vec v_p={spin_op(term.first, {i,0,0}, {layers, Lx, Ly}),spin_op(term.second, {(i+1)%layers,1,0}, {layers, Lx, Ly})};	  
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
}

  return vals; 
    

}
