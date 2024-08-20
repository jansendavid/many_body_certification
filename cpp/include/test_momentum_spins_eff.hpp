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
#include"complex_momentum_symm_eff.hpp"
#include"spin_hamiltonians_TIsym.hpp"
using namespace mosek::fusion;
using namespace monty;

// file with test functions for using translation symmetry
void test_translation()
{
  int L=3;
   
    op_vec v1={spin_op("x", {0}),spin_op("y", {1})};

    //std::cout<<print_op(v1)<<std::endl;
    // for(int i=0; i<=L;i++)
    //   {
    // 	auto O=translation(v1, i, L);
    // 	//	std::cout<<print_op(O)<<std::endl;
    //   }
    auto all_T=generate_all_translations(v1, L);
    for(auto& a: all_T)
      {
	std::cout<<print_op(a)<<std::endl;
      }
    
}
void test_permutations()
{
   int L=5;
   
   op_vec v1={spin_op("x", {0}),spin_op("y", {1}),spin_op("x", {3}),spin_op("z", {4})};

  //   //std::cout<<print_op(v1)<<std::endl;
  //   // for(int i=0; i<=L;i++)
  //   //   {
  //   // 	auto O=translation(v1, i, L);
  //   // 	//	std::cout<<print_op(O)<<std::endl;
  //   //   }
     auto all_P=generate_all_permutations_xyz(v1);
     for(auto& a: all_P)
       {
   	std::cout<<print_op(a)<<std::endl;
       }
    
}

void test_translation_2d()
{
  int Lx=5;
  int Ly=5;
   
  op_vec v1={spin_op("z", {3, 0}, Lx),spin_op("z", {3,3}, Lx)};

    //std::cout<<print_op(v1)<<std::endl;
    // for(int i=0; i<=L;i++)
    //   {
    // 	auto O=translation(v1, i, L);
    // 	//	std::cout<<print_op(O)<<std::endl;
    //   }
    auto all_T=generate_all_translations(v1, Lx);
    for(auto& a: all_T)
      {
	std::cout<<print_op(a)<<std::endl;
      }
    
}
void test_single_block()
{
  // testing translation invariance for single symmetrie sector
  int L=5;
  basis_structure states;
  std::vector<op_vec> v_block_0;
  std::vector<op_vec> v_block_1;
  std::vector<op_vec> v_block_2;
  std::vector<op_vec> v_block_3;
  states.insert({0, v_block_0});     
  std::map<std::pair<int,int>, int> map_sec;
  map_sec.insert({std::pair<int,int>(1,1), 0});
  map_sec.insert({std::pair<int,int>(1,-1), 0});
  map_sec.insert({std::pair<int,int>(-1,1), 0});
  map_sec.insert({std::pair<int,int>(-1,-1), 0});
 
       std::vector<std::string> dirs={"x", "y", "z"};
       op_vec xx={};
   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {0}),spin_op(s, {1})};
	      add_state(states, v0, map_sec);
	   }
	 
   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {0})};
	      add_state(states, v0, map_sec);
	   }
      std::vector<string_pair> occ={string_pair("x","y"),string_pair("y","x"),string_pair("y","z"),string_pair("z","y"),string_pair("x","z"),string_pair("z","x")};

      for(auto a: occ){
	    
	op_vec v0={spin_op(a.first, {0}),spin_op(a.second, {1})};
	      add_state(states, v0, map_sec);
	   }
   
	  
	
  // 	auto basis=momentum_basis(L,states);
	Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
	auto basis =momentum_basis_eff(L,states,M);
    	//auto block =

  // 	for(auto a : block.total_refs_)
  // 	  {
  // 	    //	    std::cout<< a.first<< "  "<<std::endl;
  // 	  }
   	double J=1;
   	double Delta=1;
   	auto h=define_xxz1d( basis.total_refs_,basis.TI_map_, J, Delta);
	
  // // 	// for(auto a : block.G_variables)
  // // 	//   {
  // // 	//     // std::cout<< a.first<<std::endl;
  // // 	//   }	

        basis.M_->objective(ObjectiveSense::Minimize, h);
	//   // 	  block.M_->dataReport();
	//   // M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
	   basis.M_->solve();
	  
	  
	  std::cout << "Solution : " << std::endl;
	  std::cout<<M->primalObjValue()  <<std::endl;
	  
	  double sol=M->primalObjValue();
	  
	  if(std::abs(sol+0.446701)>1e-06)
	    {std::cout<<"error, not converging properly"<<std::endl;}
	    return;}
void test_multiple_blocks()
{
  int L=5;
 basis_structure states;
  std::vector<op_vec> v_block_0;
  std::vector<op_vec> v_block_1;
  std::vector<op_vec> v_block_2;
       std::vector<op_vec> v_block_3;
       states.insert({0, v_block_0});
     states.insert({1, v_block_1});
     states.insert({2, v_block_2});
     states.insert({3, v_block_3});
    std::map<std::pair<int,int>, int> map_sec
;    map_sec.insert({std::pair<int,int>(1,1), 0});
    map_sec.insert({std::pair<int,int>(1,-1), 1});
    map_sec.insert({std::pair<int,int>(-1,1), 2});
    map_sec.insert({std::pair<int,int>(-1,-1), 3});
 
       std::vector<std::string> dirs={"x", "y", "z"};
       op_vec xx={};
   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {0}),spin_op(s, {1})};
	      add_state(states, v0, map_sec);
	   }
	 
   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {0})};
	      add_state(states, v0, map_sec);
	   }
      std::vector<string_pair> occ={string_pair("x","y"),string_pair("y","x"),string_pair("y","z"),string_pair("z","y"),string_pair("x","z"),string_pair("z","x")};

      for(auto a: occ){
	    
	op_vec v0={spin_op(a.first, {0}),spin_op(a.second, {1})};
	      add_state(states, v0, map_sec);
	   }
   
	  
    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    	auto basis =momentum_basis(L,states,M);

   	double J=1;
   	double Delta=1;
   	auto h=define_xxz1d( basis.total_refs_,basis.TI_map_, J, Delta);
	
         basis.M_->objective(ObjectiveSense::Minimize, h);
	 M->dataReport();
	 M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
	   basis.M_->solve();
	  
	  
	  std::cout << "Solution : " << std::endl;
	  std::cout<<std::setprecision(9)<<M->primalObjValue()  <<std::endl;
	  
	  double sol=M->primalObjValue();
	  

	   if(std::abs(sol+0.44670126)>1e-06)
	     {std::cout<<"error, not converging properly"<<std::endl;}
	    return;}
void test_multiple_blocks_higher_order()
{
  int L=5;
 basis_structure states;
  std::vector<op_vec> v_block_0;
  std::vector<op_vec> v_block_1;
  std::vector<op_vec> v_block_2;
       std::vector<op_vec> v_block_3;
       states.insert({0, v_block_0});
     states.insert({1, v_block_1});
     states.insert({2, v_block_2});
     states.insert({3, v_block_3});
    std::map<std::pair<int,int>, int> map_sec
;    map_sec.insert({std::pair<int,int>(1,1), 0});
    map_sec.insert({std::pair<int,int>(1,-1), 1});
    map_sec.insert({std::pair<int,int>(-1,1), 2});
    map_sec.insert({std::pair<int,int>(-1,-1), 3});
 
       std::vector<std::string> dirs={"x", "y", "z"};
       op_vec xx={};
   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {0}),spin_op(s, {1})};
	      add_state(states, v0, map_sec);
	   }
      for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {0}),spin_op(s, {2})};
	      add_state(states, v0, map_sec);
	   }
	 
   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {0})};
	      add_state(states, v0, map_sec);
	   }
      std::vector<string_pair> occ={string_pair("x","y"),string_pair("y","x"),string_pair("y","z"),string_pair("z","y"),string_pair("x","z"),string_pair("z","x")};

      for(auto a: occ){
	    
	op_vec v0={spin_op(a.first, {0}),spin_op(a.second, {1})};
	      add_state(states, v0, map_sec);
	   }
   
      for(auto a: occ){
	    
	op_vec v0={spin_op(a.first, {0}),spin_op(a.second, {2})};
	      add_state(states, v0, map_sec);
	   }

      
    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    	auto basis =momentum_basis(L,states,M);

  // // 	for(auto a : block.total_refs_)
  // 	  {
  // 	    //	    std::cout<< a.first<< "  "<<std::endl;
  // 	  }
   	double J=1;
   	double Delta=1;
   	auto h=define_xxz1d( basis.total_refs_,basis.TI_map_, J, Delta);
	
  // // 	// for(auto a : block.G_variables)
  // // 	//   {
  // // 	//     // std::cout<< a.first<<std::endl;
  // // 	//   }	

         basis.M_->objective(ObjectiveSense::Minimize, h);
	//   // 	  block.M_->dataReport();
	//   // M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
	   basis.M_->solve();
	  
	  
	  std::cout << "Solution : " << std::endl;
	  std::cout<<std::setprecision(9)<<M->primalObjValue()  <<std::endl;
	  std::cout<< "true ground state value "<<-0.37360679774997907<<std::endl;
	  
	  double sol=M->primalObjValue();
	  

	    if(std::abs(sol+0.389706358)>1e-06)
	      {std::cout<<"error, not converging properly"<<std::endl;}
	    return;}
// two dimensional chain
void test_multiple_blocks_2d()
{
  int Lx=6;
  int Ly=6;
 basis_structure states;
  std::vector<op_vec> v_block_0;
  std::vector<op_vec> v_block_1;
  std::vector<op_vec> v_block_2;
  std::vector<op_vec> v_block_3;
  states.insert({0, v_block_0});
  states.insert({1, v_block_1});
  states.insert({2, v_block_2});
  states.insert({3, v_block_3});
  std::map<std::pair<int,int>, int> map_sec;
      map_sec.insert({std::pair<int,int>(1,1), 0});
  map_sec.insert({std::pair<int,int>(1,-1), 1});
  map_sec.insert({std::pair<int,int>(-1,1), 2});
    map_sec.insert({std::pair<int,int>(-1,-1), 3});
 
       std::vector<std::string> dirs={"x", "y", "z"};
       for(int i=0; i<Ly;i++)
	 {
	   auto mn=std::min(i, (i+1)%Ly);
	   auto mx=std::max(i, (i+1)%Ly);
	   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {i,0}, Lx),spin_op(s, {i,1}, Lx)};
	     add_state(states, v0, map_sec);
	     op_vec v1={spin_op(s, {mn,0}, Lx),spin_op(s, {mx,0}, Lx)};
	     add_state(states, v1, map_sec);

	    }
	 }

for(int i=0; i<Ly;i++)
	 {

	   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {i,0}, Lx)};
	     add_state(states, v0, map_sec);


	    }
	 }


       
      std::vector<string_pair> occ={string_pair("x","y"),string_pair("y","x"),string_pair("y","z"),string_pair("z","y"),string_pair("x","z"),string_pair("z","x")};   
	        for(int i=0; i<Ly;i++)
	 {
	   auto mn=std::min(i, (i+1)%Ly);
	   auto mx=std::max(i, (i+1)%Ly);
	   for(auto a: occ){
	     op_vec v0={spin_op(a.first, {i,0}, Lx),spin_op(a.second, {i,1}, Lx)};
	     add_state(states, v0, map_sec);
	     op_vec v1={spin_op(a.first, {mn,0}, Lx),spin_op(a.second, {mx,0}, Lx)};
	     add_state(states, v1, map_sec);
	   }
	 }
    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    auto basis =momentum_basis_eff(Lx,states,M, "xyz");
    // for(auto a: basis.TI_map_)
    //   {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
    double J=1;
    double Delta=1;
    auto h=define_xxz2d( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
    
    basis.M_->objective(ObjectiveSense::Minimize, h);
	 	  M->dataReport();
	 M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
    basis.M_->solve();
	  
	  
    std::cout << "Solution : " << std::endl;
    std::cout<<std::setprecision(9)<<M->primalObjValue()/Ly  <<std::endl;
	  
     double sol=M->primalObjValue()/Ly;
	  

	  	   if(std::abs(sol+0.721905655)>1e-06)
	   {std::cout<<"error, not converging properly"<<std::endl;}
	    return;
}

void test_multiple_blocks_higher_order_2d()
{
  std::cout<< "WARNING! Takes a lot of memory"<<std::endl;
  int Lx=4;
  int Ly=4;
 basis_structure states;
  std::vector<op_vec> v_block_0;
  std::vector<op_vec> v_block_1;
  std::vector<op_vec> v_block_2;
  std::vector<op_vec> v_block_3;
  states.insert({0, v_block_0});
  states.insert({1, v_block_1});
  states.insert({2, v_block_2});
  states.insert({3, v_block_3});
  std::map<std::pair<int,int>, int> map_sec;
      map_sec.insert({std::pair<int,int>(1,1), 0});
  map_sec.insert({std::pair<int,int>(1,-1), 1});
  map_sec.insert({std::pair<int,int>(-1,1), 2});
    map_sec.insert({std::pair<int,int>(-1,-1), 3});
 
       std::vector<std::string> dirs={"x", "y", "z"};
       for(int i=0; i<Ly;i++)
	 {
	   auto mn=std::min(i, (i+1)%Ly);
	   auto mx=std::max(i, (i+1)%Ly);
	   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {i,0}, Lx),spin_op(s, {i,1}, Lx)};
	     add_state(states, v0, map_sec);
	     op_vec v1={spin_op(s, {mn,0}, Lx),spin_op(s, {mx,0}, Lx)};
	     add_state(states, v1, map_sec);

	    }
	 }
       for(int i=0; i<Ly;i++)
	 {
	   auto mn=std::min(i, (i+2)%Ly);
	   auto mx=std::max(i, (i+2)%Ly);
	   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {i,0}, Lx),spin_op(s, {i,2}, Lx)};
	     add_state(states, v0, map_sec);
	     op_vec v1={spin_op(s, {mn,0}, Lx),spin_op(s, {mx,0}, Lx)};
	     add_state(states, v1, map_sec);

	    }
	 }

       
for(int i=0; i<Ly;i++)
	 {

	   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {i,0}, Lx)};
	     add_state(states, v0, map_sec);


	    }
	 }


       
      std::vector<string_pair> occ={string_pair("x","y"),string_pair("y","x"),string_pair("y","z"),string_pair("z","y"),string_pair("x","z"),string_pair("z","x")};   
	        for(int i=0; i<Ly;i++)
	 {
	   auto mn=std::min(i, (i+1)%Ly);
	   auto mx=std::max(i, (i+1)%Ly);
	   for(auto a: occ){
	     op_vec v0={spin_op(a.first, {i,0}, Lx),spin_op(a.second, {i,1}, Lx)};
	     add_state(states, v0, map_sec);
	     op_vec v1={spin_op(a.first, {mn,0}, Lx),spin_op(a.second, {mx,0}, Lx)};
	     add_state(states, v1, map_sec);
	   }
	 }
			        for(int i=0; i<Ly;i++)
	 {
	   auto mn=std::min(i, (i+2)%Ly);
	   auto mx=std::max(i, (i+2)%Ly);
	   for(auto a: occ){
	     op_vec v0={spin_op(a.first, {i,0}, Lx),spin_op(a.second, {i,2}, Lx)};
	     add_state(states, v0, map_sec);
	     op_vec v1={spin_op(a.first, {mn,0}, Lx),spin_op(a.second, {mx,0}, Lx)};
	     add_state(states, v1, map_sec);
	   }
	 }
    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    auto basis =momentum_basis_eff(Lx,states,M,"xyz");
    // for(auto a: basis.TI_map_)
    //   {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
    double J=1;
    double Delta=1;
    auto h=define_xxz2d( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
    
    basis.M_->objective(ObjectiveSense::Minimize, h);
	//   // 	  block.M_->dataReport();
	//   // M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
    basis.M_->solve();
	  
	  
    std::cout << "Solution : " << std::endl;
    std::cout<<std::setprecision(9)<<M->primalObjValue()/Ly  <<std::endl;
	  
    // double sol=M->primalObjValue(); -0.720847431
	  

	  //	   if(std::abs(sol+0.44670126)>1e-06)
	  // {std::cout<<"error, not converging properly"<<std::endl;}
	    return;
}


void test_J1J2_2d()
{

  int Lx=6;
  int Ly=6;
 basis_structure states;
  std::vector<op_vec> v_block_0;
  std::vector<op_vec> v_block_1;
  std::vector<op_vec> v_block_2;
  std::vector<op_vec> v_block_3;
  states.insert({0, v_block_0});
  states.insert({1, v_block_1});
  states.insert({2, v_block_2});
  states.insert({3, v_block_3});
std::map<std::pair<int,int>, int> map_sec;

 map_sec.insert({std::pair<int,int>(1,1), 0});
  map_sec.insert({std::pair<int,int>(1,-1), 1});
  map_sec.insert({std::pair<int,int>(-1,1), 2});
    map_sec.insert({std::pair<int,int>(-1,-1), 3});
 
       std::vector<std::string> dirs={"x", "y", "z"};
       for(int i=0; i<Ly;i++)
	 {
	   auto mn=std::min(i, (i+1)%Ly);
	   auto mx=std::max(i, (i+1)%Ly);
	   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {i,0}, Lx),spin_op(s, {i,1}, Lx)};
	     add_state(states, v0, map_sec);
	     op_vec v1={spin_op(s, {mn,0}, Lx),spin_op(s, {mx,0}, Lx)};
	     //add_state(states, v1, map_sec);

	    }
	 }
              for(int i=0; i<Ly;i++)
	 {
	   auto mn=std::min(i, (i+2)%Ly);
	   auto mx=std::max(i, (i+2)%Ly);
	   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {i,0}, Lx),spin_op(s, {i,2}, Lx)};
	     add_state(states, v0, map_sec);
	     op_vec v1={spin_op(s, {mn,0}, Lx),spin_op(s, {mx,0}, Lx)};
	     add_state(states, v1, map_sec);

	    }
	 }
for(int i=0; i<Ly;i++)
	 {

	   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {i,0}, Lx)};
	     add_state(states, v0, map_sec);


	    }
	 }

    std::vector<string_pair> occ={string_pair("x","y"),string_pair("y","x"),string_pair("y","z"),string_pair("z","y"),string_pair("x","z"),string_pair("z","x")};   
	 //        for(int i=0; i<Ly;i++)
	 // {
	 //   auto mn=std::min(i, (i+1)%Ly);
	 //   auto mx=std::max(i, (i+1)%Ly);
	 //   for(auto a: occ){
	 //     op_vec v0={spin_op(a.first, {i,0}, Lx),spin_op(a.second, {i,1}, Lx)};
	 //     add_state(states, v0, map_sec);
	 //     op_vec v1={spin_op(a.first, {mn,0}, Lx),spin_op(a.second, {mx,0}, Lx)};
	 //     add_state(states, v1, map_sec);
	 //   }
	 // }
    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    auto basis =momentum_basis_eff(Lx,states,M,"xyz");
    // for(auto a: basis.TI_map_)
    //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
    double J1=1;
    double J2=1.0;
    
     auto h=define_J1J2_2d( basis.total_refs_,basis.TI_map_, J1, J2, Ly, Lx);
    
    basis.M_->objective(ObjectiveSense::Minimize, h);
		  basis.M_->dataReport();
	  M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
    basis.M_->solve();
	  
	  
    std::cout << "Solution : " << std::endl;
    std::cout<<std::setprecision(9)<<M->primalObjValue()  /Ly<<std::endl;
	  
    // double sol=M->primalObjValue(); -0.720847431
	  

	  //	   if(std::abs(sol+0.44670126)>1e-06)
	  // {std::cout<<"error, not converging properly"<<std::endl;}
	    return;
}
