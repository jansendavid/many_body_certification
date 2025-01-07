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
#include"sos/complex_momentum_dual.hpp"
#include"spin_hamiltonians_TIsym.hpp"
#include"functions.hpp"
#include"reduced_dms.hpp"
using namespace mosek::fusion;
using namespace monty;


void test_multiple_blocks_higher_order_2d_rdm()
{
  std::cout<< "WARNING! Takes a lot of memory"<<std::endl;
  int Lx=4;
  int Ly=4;
 basis_structure states=get_basis_2d(Lx, 1, 0, true);

 
for(auto a: states)
{

	std::cout<<"sec "<< a.first<< " and size "<< a.second.size()<< " and "<<a.second.size()/(Lx*Lx)<<std::endl;
	for(auto l: a.second)
	{
	//	std::cout<< print_op(l)<<std::endl;
	}
}

    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    auto basis =momentum_symmetry_solver_dual(Lx,states,M,"xyz");
    //for(auto a: basis.TI_map_)
     //{std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
     double J=1;
     double Delta=1.;
    
     auto b=define_xxz2d_sos( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
 
        basis.set_b(b);
        basis.fix_constrains();
       auto h=basis.get_costfunction();
	//   std::cout<<C->toString()<<std::endl;
	  
	  {std::pair<int,int> a(0,0);
	  std::pair<int,int> b(0,1);
	//generate_rmds_primal({a, b},basis.total_refs_,basis.TI_map_ , basis.y_, Lx, M);
	  }
	// // //     {std::pair<int,int> a(0,0);
	// // //   std::pair<int,intget_basis_2d(Lx, 3, -3, true);> b(0,1);
	// // //    std::pair<int,int> c(0,2);
	// // // //generate_rmds_primal({a, b, c},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
	// // //   }
    basis.M_->objective(ObjectiveSense::Minimize, h);
		  basis.M_->dataReport();
	  M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
    basis.M_->solve();
	  
	  
    std::cout << "Solution : " << std::endl;
    std::cout<<std::setprecision(9)<<M->primalObjValue()<<std::endl;
	  
   //double sol=M->primalObjValue(); //-0.720847431
	  

	  //	   if(std::abs(sol+0.44670126)>1e-06)
	  // {std::cout<<"error, not converging properly"<<std::endl;}
	    return;
}

void test_multiple_blocks_higher_order_2d_rdm_sos()
{
  std::cout<< "WARNING! Takes a lot of memory"<<std::endl;
  int Lx=4;
  int Ly=4;

    basis_structure states=get_basis_2d(Lx, 1, 0, true);


for(auto a: states)
{

	std::cout<<"sec "<< a.first<< " and size "<< a.second.size()<< " and "<<a.second.size()/(Lx*Lx)<<std::endl;
	  //for(auto n:a.second)
     //{std::cout<<print_op(n)<<std::endl;}
}

    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    auto basis =momentum_symmetry_solver_sos(Lx,states,M,"xyz");
    // for(auto a: basis.TI_map_)
    //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
  for(auto k: basis.total_refs_)
  {

//	std::cout<<k.first<<std::endl; 
  }
     double J=1;
     double Delta=1.;
    
     auto b=define_xxz2d_sos( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
 std::cout<< "x"<<std::endl;
        basis.set_b(b);
        
        basis.fix_constrains();
       auto h=basis.get_costfunction();
	// //   std::cout<<C->toString()<<std::endl;
	  
	  {std::pair<int,int> a(0,0);
	  std::pair<int,int> b(1,0);
	//generate_rmds_primal({a, b},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
	  }
	    {std::pair<int,int> a(0,0);
	  std::pair<int,int> b(0,1);
	   std::pair<int,int> c(0,2);
	//generate_rmds_primal({a, b, c},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
	  }

  std::cout<< "starting solving SDP"<<std::endl;
    basis.M_->objective(ObjectiveSense::Maximize, h);
		  basis.M_->dataReport();
	  M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
    basis.M_->solve();
	  auto cons=M->getConstraint(0)->dual();
    std::cout<<cons<<std::endl;
	  
    std::cout << "Solution : " << std::endl;
    std::cout<<std::setprecision(9)<<M->primalObjValue()<<std::endl;
	  
// //    //double sol=M->primalObjValue(); //-0.720847431
	  

// // 	  //	   if(std::abs(sol+0.44670126)>1e-06)
// // 	  // {std::cout<<"error, not converging properly"<<std::endl;}
	    return;
}
void test_x()
{
   
 int Lx=4;
  int Ly=4;

    auto states=get_basis_2d_new(Lx, 2);
	std::vector<int> sizes(4,0);
	for(auto sign: states)
{
// 	std::cout<< "sign "<<sign.first<<std::endl;
// 	for(auto b: sign.second)
// 	{
// 		for(auto g: b)
// 		{
// 			std::cout<<print_op(g)<<std::endl;
// 		}
 		sizes[sign.first]+=sign.second.size();
 	}
// }
// std::cout<< "sizes "<<std::endl;
for(auto b: sizes)
{std::cout<<b<<std::endl;}

	    return;
}

void test_y()
{
  int L=4;
 op_vec vec={spin_op("x", {2,2}, {L,L}),spin_op("x", {0,0}, {L,L})};

std::cout<< print_op(vec)<<std::endl;
 auto [fac, nf] =get_normal_form(vec);
 std::cout<< fac << " and "<< print_op(nf)<<std::endl;

	    return;
}

void test_d8_symm()
{
  int L=6;
//  op_vec vec={spin_op("x", {0,0}, {L,L}),spin_op("x", {1,1}, {L,L}),spin_op("y", {2,1}, {L,L})};
// std::cout<< print_op(vec)<<std::endl;
// auto all_d8=generate_all_d8(vec,  L);
// std::cout<< "start "<<std::endl;
// for(auto b: all_d8)
// {
// std::cout<<print_op(b)<<std::endl;

// }
//s_[x,(1,0)]s_[y,(2,0)]s_[x,(4,3)]s_[y,(5,3)]
// s_[x,(1,0)]s_[y,(2,0)]s_[x,(4,3)]s_[y,(5,3)]
// s_[x,(0,5)]s_[y,(0,4)]s_[x,(3,2)]s_[y,(3,1)]
// s_[x,(5,0)]s_[y,(4,0)]s_[x,(2,3)]s_[y,(1,3)]
// s_[x,(0,1)]s_[y,(0,2)]s_[x,(3,4)]s_[y,(3,5)]
// s_[x,(5,0)]s_[y,(4,0)]s_[x,(2,3)]s_[y,(1,3)]
// s_[x,(0,1)]s_[y,(0,2)]s_[x,(3,4)]s_[y,(3,5)]
// s_[x,(1,0)]s_[y,(2,0)]s_[x,(4,3)]s_[y,(5,3)]
// s_[x,(0,5)]s_[y,(0,4)]s_[x,(3,2)]s_[y,(3,1)]
op_vec vec={spin_op("x", {1,0}, {L,L}),spin_op("y", {2,0}, {L,L}),spin_op("x", {4,3}, {L,L}),spin_op("y", {5,3}, {L,L})};
std::cout<<print_op(vec)<<std::endl;
auto all_d8=generate_all_d8(vec,  L);
for(auto b: all_d8)
{
std::cout<<print_op(b)<<std::endl;

}
	    return;
}

// void test_J1J2_2d()
// {

//   int Lx=6;
//   int Ly=6;
//  basis_structure states;
//   std::vector<op_vec> v_block_0;
//   std::vector<op_vec> v_block_1;
//   std::vector<op_vec> v_block_2;
//   std::vector<op_vec> v_block_3;
//   states.insert({0, v_block_0});
//   states.insert({1, v_block_1});
//   states.insert({2, v_block_2});
//   states.insert({3, v_block_3});
// std::map<std::pair<int,int>, int> map_sec;

//  map_sec.insert({std::pair<int,int>(1,1), 0});
//   map_sec.insert({std::pair<int,int>(1,-1), 1});
//   map_sec.insert({std::pair<int,int>(-1,1), 2});
//     map_sec.insert({std::pair<int,int>(-1,-1), 3});
 
//        std::vector<std::string> dirs={"x", "y", "z"};
//        for(int i=0; i<Ly;i++)
// 	 {
// 	   auto mn=std::min(i, (i+1)%Ly);
// 	   auto mx=std::max(i, (i+1)%Ly);
// 	   for(auto s: dirs){
	    
// 	     op_vec v0={spin_op(s, {i,0}, Lx),spin_op(s, {i,1}, Lx)};
// 	     add_state(states, v0, map_sec);
// 	     op_vec v1={spin_op(s, {mn,0}, Lx),spin_op(s, {mx,0}, Lx)};
// 	     //add_state(states, v1, map_sec);

// 	    }
// 	 }
//               for(int i=0; i<Ly;i++)
// 	 {
// 	   auto mn=std::min(i, (i+2)%Ly);
// 	   auto mx=std::max(i, (i+2)%Ly);
// 	   for(auto s: dirs){
	    
// 	     op_vec v0={spin_op(s, {i,0}, Lx),spin_op(s, {i,2}, Lx)};
// 	     add_state(states, v0, map_sec);
// 	     op_vec v1={spin_op(s, {mn,0}, Lx),spin_op(s, {mx,0}, Lx)};
// 	     add_state(states, v1, map_sec);

// 	    }
// 	 }
// for(int i=0; i<Ly;i++)
// 	 {

// 	   for(auto s: dirs){
	    
// 	     op_vec v0={spin_op(s, {i,0}, Lx)};
// 	     add_state(states, v0, map_sec);


// 	    }
// 	 }

//     std::vector<string_pair> occ={string_pair("x","y"),string_pair("y","x"),string_pair("y","z"),string_pair("z","y"),string_pair("x","z"),string_pair("z","x")};   
// 	        for(int i=0; i<Ly;i++)
// 	 {
// 	   auto mn=std::min(i, (i+1)%Ly);
// 	   auto mx=std::max(i, (i+1)%Ly);
// 	   for(auto a: occ){
// 	     op_vec v0={spin_op(a.first, {i,0}, Lx),spin_op(a.second, {i,1}, Lx)};
// 	     add_state(states, v0, map_sec);
// 	     op_vec v1={spin_op(a.first, {mn,0}, Lx),spin_op(a.second, {mx,0}, Lx)};
// 	     add_state(states, v1, map_sec);
// 	   }
// 	 }
//     Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
//     auto basis =momentum_basis_eff(Lx,states,M,"xyz");
//     // for(auto a: basis.TI_map_)
//     //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
//     double J1=1;
//     double J2=0.2;
    
//      auto C=define_J1J2_2d_dual( basis.total_refs_,basis.TI_map_, J1, J2, Ly, Lx);
 
//        basis.set_C(C);
//       auto h=basis.get_costfunction();


//     basis.M_->objective(ObjectiveSense::Minimize, h);
// 		  basis.M_->dataReport();
// 	  M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
//     basis.M_->solve();
	  
	  
//     std::cout << "Solution : " << std::endl;
//     std::cout<<std::setprecision(9)<<M->primalObjValue()  /Ly<<std::endl;
	  
//     // double sol=M->primalObjValue(); -0.720847431
	  

// 	  //	   if(std::abs(sol+0.44670126)>1e-06)
// 	  // {std::cout<<"error, not converging properly"<<std::endl;}
// 	    return;
// }
// void test_load_basis()
// {
//   int Lx=4;
//   auto states=load_basis_from_file("test_basis.txt", Lx);
//   for(auto s: states)
//     {std::cout<<print_op(s)<<std::endl;}
// }
