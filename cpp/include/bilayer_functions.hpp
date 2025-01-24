#pragma once
#include "fusion.h"
#include"spins.hpp"
#include<unordered_map>
#include <memory>
#include"symmetries.hpp"
using namespace mosek::fusion;
using namespace monty;

void get_order_one_monomials_bilayer(basis_structure& states, std::map<std::pair<int,int>, int>& map_sec, int layers, int Ly,int Lx,  bool use_symm)
{

assert(Lx==Ly);
   std::vector<std::string> dirs={"x","y","z"};
       
 
for(int i=0; i<layers; i++)
{
	   for(auto s: dirs){
	    
 	     op_vec v0={spin_op(s, {i,0,0}, {layers, Ly,Lx})};

         auto sign=get_sec( v0);

	         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, Lx);
        }
        else{
          add_state(states, v0, map_sec);
        }
 


     }
 	    }
}

void get_order_two_monomials_bilayer(basis_structure& states, std::map<std::pair<int,int>, int>& map_sec,int layers, int Ly,int Lx,  int r, int start,bool use_symm)
{
std::vector<std::string> dirs={"x","y","z"};
assert(Lx==Ly);
   //    dirs={"x", "y", "z"};
    //  std::vector<std::pair<int, int>> rvals={{1,0}, {0,1}, {1,1}, {1,3}, {2,0}, {0,2}, {2,1}, {1,2}, {2,2}, {2,3}};
    //  for(auto s1: dirs){
	     
    //    for(auto s2: dirs)
    //    {
    //    //if(s1!=s2)
    //    {
    //     for(auto b: rvals)
    //     {
    //               op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {b.first, b.second}, {L,L})};

	
    //     auto [fac, vec] =get_normal_form(v0);
    //     add_state(states, vec, map_sec);
    //     }
    //    }}
    //  }
//       int SS=0;

for(int layer=0; layer<layers; layer++)
{
      for(int i=start; i<=r; i++)
      {
     
          for(int j=start; j<=r; j++)
      {
        for(auto s1: dirs){
	     
       for(auto s2: dirs)
       {
   
        if(i!=0 or j!=0)
        {
          int ind1=(Ly+i)%Ly;
          int ind2=(Lx+j)%Lx;
        op_vec v0={spin_op(s1, {0,0,0}, {layers, Ly,Lx}),spin_op(s2, {layer, ind1, ind2}, {layers, Ly,Lx})};

	
        auto [fac, vec] =get_normal_form(v0);
   
   
    
         if(use_symm)
        {
add_state_with_symmetries(states, vec, map_sec, Lx);
        }
        else{
          add_state(states, vec, map_sec);
        }
        
      

	     
        }
		 }
		   
	    }
     }
      }
}
}
// void get_order_three_monomials(basis_structure& states, std::map<std::pair<int,int>, int>& map_sec, int L, bool use_symm)
// {

// std::vector<std::string> dirs={"x","y","z"};
//        	   for(auto s1: dirs){
// 	      for(auto s2: dirs){
//            for(auto s3: dirs){
	     
//        {
//         op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {0,1},{L,L}),spin_op(s3, {1,1}, {L,L})};
//         auto [fac, vec] =get_normal_form(v0);
	       
//        {
      
//          if(use_symm)
//         {
// add_state_with_symmetries(states, v0, map_sec, L);
//         }
//         else{
//           add_state(states, v0, map_sec);
//         }
// 		 }
//         }
//       {
//         op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {0,1}, {L,L}),spin_op(s3, {L-1,1}, {L,L})};
// 	     auto [fac, vec] =get_normal_form(v0);
	    
//        {
     
//          if(use_symm)
//         {
// add_state_with_symmetries(states, v0, map_sec, L);
//         }
//         else{
//           add_state(states, v0, map_sec);
//         }
//        }
// 		 }
//      {
//           op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {1,0}, {L,L}),spin_op(s3, {1,1}, {L,L})};
// 	     auto [fac, vec] =get_normal_form(v0);
	     
//        {
      
//          if(use_symm)
//         {
// add_state_with_symmetries(states, v0, map_sec, L);
//         }
//         else{
//           add_state(states, v0, map_sec);
//         }
//        }
// 		 }
//         {
//           op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {L-1,0}, {L,L}),spin_op(s3, {L-1,1}, {L,L})};
// 	     auto [fac, vec] =get_normal_form(v0);
	       
//        {
       
//          if(use_symm)
//         {
// add_state_with_symmetries(states, v0, map_sec, L);
//         }
//         else{
//           add_state(states, v0, map_sec);
//         }
//         }
// 		 }
//       {
//           op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {1,0}, {L,L}),spin_op(s3, {2,0}, {L,L})};
// 	     auto [fac, vec] =get_normal_form(v0);
	       
//        {
      
//          if(use_symm)
//         {
// add_state_with_symmetries(states, v0, map_sec, L);
//         }
//         else{
//           add_state(states, v0, map_sec);
//         }
//         }
// 		 }
//         {
//           op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {0,1}, {L,L}),spin_op(s3, {0,2}, {L,L})};
// 	     auto [fac, vec] =get_normal_form(v0);
	      
//        {
      
//          if(use_symm)
//         {
// add_state_with_symmetries(states, v0, map_sec, L);
//         }
//         else{
//           add_state(states, v0, map_sec);
//         }
//         }
// 		 }
		     
     
// 	    }
//         }
//      }
// }
  void get_order_four_monomials_bilayer(basis_structure& states, std::map<std::pair<int,int>, int>& map_sec, int layers, int Ly,int Lx,  bool use_symm)
{
std::vector<std::string> dirs={"x","y","z"};
std::vector<int> offsets_={layers, Ly, Lx};
  for(auto s1: dirs){
	      for(auto s2: dirs){
           for(auto s3: dirs){
	     
       {
            for(auto s4: dirs){
	     
       {
        op_vec v0={spin_op(s1, {0,0,0}, offsets_),spin_op(s2, {1,1,0}, offsets_),spin_op(s3, {0,0,1}, offsets_),spin_op(s4, {1,1,1}, offsets_)};
        auto [fac, vec] =get_normal_form(v0);
	 
         if(use_symm)
        {
add_state_with_symmetries(states, vec, map_sec, Lx);
        }
        else{
          add_state(states, vec, map_sec);
        }
		 }
          {
        op_vec v0={spin_op(s1, {0,0,0}, offsets_),spin_op(s2, {1,1,0}, offsets_),spin_op(s3, {0,2,0},offsets_),spin_op(s4, {1,3,0}, {offsets_})};
        auto [fac, vec] =get_normal_form(v0);
	 
         if(use_symm)
        {
add_state_with_symmetries(states, vec, map_sec, Lx);
        }
        else{
          add_state(states, vec, map_sec);
        }
		 }
            } }}}}
}
rdms_struct  get_rdms_bilayer(int Lx, int dim)
{
  rdms_struct data;
     
        rdm_operator newstate({{0,0,0}, {0,0,1}});
data.add_operator(newstate);
    rdm_operator newstate_2({{1,0,0}, {1,0,1}});
data.add_operator(newstate_2);
    rdm_operator newstate_3({{0,0,0}, {1,0,1}});
data.add_operator(newstate_3);
  rdm_operator newstate_4({{0,0,0}, {0,1,1}});
data.add_operator(newstate_4);
 rdm_operator newstate_5({{0,0,0}, {1,1,1}});
data.add_operator(newstate_5);
 rdm_operator newstate_6({{1,0,0}, {1,1,1}});
data.add_operator(newstate_6);

if(dim>=4)
{
     rdm_operator newstate({{0,0,0}, {0,0,1}, {0,1,0}});
 data.add_operator(newstate);

   rdm_operator newstate_1({{0,0,0}, {0,1,0}, {0,0,1}});
 data.add_operator(newstate_1);
 rdm_operator newstate_2({{0,0,0}, {0,1,0}, {0,1,1}});
 data.add_operator(newstate_2);

  rdm_operator newstate_3({{0,0,0}, {1,1,0}, {1,1,1}});
 data.add_operator(newstate_3);

  rdm_operator newstate_4({{0,0,0}, {1,1,0}, {0,1,1}});
 data.add_operator(newstate_4);
  rdm_operator newstate_5({{1,0,0}, {0,1,0}, {0,1,1}});
 data.add_operator(newstate_5);


    rdm_operator newstate_x({{0,0,0}, {0,0,1}, {0,1,0},{1,1,0}});
 data.add_operator(newstate_x);

   rdm_operator newstate_x1({{0,0,0}, {0,1,0}, {0,0,1},{1,0,1}});
 data.add_operator(newstate_x1);
 rdm_operator newstate_x2({{0,0,0}, {0,1,0}, {0,1,1},{1,1,1}});
 data.add_operator(newstate_x2);

  rdm_operator newstate_x3({{0,0,0}, {1,1,0}, {1,1,1}, {1,2,1}});
 data.add_operator(newstate_x3);

  rdm_operator newstate_x4({{0,0,0}, {1,1,0}, {0,1,1},{0,1,2}});
 data.add_operator(newstate_x4);
  rdm_operator newstate_x5({{1,0,0}, {0,1,0}, {0,1,1},{1,2,2}});
 data.add_operator(newstate_x5);


 }
// if(dim>=6)
// {
//     rdm_operator newstate({{0,0}, {0,1}, {0,2},{0,3},{0,4}});
// data.add_operator(newstate);
//  rdm_operator newstate_1({{0,0}, {0,1}, {0,2},{0,3},{0,4},{0,5}});
// data.add_operator(newstate_1);
// }
  
return data;

}