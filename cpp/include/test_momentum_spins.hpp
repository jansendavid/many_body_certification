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
#include"complex_momentum_symm.hpp"
#include"spin_hamiltonians_TIsym.hpp"
using namespace mosek::fusion;
using namespace monty;
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


void test_single_block()
{
  int L=5;
      basis_structure states;
      op_vec v1={spin_op("x", {0})};
      op_vec v2={spin_op("y", {0})};
     op_vec v3={spin_op("z", {0})};
    op_vec v4={spin_op("x", {0}),spin_op("x", {1})};
     op_vec v5={spin_op("y", {0}),spin_op("y", {1})};
     op_vec v6={spin_op("z", {0}),spin_op("z", {1})};
     
     op_vec v7={spin_op("x", {0}),spin_op("y", {1})};
     op_vec v8={spin_op("y", {0}),spin_op("x", {1})};
     op_vec v9={spin_op("y", {0}),spin_op("z", {1})};
     op_vec v10={spin_op("z", {0}),spin_op("y", {1})};
     op_vec v11={spin_op("x", {0}),spin_op("z", {1})};
     op_vec v12={spin_op("z", {0}),spin_op("x", {1})};
  //   op_vec v5={electron_op("up", {0}, false),electron_op("up", {1}, true)};
  //   //    op_vec v5={electron_op("up", {1}, true),electron_op("up", {0}, false)};

    
    std::vector<op_vec> v_tot;
    v_tot.push_back(v1);
    v_tot.push_back(v2);
    v_tot.push_back(v3);
    v_tot.push_back(v4);
    v_tot.push_back(v5);
    v_tot.push_back(v6);
    v_tot.push_back(v7);
    v_tot.push_back(v8);
    v_tot.push_back(v9);
    v_tot.push_back(v10);
    v_tot.push_back(v11);
    v_tot.push_back(v12);
    states.insert({0, v_tot});
	
  // 	auto basis=momentum_basis(L,states);
	Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
	auto basis =momentum_basis(L,states,M);
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
	  
	  if(std::abs(sol+0.464015)>1e-06)
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

     
    op_vec v0={spin_op("x", {0}),spin_op("x", {1})};
     op_vec v1={spin_op("y", {0}),spin_op("y", {1})};
     op_vec v2={spin_op("z", {0}),spin_op("z", {1})};
     v_block_0.push_back(v0);
     v_block_0.push_back(v1);
     v_block_0.push_back(v2);
     op_vec v3={spin_op("z", {0})};
      op_vec v4={spin_op("x", {0})};
      op_vec v5={spin_op("y", {0})};
      v_block_1.push_back(v3);
       v_block_2.push_back(v4);
       v_block_3.push_back(v5);
       op_vec v6={spin_op("x", {0}),spin_op("y", {1})};
      op_vec v7={spin_op("y", {0}),spin_op("x", {1})};
       v_block_1.push_back(v6);
      v_block_1.push_back(v7);
      op_vec v8={spin_op("y", {0}),spin_op("z", {1})};
      op_vec v9={spin_op("z", {0}),spin_op("y", {1})};
      v_block_2.push_back(v8);
      v_block_2.push_back(v9);
      
      op_vec v10={spin_op("x", {0}),spin_op("z", {1})};
      op_vec v11={spin_op("z", {0}),spin_op("x", {1})};
      v_block_3.push_back(v10);
      v_block_3.push_back(v11);
    
    
    states.insert({0, v_block_0});
    states.insert({1, v_block_1});
    states.insert({2, v_block_2});
    states.insert({3, v_block_3});
    
    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    	auto basis =momentum_basis(L,states,M);

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
	  std::cout<<std::setprecision(9)<<M->primalObjValue()  <<std::endl;
	  
	  double sol=M->primalObjValue();
	  

	  // if(std::abs(sol+0.464015)>1e-06)
	  //   {std::cout<<"error, not converging properly"<<std::endl;}
	    return;}
