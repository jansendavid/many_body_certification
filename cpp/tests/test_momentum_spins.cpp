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
Expression::t define_ham( std::map<std::string, mom_ref> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, double J, double Delta)
{
  std::vector<Expression::t> expressions={};
     op_vec v0={spin_op("x", {0}),spin_op("x", {1})};
     auto [fac0, nf0] =get_normal_form(v0);
  auto [ key0,coeff0_map]=map.at(print_op(nf0));  
  auto el0=refs.at(key0);

  
       op_vec v1={spin_op("y", {0}),spin_op("y", {1})};
     auto [fac1, nf1] =get_normal_form(v1);
  auto [ key1,coeff1_map]=map.at(print_op(nf1));
  auto el1=refs.at(key1);
    std::cout<<"key 1"<< fac1 << " "<<coeff1_map<<std::endl;

       op_vec v2={spin_op("z", {0}),spin_op("z", {1})};
     auto [fac2, nf2] =get_normal_form(v2);
  auto [ key2,coeff2_map]=map.at(print_op(nf2));
    std::cout<<"key 2"<< fac2 << " "<<coeff2_map<<std::endl;
   auto el2=refs.at(key2);
   
   Expression::t ham=Expr::add(Expr::mul(Delta*fac2.real()*coeff2_map.real()/4., el2.var_),Expr::add(Expr::mul(J*fac1.real()*coeff1_map.real()/4., el1.var_),Expr::mul(J*fac0.real()*coeff0_map.real()/4., el0.var_)));
  

  return ham;
}

// Expression::t particle_number( std::unordered_map<std::string, mom_ref> refs)
// {
//   std::vector<Expression::t> expressions={};
  
//   op_vec v4={electron_op("up", {0}, true),electron_op("up", {0}, false)};
//   auto coeff4= bubbleSort(v4, v4.size());
  
//   auto el4=refs.at(print_op(v4));
  
//   Expression::t ham=Expr::mul(1, el4.var_);

//   return ham;
// }


int main()
{
  //test_translation();
    int L=5;
      basis_structure states;
      op_vec v1={spin_op("x", {0})};
      op_vec v2={spin_op("y", {0})};
     op_vec v3={spin_op("z", {0})};
    op_vec v4={spin_op("x", {0}),spin_op("x", {1})};
     op_vec v5={spin_op("y", {0}),spin_op("y", {1})};
     op_vec v6={spin_op("z", {0}),spin_op("z", {1})};
     
     // op_vec v7={spin_op("x", {0}),spin_op("y", {1})};
     // op_vec v8={spin_op("y", {0}),spin_op("z", {1})};
     // op_vec v9={spin_op("x", {0}),spin_op("z", {1})};
  //   op_vec v5={electron_op("up", {0}, false),electron_op("up", {1}, true)};
  //   //    op_vec v5={electron_op("up", {1}, true),electron_op("up", {0}, false)};

    
    std::vector<op_vec> v_tot;
    v_tot.push_back(v1);
    v_tot.push_back(v2);
    v_tot.push_back(v3);
     v_tot.push_back(v4);
     v_tot.push_back(v5);
     v_tot.push_back(v6);
    // v_tot.push_back(v7);
    // v_tot.push_back(v8);
    // v_tot.push_back(v9);
    states.insert({0, v_tot});
	
  // // 	auto basis=momentum_basis(L,states);
   	Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    	auto block =momentum_zero_block(L,v_tot,M);

	for(auto a : block.total_refs_)
	  {
	    //	    std::cout<< a.first<< "  "<<std::endl;
	  }
   	double J=1;
   	double Delta=1;
   	auto h=define_ham( block.total_refs_,block.TI_map, J, Delta);
	
  // 	// for(auto a : block.G_variables)
  // 	//   {
  // 	//     // std::cout<< a.first<<std::endl;
  // 	//   }	

        block.M_->objective(ObjectiveSense::Minimize, h);
	//   // 	  block.M_->dataReport();
	//   // M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
	   block.M_->solve();
	  
	  
	  std::cout << "Solution : " << std::endl;
	  std::cout<<M->primalObjValue()  <<std::endl;
 
  // auto x=Matrix::dense(block.total_refs_.size(), 1, block.XT_->level());
		    // 	  for(auto it=block.total_refs_.begin(); it!=block.total_refs_.end();++it)
		    // {
		    //   std::cout <<"  "<<it->first<< "  " <<x->get(it->second.i1_,0)<< " index "<< it->second.i1_<<std::endl;
		      
		      
		    // }
    
}
