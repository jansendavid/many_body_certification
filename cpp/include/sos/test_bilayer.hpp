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
void test_pos()
{
  int Lx=4;
  int Ly=4;
  int layers=2;
  std::vector<int> pos;
 
//   for(int i=0; i<Lx; i++)
//   {

//      for(int j=0; j<Ly; j++)
//   {

//     for(int n=0; n<layers; n++)
//   {

//     op_vec vec={spin_op("x", {n,i,j}, {layers, Lx,Lx})};
//     int cnt = count(pos.begin(), pos.end(), vec[0].pos());
// if(cnt>0)
// {
//   std::cout<< "mapping does not work"<<std::endl;
// }
//     pos.push_back(vec[0].pos());
//   }  
//   }

//   }
// input s_[x,(0,0,0)]s_[z,(0,0,1)]s_[z,(0,3,1)]s_[x,(1,3,1)]
// s_[x,(0,2,0)]s_[z,(0,2,1)]s_[z,(0,1,1)]s_[x,(1,1,1)] fac (0,1) s_[x,(0,2,0)]s_[z,(0,1,1)]s_[y,(0,2,1)]
//s_[x,(0,0,0)]s_[y,(0,3,0)]s_[x,(0,0,1)]s_[y,(1,3,1)]
op_vec op={spin_op("x", {0,0,0}, {layers, Lx,Lx}),spin_op("y", {0,3,0}, {layers, Lx,Lx}),spin_op("x", {0,0,1}, {layers, Lx,Lx}),spin_op("y", {1,3,1}, {layers, Lx,Lx})};
std::cout<< print_op(op)<<std::endl;
int inc=1;
  for(int i=1; i<Lx; i+=inc)
    {
      auto new_op=translation_y(op, i,Lx);
      //      auto coeff=bubbleSort(new_op, new_op.size());
      auto [fac, vec] =get_normal_form(new_op);
      for(auto l: vec)
      {
        std::cout<< l.pos()<<std::endl;
      }
  std::cout<< fac << " and "<< print_op(vec)<<std::endl;
    }
  // auto [fac, nf] =get_normal_form(vec);

  //  generate_all_translations_y(vec,  Lx, 1);

}



void test_bilayer_sos()
{
  std::cout<< "WARNING! Takes a lot of memory"<<std::endl;
  int Lx=4;
  int Ly=4;
  int layers=2;

    basis_structure states=get_basis_bilayer_2d(Lx,Ly, layers, 2, 0,true);

   for(auto a: states)
     {std::cout<< a.first << " -> "<<a.second.size<<std::endl;}



    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
   
 auto basis =momentum_symmetry_solver_sos(Lx,states,M,"xyz");
  
      double J=1;
      double Delta=1.;
    double J_perpendicular=1;
    double J_parallel=1.8;
     double J_x=0.6;
      auto b=define_heisenberg_bilayer_sos( basis.total_refs_,basis.TI_map_, J_perpendicular, J_parallel, J_x , Ly, Lx, layers);
//  std::cout<< "x"<<std::endl;
        basis.set_b(b);
        
        basis.fix_constrains();
       auto h=basis.get_costfunction();
// 	// //   std::cout<<C->toString()<<std::endl;
	  
// 	// // //   {std::pair<int,int> a(0,0);
// 	// // //   std::pair<int,int> b(1,0);
// 	// // // //generate_rmds_primal({a, b},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
// 	// // //   }
// 	// // //     {std::pair<int,int> a(0,0);
// 	// // //   std::pair<int,int> b(0,1);
// 	// // //    std::pair<int,int> c(0,2);
// 	// // // //generate_rmds_primal({a, b, c},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
// 	// // //   }

  std::cout<< "starting solving SDP"<<std::endl;
    basis.M_->objective(ObjectiveSense::Maximize, h);
		  basis.M_->dataReport();
	  M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
    basis.M_->solve();
	  
	  
//     std::cout << "Solution : " << std::endl;
    std::cout<<std::setprecision(9)<<M->primalObjValue()<<std::endl;
	  
//    //double sol=M->primalObjValue(); //-0.720847431
	  

// 	  //	   if(std::abs(sol+0.44670126)>1e-06)
// 	  // {std::cout<<"error, not converging properly"<<std::endl;}
	    return;
}
