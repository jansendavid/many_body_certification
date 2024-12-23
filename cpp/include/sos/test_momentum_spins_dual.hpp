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
 basis_structure states=get_basis_2d(Lx, 2);

 
for(auto a: states)
{

	std::cout<<"sec "<< a.first<< " and size "<< a.second.size()<< " and "<<a.second.size()/(Lx*Lx)<<std::endl;
}

    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    auto basis =momentum_symmetry_solver_dual(Lx,states,M,"xyz");
  //   // for(auto a: basis.TI_map_)
  //   //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
     double J=1;
     double Delta=1.;
    
  //    auto b=define_xxz2d_sos( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
 
  //       basis.set_b(b);
  //       basis.fix_constrains();
  //      auto h=basis.get_costfunction();
	// //   std::cout<<C->toString()<<std::endl;
	  
	// // //   {std::pair<int,int> a(0,0);
	// // //   std::pair<int,int> b(1,0);
	// // // //generate_rmds_primal({a, b},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
	// // //   }
	// // //     {std::pair<int,int> a(0,0);
	// // //   std::pair<int,int> b(0,1);
	// // //    std::pair<int,int> c(0,2);
	// // // //generate_rmds_primal({a, b, c},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
	// // //   }
  //   basis.M_->objective(ObjectiveSense::Minimize, h);
	// 	  basis.M_->dataReport();
	//   M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
  //   basis.M_->solve();
	  
	  
  //   std::cout << "Solution : " << std::endl;
  //   std::cout<<std::setprecision(9)<<M->primalObjValue()<<std::endl;
	  
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

    basis_structure states=get_basis_2d(Lx, 3, -3, false);

 
for(auto a: states)
{

	std::cout<<"sec "<< a.first<< " and size "<< a.second.size()<< " and "<<a.second.size()/(Lx*Lx)<<std::endl;
}

    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
    auto basis =momentum_symmetry_solver_sos(Lx,states,M,"xyz");
  //   // for(auto a: basis.TI_map_)
  //   //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
     double J=1;
     double Delta=1.;
    
     auto b=define_xxz2d_sos( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
 std::cout<< "x"<<std::endl;
        basis.set_b(b);
        
        basis.fix_constrains();
       auto h=basis.get_costfunction();
	// //   std::cout<<C->toString()<<std::endl;
	  
	// // //   {std::pair<int,int> a(0,0);
	// // //   std::pair<int,int> b(1,0);
	// // // //generate_rmds_primal({a, b},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
	// // //   }
	// // //     {std::pair<int,int> a(0,0);
	// // //   std::pair<int,int> b(0,1);
	// // //    std::pair<int,int> c(0,2);
	// // // //generate_rmds_primal({a, b, c},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
	// // //   }

  std::cout<< "starting solving SDP"<<std::endl;
    basis.M_->objective(ObjectiveSense::Maximize, h);
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
void test_x()
{
  int Lx=4;
  int Ly=4;
 basis_structure states=get_basis_2d(Lx, 1,0);
 basis_structure states_2=get_basis_2d(Lx, 1, -1);
//  basis_structure states;
//   std::vector<op_vec> v_block_0;
//   std::vector<op_vec> v_block_1;
//   std::vector<op_vec> v_block_2;
//   std::vector<op_vec> v_block_3;
//   states.insert({0, v_block_0});
//   states.insert({1, v_block_1});
//   states.insert({2, v_block_2});
//   states.insert({3, v_block_3});

//   op_vec v0={spin_op("x", {0,0}, Lx),spin_op("x", {0, 3}, Lx)};
//     op_vec v1={spin_op("x", {0,0}, Lx),spin_op("x", {1, 3}, Lx)};
//  states[0].push_back(v0);
//   states[0].push_back(v1);
//    basis_structure states_2;
//   std::vector<op_vec> v_block_0_x;
//   std::vector<op_vec> v_block_1_x;
//   std::vector<op_vec> v_block_2_x;
//   std::vector<op_vec> v_block_3_x;
//   states_2.insert({0, v_block_0_x});
//   states_2.insert({1, v_block_1_x});
//   states_2.insert({2, v_block_2_x});
//   states_2.insert({3, v_block_3_x});
//    v0={spin_op("x", {0,0}, Lx),spin_op("x", {0, 3}, Lx)};
//     v1={spin_op("x", {0,0}, Lx),spin_op("x", {3, 1}, Lx)};
//  states_2[0].push_back(v0);
//   states_2[0].push_back(v1);

for(auto a: states)
{

	std::cout<<"sec "<< a.first<< " and size "<< a.second.size()<< " and "<<a.second.size()/(Lx*Lx)<<std::endl;
  for(auto x: a.second)
  {
    std::cout<<print_op(x)<<std::endl;
  }
}
for(auto a: states_2)
{

	std::cout<<"sec "<< a.first<< " and size "<< a.second.size()<< " and "<<a.second.size()/(Lx*Lx)<<std::endl;
   for(auto x: a.second)
  {
    std::cout<<print_op(x)<<std::endl;
  }
}


    Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
     Model::t M_2= new Model("sdo1"); auto _M_2 = finally([&]() { M_2->dispose(); });
    auto basis =momentum_symmetry_solver_sos(Lx,states,M,"xyz");
    std::cout<<"##########################################################################################"<<std::endl;
     auto basis_2 =momentum_symmetry_solver_sos(Lx,states_2,M_2,"xyz");
     int i=0;



// for(auto a: basis.total_refs_)
// {
//   bool in=false;
//  for(auto b: basis_2.total_refs_)
// {
//   if(a.first==b.first)
//   {

//     in=true;
 
//   }
// }
//   if(!in)
// {
// std::cout<<a.first<<std::endl;
//     i+=1;
// }

// } 

// std::cout<< "i "<<i<<std::endl;

//     Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
//     int tot=3;
// matrix_organizer C;
// matrix_organizer Ctilde;
// for(int i=0; i<tot; i++)
// {
// C.add_values({i,i},1);
// Ctilde.add_values({i,i},1);
// }

// matrix_organizer A;
// A.add_values({1,0},1);
// A.add_values({0,1},1);

// matrix_organizer B;
// B.add_values({2,0},1);
// B.add_values({0,2},1);
// std::vector<double> b;
// b.push_back(0);
// b.push_back(1);
// b.push_back(1);
// b.push_back(0);
// std::vector<Expression::t> Xs;
// auto x1=M->variable("1", Domain::inPSDCone(tot));
// Xs.push_back(Expr::neg(x1));
// auto x2=M->variable("2", Domain::inPSDCone(tot));
// Xs.push_back(Expr::neg(x2));

//  Expression::t expressions_1 =Expr::constTerm( 0.);
//     //expressions_1->index(0)=Expr::add(expressions_1->index(0), Expr::dot(A.make_matrix(tot,tot),(x1) ));
//     expressions_1=Expr::add(expressions_1, Expr::dot(A.make_matrix(tot,tot),(x1) ));

// //    expressions_1->index(1)=Expr::add(expressions_1->index(1), Expr::dot(B.make_matrix(tot,tot),(x2) ));

// //    expressions_1->index(1)=Expr::add(expressions_1->index(1), Expr::dot(B.make_matrix(tot,tot),(x2) ));
// M->constraint( expressions_1, Domain::equalsTo(b[1]));
// //M->constraint(  Expr::dot(A.make_matrix(tot,tot),(x1) ), Domain::equalsTo(b[1]));
// //M->constraint( expressions_1->index(1), Domain::equalsTo(b[2]));

//  Expression::t ee=Expr::constTerm(0.);
//         ee=Expr::sub(ee, Expr::dot(C.make_matrix(tot,tot),(x1)));
//         //ee=Expr::add(ee, Expr::dot(Ctilde.make_matrix(tot,tot),(Xs[1])));

//   //   auto basis =momentum_symmetry_solver_sos(Lx,states,M,"xyz");

//   //    double J=1;
//   //    double Delta=1.;
    
//   //    auto b=define_xxz2d_sos( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
 
//   //       basis.set_b(b);
//   //       basis.fix_constrains();
//   //      auto h=basis.get_costfunction();
// 	// // //   std::cout<<C->toString()<<std::endl;
	  

//     M->objective(ObjectiveSense::Maximize, ee);
// 		 M->dataReport();
// 	  M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
//    M->solve();
	  
	  
//     std::cout << "Solution : " << std::endl;
//     std::cout<<std::setprecision(9)<<M->primalObjValue()<<std::endl;
	
	    return;
}
void test_y()
{
  int Lx=4;
  //basis_structure states=get_basis_1d(Lx, int(Lx/2),0);
   basis_structure states=get_basis_2d(Lx, 3, -3);

for(auto a: states)
{

 	std::cout<<"sec "<< a.first<< " and size "<< a.second.size()<< " and "<<a.second.size()/(Lx)<<std::endl;

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
