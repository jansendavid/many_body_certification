#pragma once
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>
#include "spins.hpp"
#include <unordered_map>
#include <Eigen/Dense>
#include "sos/complex_momentum_dual.hpp"
#include "spin_hamiltonians_TIsym.hpp"
#include "functions.hpp"
#include "reduced_dms.hpp"
using namespace mosek::fusion;
using namespace monty;

void test_multiple_blocks_bounding_observables_2d_rdm_sos()
{
  std::cout << "WARNING! Takes a lot of memory" << std::endl;
  int Lx = 4;
  int Ly = 4;
  auto map_sec = get_sector_map();
  basis_structure states = get_states();
  get_order_one_monomials(states, map_sec, Lx, true);
  get_order_two_monomials(states, map_sec, Lx, 3, -3, true);
  get_order_three_monomials(states, map_sec, Lx, true);

  get_order_four_monomials(states, map_sec, Lx, true);

  double E_upper = -0.7017777;
  double E_lower = -0.70305078;
  for (auto a : states)
  {

    std::cout << "sec " << a.first << " and size " << a.second.size() << " and " << a.second.size() / (Lx * Lx) << std::endl;
    for (auto l : a.second)
    {
      //	std::cout<< print_op(l)<<std::endl;
    }
  }
  auto data = get_rdms(Lx, 2);

  rdms_struct rdms(data);
  Model::t M = new Model("sdo1");
  auto _M = finally([&]()
                    { M->dispose(); });
  bool maximize = false;
  bool bilayer = false;
  auto basis = momentum_symmetry_solver_sos(Lx, states, M, rdms, "xyz", bilayer, maximize);

  // for(auto a: basis.TI_map_)
  //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
  double J = 1;
  double Delta = 1.;
  std::cout << "done 1" << std::endl;
  auto energy_vec = define_xxz2d_sos(basis.total_refs_, basis.TI_map_, J, Delta, Ly, Lx);
  std::cout << "done 2" << std::endl;

  basis.set_energy_vec(energy_vec, E_upper, E_lower);
  std::cout << "done 3" << std::endl;
  std::pair<std::vector<int>, std::vector<int>> pos{{0, 0}, {2, 2}};
  std::pair<std::string, std::string> dirs{"x", "x"};

  auto corr_func = define_correlation_function_sos(basis.total_refs_, basis.TI_map_, dirs, pos, Ly, Lx);
  // define_xxz2d_sos( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
  //
  basis.set_b(corr_func);
  std::cout << "done 3" << std::endl;
  basis.fix_constrains();

  auto h = basis.get_costfunction();

  if (maximize)
  {
    basis.M_->objective(ObjectiveSense::Maximize, h);
  }
  else
  {
    basis.M_->objective(ObjectiveSense::Minimize, h);
  }
  basis.M_->dataReport();
  M->setLogHandler([=](const std::string &msg)
                   { std::cout << msg << std::flush; });
  basis.M_->solve();

  std::cout << "Solution : " << std::endl;
  std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;

  //  //double sol=M->primalObjValue(); //-0.720847431

  //   //	   if(std::abs(sol+0.44670126)>1e-06)
  //   // {std::cout<<"error, not converging properly"<<std::endl;}
  //     return;
}

void test_multiple_blocks_bounding_observables_2d_rdm()
{
  std::cout << "WARNING! Takes a lot of memory" << std::endl;
  int Lx = 4;
  int Ly = 4;
  auto map_sec = get_sector_map();
  basis_structure states = get_states();
  get_order_one_monomials(states, map_sec, Lx, true);
  get_order_two_monomials(states, map_sec, Lx, 3, -3, true);
  // get_order_three_monomials(states, map_sec, Lx, true);

  // get_order_four_monomials(states, map_sec, Lx, true);

  double E_upper = -0.7017777;
  double E_lower = -0.70305078;
  for (auto a : states)
  {

    std::cout << "sec " << a.first << " and size " << a.second.size() << " and " << a.second.size() / (Lx * Lx) << std::endl;
    for (auto l : a.second)
    {
      //	std::cout<< print_op(l)<<std::endl;
    }
  }
  auto data = get_rdms(Lx, 4);

  rdms_struct rdms(data);
  Model::t M = new Model("sdo1");
  auto _M = finally([&]()
                    { M->dispose(); });
  auto basis = momentum_symmetry_solver_dual(Lx, states, M, rdms, "xyz");

  // for(auto a: basis.TI_map_)
  //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
  double J = 1;
  double Delta = 1.;
  std::cout << "done 1" << std::endl;
  auto energy_vec = define_xxz2d_sos(basis.total_refs_, basis.TI_map_, J, Delta, Ly, Lx);
  std::cout << "done 2" << std::endl;

  basis.set_energy_vec(energy_vec, E_upper, E_lower);
  std::cout << "done 3" << std::endl;
  std::pair<std::vector<int>, std::vector<int>> pos{{0, 0}, {2, 2}};
  std::pair<std::string, std::string> dirs{"x", "x"};
  auto corr_func = define_correlation_function_sos(basis.total_refs_, basis.TI_map_, dirs, pos, Ly, Lx);
  // define_xxz2d_sos( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
  //
  basis.set_b(corr_func);
  std::cout << "done 3" << std::endl;
  basis.fix_constrains();

  auto h = basis.get_costfunction();

  //  basis.M_->objective(ObjectiveSense::Maximize, h);
  basis.M_->objective(ObjectiveSense::Minimize, h);
  basis.M_->dataReport();
  M->setLogHandler([=](const std::string &msg)
                   { std::cout << msg << std::flush; });
  basis.M_->solve();

  std::cout << "Solution : " << std::endl;
  std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;

  //  //double sol=M->primalObjValue(); //-0.720847431

  //   //	   if(std::abs(sol+0.44670126)>1e-06)
  //   // {std::cout<<"error, not converging properly"<<std::endl;}
  //     return;
}
void test_multiple_blocks_higher_order_2d_rdm()
{
  std::cout << "WARNING! Takes a lot of memory" << std::endl;
  int Lx = 4;
  int Ly = 4;
  auto map_sec = get_sector_map();
  basis_structure states = get_states();
  get_order_one_monomials(states, map_sec, Lx, true);
  get_order_two_monomials(states, map_sec, Lx, 1, -1, true);
  // get_order_three_monomials(states, map_sec, Lx, true);

  // get_order_four_monomials(states, map_sec, Lx, true);

  for (auto a : states)
  {

    std::cout << "sec " << a.first << " and size " << a.second.size() << " and " << a.second.size() / (Lx * Lx) << std::endl;
    for (auto l : a.second)
    {
      //	std::cout<< print_op(l)<<std::endl;
    }
  }
  auto data = get_rdms(Lx, Lx);

  rdms_struct rdms(data);
  Model::t M = new Model("sdo1");
  auto _M = finally([&]()
                    { M->dispose(); });
  auto basis = momentum_symmetry_solver_dual(Lx, states, M, rdms, "xyz");
  std::cout << "here" << std::endl;

  // for(auto a: basis.TI_map_)
  //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
  double J = 1;
  double Delta = 1.;

  auto b = define_xxz2d_sos(basis.total_refs_, basis.TI_map_, J, Delta, Ly, Lx);

  basis.set_b(b);
  basis.fix_constrains();
  auto h = basis.get_costfunction();

  // // //     {std::pair<int,int> a(0,0);
  // // //   std::pair<int,intget_basis_2d(Lx, 3, -3, true);> b(0,1);
  // // //    std::pair<int,int> c(0,2);
  // // // //generate_rmds_primal({a, b, c},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
  // // //   }
  basis.M_->objective(ObjectiveSense::Minimize, h);
  basis.M_->dataReport();
  M->setLogHandler([=](const std::string &msg)
                   { std::cout << msg << std::flush; });
  basis.M_->solve();

  std::cout << "Solution : " << std::endl;
  std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;
}

void test_multiple_blocks_higher_order_2d_rdm_sos()
{
  std::cout << "WARNING! Takes a lot of memory" << std::endl;
  int Lx = 4;
  int Ly = 4;
  auto map_sec = get_sector_map();
  basis_structure states = get_states();
  get_order_one_monomials(states, map_sec, Lx, true);
  get_order_two_monomials(states, map_sec, Lx, 3, -3, true);
  get_order_three_monomials(states, map_sec, Lx, true);
  basis_structure states_2;
  states_2[0] = states[0];
  states_2[1] = states[1];
  states = states_2;
  // get_order_four_monomials(states, map_sec, Lx, true);

  auto data = get_rdms(Lx, Lx);
  // translation_invariant_rdms_4th(Lx, Ly);
  //   rdms_struct data;
  //   rdm_operator newstate({{0,0}, {1,1}});
  //   data.add_operator(newstate);
  rdms_struct rdms(data);

  Model::t M = new Model("sdo1");
  auto _M = finally([&]()
                    { M->dispose(); });
  auto basis = momentum_symmetry_solver_sos(Lx, states, M, rdms, "xyz");
  //     // for(auto a: basis.TI_map_)
  //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
  //   for(auto k: basis.total_refs_)
  //   {

  // //	std::cout<<k.first<<std::endl;
  //   }
  double J = 1;
  double Delta = 1.;

  auto b = define_xxz2d_sos(basis.total_refs_, basis.TI_map_, J, Delta, Ly, Lx);

  basis.set_b(b);

  basis.fix_constrains();
  auto h = basis.get_costfunction();
  // //   std::cout<<C->toString()<<std::endl;

  std::cout << "starting solving SDP" << std::endl;
  basis.M_->objective(ObjectiveSense::Maximize, h);
  basis.M_->dataReport();
  M->setLogHandler([=](const std::string &msg)
                   { std::cout << msg << std::flush; });
  basis.M_->solve();
  auto cons = M->getConstraint(0)->dual();
  std::cout << cons << std::endl;

  std::cout << "Solution : " << std::endl;
  std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;
  // thord order -0.703778105
  return;
}

// void test_y()
// {
//   int L=4;
//  op_vec vec={spin_op("x", {2,2}, {L,L}),spin_op("x", {0,0}, {L,L})};

// std::cout<< print_op(vec)<<std::endl;
//  auto [fac, nf] =get_normal_form(vec);
//  std::cout<< fac << " and "<< print_op(nf)<<std::endl;

// 	    return;
// }

// void test_d8_symm()
// {
//   int L=6;
// //  op_vec vec={spin_op("x", {0,0}, {L,L}),spin_op("x", {1,1}, {L,L}),spin_op("y", {2,1}, {L,L})};
// // std::cout<< print_op(vec)<<std::endl;
// // auto all_d8=generate_all_d8(vec,  L);
// // std::cout<< "start "<<std::endl;
// // for(auto b: all_d8)
// // {
// // std::cout<<print_op(b)<<std::endl;

// // }
// //s_[x,(1,0)]s_[y,(2,0)]s_[x,(4,3)]s_[y,(5,3)]
// // s_[x,(1,0)]s_[y,(2,0)]s_[x,(4,3)]s_[y,(5,3)]
// // s_[x,(0,5)]s_[y,(0,4)]s_[x,(3,2)]s_[y,(3,1)]
// // s_[x,(5,0)]s_[y,(4,0)]s_[x,(2,3)]s_[y,(1,3)]
// // s_[x,(0,1)]s_[y,(0,2)]s_[x,(3,4)]s_[y,(3,5)]
// // s_[x,(5,0)]s_[y,(4,0)]s_[x,(2,3)]s_[y,(1,3)]
// // s_[x,(0,1)]s_[y,(0,2)]s_[x,(3,4)]s_[y,(3,5)]
// // s_[x,(1,0)]s_[y,(2,0)]s_[x,(4,3)]s_[y,(5,3)]
// // s_[x,(0,5)]s_[y,(0,4)]s_[x,(3,2)]s_[y,(3,1)]
// op_vec vec={spin_op("x", {1,0}, {L,L}),spin_op("y", {2,0}, {L,L}),spin_op("x", {4,3}, {L,L}),spin_op("y", {5,3}, {L,L})};
// std::cout<<print_op(vec)<<std::endl;
// auto all_d8=generate_all_d8(vec,  L);
// for(auto b: all_d8)
// {
// std::cout<<print_op(b)<<std::endl;

// }
// 	    return;
// }
//  void test_primal_and_dual()
//  {
//   int Lx=4;
//   int Ly=4;
//     auto map_sec=get_sector_map();
//     basis_structure states1=get_states();
//     get_order_one_monomials(states1, map_sec, Lx, true);
//      get_order_two_monomials(states1, map_sec, Lx,3, -3,true);
//    basis_structure states2=get_states();
//     get_order_one_monomials(states2, map_sec, Lx, true);
//      get_order_two_monomials(states2, map_sec, Lx,3, -3,true);

//          Model::t M1 = new Model("sdo1"); auto _M1 = finally([&]() { M1->dispose(); });

//    Model::t M2 = new Model("sdo1"); auto _M2 = finally([&]() { M2->dispose(); });

//     auto basis1 =momentum_symmetry_solver_sos(Lx,states1,M1,"xyz");
//  auto basis2 =momentum_symmetry_solver_dual(Lx,states2,M2,"xyz");
//   double J=1;
//      double Delta=1.;

//      auto b1=define_xxz2d_sos( basis1.total_refs_,basis1.TI_map_, J, Delta, Ly, Lx);

//         basis1.set_b(b1);

//         basis1.fix_constrains();
//        auto h1=basis1.get_costfunction();

//           auto b2=define_xxz2d_sos( basis2.total_refs_,basis2.TI_map_, J, Delta, Ly, Lx);

//         basis2.set_b(b2);

//         basis2.fix_constrains();
//        auto h2=basis2.get_costfunction();

//         basis1.M_->objective(ObjectiveSense::Maximize, h1);
// 		  basis1.M_->dataReport();
// 	  M1->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
//     basis1.M_->solve();
// 	  //auto cons=M->getConstraint(0)->dual();
//     //std::cout<<cons<<std::endl;
// 	   basis2.M_->objective(ObjectiveSense::Minimize, h2);
// 		  basis2.M_->dataReport();
// 	  M2->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
//     basis2.M_->solve();

//     std::cout << "Solution : " << std::endl;
//     std::cout<<std::setprecision(9)<<M1->primalObjValue()<<std::endl;

//     std::cout << "Solution : " << std::endl;
//     std::cout<<std::setprecision(9)<<M2->primalObjValue()<<std::endl;
//     auto y_fusion=Matrix::dense(basis2.total_refs_.size(), 1, basis2.y_->level());
//     int i=0;
//     double max=0.;
//     for(auto val : basis2.total_refs_)
//     {
//       if(val.first!="1" and val.first!="0")
//       {
//         //auto c_fusion=Matrix::dense(basis2.total_refs_.size(), 1, basis2.y_->level());
//    //  std::cout<< (*(basis2.y_->level()))[val.second]<< " and "<< -1.*(*(M1->getConstraint(i)->dual()))[0]<<std::endl;
//      max=std::max(max, std::abs((*(basis2.y_->level()))[val.second]-(-1.*(*(M1->getConstraint(i)->dual()))[0])));
//      i++;
//       }
//    }
// //std::cout<< max << " and "<<std::abs(M1->primalObjValue()-M2->primalObjValue())<<std::endl;
// assert(max<1e-3);
// assert(std::abs(M1->primalObjValue()-M2->primalObjValue())<1e-6);

//   return;
//  }
void test_multiple_blocks_higher_order_J1J22d_rdm_sos()
{
  std::cout << "WARNING! Takes a lot of memory" << std::endl;
  int Lx = 4;
  int Ly = 4;
  auto map_sec = get_sector_map();
  basis_structure states = get_states();
  get_order_one_monomials(states, map_sec, Lx, true);
  get_order_two_monomials(states, map_sec, Lx, 3, -3, true);
  get_order_three_monomials(states, map_sec, Lx, true);

  get_order_four_monomials(states, map_sec, Lx, true);
  auto data = get_rdms(Lx, Lx);
  // translation_invariant_rdms_4th(Lx, Ly);
  //   rdms_struct data;
  //   rdm_operator newstate({{0,0}, {1,1}});
  //   data.add_operator(newstate);
  rdms_struct rdms(data);

  for (auto a : states)
  {

    std::cout << "sec " << a.first << " and size " << a.second.size() << " and " << a.second.size() / (Lx * Lx) << std::endl;

    // for(auto n:a.second)
    //  {std::cout<<print_op(n)<<std::endl;}
    std::cout << "###################################" << std::endl;
  }

  Model::t M = new Model("sdo1");
  auto _M = finally([&]()
                    { M->dispose(); });
  auto basis = momentum_symmetry_solver_sos(Lx, states, M, rdms, "xyz");
  //     // for(auto a: basis.TI_map_)
  //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
  //   for(auto k: basis.total_refs_)
  //   {

  // //	std::cout<<k.first<<std::endl;
  //   }
  double J1 = 1;
  double J2 = 1.;

  auto b = define_J1J22d_sos(basis.total_refs_, basis.TI_map_, J1, J2, Ly, Lx);

  basis.set_b(b);

  basis.fix_constrains();
  auto h = basis.get_costfunction();
  // //   std::cout<<C->toString()<<std::endl;

  std::cout << "starting solving SDP" << std::endl;
  basis.M_->objective(ObjectiveSense::Maximize, h);
  basis.M_->dataReport();
  M->setLogHandler([=](const std::string &msg)
                   { std::cout << msg << std::flush; });
  basis.M_->solve();
  auto cons = M->getConstraint(0)->dual();
  std::cout << cons << std::endl;

  std::cout << "Solution : " << std::endl;
  std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;

  return;
}
