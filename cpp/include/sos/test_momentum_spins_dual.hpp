#pragma once
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"

//#include "spins.hpp"
#include <unordered_map>
#include <Eigen/Dense>
#include "sos/complex_momentum_dual.hpp"
#include "spin_hamiltonians_TIsym.hpp"
#include "functions.hpp"
//#include "reduced_dms.hpp"
#include "symmetries.hpp"
#include "lattices.hpp"
using namespace mosek::fusion;
using namespace monty;
void test_counting()
{
  // making sur ethe spin counting is working
  // {
  //   int Lx = 4;
  //   int Ly = 4;
  //   int layer = 2;
  //   std::vector<int> offset_vector = {layer, Ly, Lx};
  //   std::set<int> exists;
  //   for (int i = 0; i < layer; i++)
  //   {
  //     for (int j = 0; j < Ly; j++)
  //     {

  //       for (int k = 0; k < Lx; k++)
  //       {
  //         auto v_p = spin_op("x", {i, j, k}, offset_vector);
  //         exists.insert(v_p.pos());
  //       }
  //     }
  //   }

  //   assert(exists.size() == Lx * Ly * layer);
  // }

  // {
  //   int Lx = 50;
  //   int Ly = 1;
  //   int layer = 1;
  //   std::vector<int> offset_vector = {layer, Ly, Lx};
  //   std::set<int> exists;
  //   for (int i = 0; i < layer; i++)
  //   {
  //     for (int j = 0; j < Ly; j++)
  //     {

  //       for (int k = 0; k < Lx; k++)
  //       {
  //         auto v_p = spin_op("x", {i, j, k}, offset_vector);
  //         exists.insert(v_p.pos());
  //       }
  //     }
  //   }

  //   assert(exists.size() == Lx * Ly * layer);
  // }

  // {
  //   int Lx = 12;
  //   int Ly = 1;
  //   int layer = 4;
  //   std::vector<int> offset_vector = {layer, Ly, Lx};
  //   std::set<int> exists;
  //   for (int i = 0; i < layer; i++)
  //   {
  //     for (int j = 0; j < Ly; j++)
  //     {

  //       for (int k = 0; k < Lx; k++)
  //       {
  //         auto v_p = spin_op("x", {i, j, k}, offset_vector);
  //         exists.insert(v_p.pos());
  //       }
  //     }
  //   }

  //   assert(exists.size() == Lx * Ly * layer);
  // }
}
// void test_multiple_blocks_bounding_observables_2d_rdm_sos()
// {
//   std::cout << "WARNING! Takes a lot of memory" << std::endl;
//   int Lx = 4;
//   int Ly = 4;
//   auto lattice = SquareLattice(Ly, Lx, true, false, "xyz", "xyz");
//   auto map_sec = get_sector_map();
//   basis_structure states = get_states();
//   get_order_one_monomials(states, map_sec, Ly, Lx, true);
//   get_order_two_monomials(states, map_sec, Ly, Lx, 3, 3, -3, -3, true);
//   get_order_three_monomials(states, map_sec, Ly, Lx, true);

//   get_order_four_monomials(states, map_sec, Ly, Lx, true);

//   double E_upper = -0.7017777;
//   double E_lower = -0.70305078;
//   for (auto a : states)
//   {

//     std::cout << "sec " << a.first << " and size " << a.second.size() << " and " << a.second.size() / (Lx * Lx) << std::endl;
//     for (auto l : a.second)
//     {
//       //	std::cout<< print_op(l)<<std::endl;
//     }
//   }
//   auto data = get_rdms(Lx, 2);

//   rdms_struct rdms(data);
//   Model::t M = new Model("sdo1");
//   auto _M = finally([&]()
//                     { M->dispose(); });
//   bool maximize = false;
//   bool bilayer = false;
//   auto basis = momentum_symmetry_solver_sos(lattice, states, M, rdms, maximize);

//   // for(auto a: basis.TI_map_)
//   //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
//   double J = 1;
//   double Delta = 1.;
//   std::cout << "done 1" << std::endl;
//   auto energy_vec = define_xxz2d_sos(basis.total_refs_, lattice, J, Delta);
//   std::cout << "done 2" << std::endl;

//   basis.set_energy_vec(energy_vec, E_upper, E_lower);
//   std::cout << "done 3" << std::endl;
//   std::pair<std::vector<int>, std::vector<int>> pos{{0, 0}, {2, 2}};
//   std::pair<std::string, std::string> dirs{"x", "x"};

//   auto corr_func = define_correlation_function_sos(basis.total_refs_, lattice, dirs, pos);
//   // define_xxz2d_sos( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
//   //
//   basis.set_b(corr_func);
//   std::cout << "done 3" << std::endl;
//   basis.fix_constrains();

//   auto h = basis.get_costfunction();

//   if (maximize)
//   {
//     basis.M_->objective(ObjectiveSense::Maximize, h);
//   }
//   else
//   {
//     basis.M_->objective(ObjectiveSense::Minimize, h);
//   }
//   basis.M_->dataReport();
//   M->setLogHandler([=](const std::string &msg)
//                    { std::cout << msg << std::flush; });
//   basis.M_->solve();

//   std::cout << "Solution : " << std::endl;
//   std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;

//   //   //	   if(std::abs(sol+0.44670126)>1e-06)
//   //   // {std::cout<<"error, not converging properly"<<std::endl;}
//   //     return;
// }

// // void test_multiple_blocks_bounding_observables_2d_rdm()
// // {
// //   std::cout << "WARNING! Takes a lot of memory" << std::endl;
// //   int Lx = 4;
// //   int Ly = 4;
// //   auto map_sec = get_sector_map();
// //   basis_structure states = get_states();
// //   get_order_one_monomials(states, map_sec, Lx, true);
// //   get_order_two_monomials(states, map_sec, Lx, 3, -3, true);
// //   // get_order_three_monomials(states, map_sec, Lx, true);

// //   // get_order_four_monomials(states, map_sec, Lx, true);

// //   double E_upper = -0.7017777;
// //   double E_lower = -0.70305078;
// //   for (auto a : states)
// //   {

// //     std::cout << "sec " << a.first << " and size " << a.second.size() << " and " << a.second.size() / (Lx * Lx) << std::endl;
// //     for (auto l : a.second)
// //     {
// //       //	std::cout<< print_op(l)<<std::endl;
// //     }
// //   }
// //   auto data = get_rdms(Lx, 4);

// //   rdms_struct rdms(data);
// //   Model::t M = new Model("sdo1");
// //   auto _M = finally([&]()
// //                     { M->dispose(); });
// //   auto basis = momentum_symmetry_solver_dual(Lx, states, M, rdms, "xyz");

// //   // for(auto a: basis.TI_map_)
// //   //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
// //   double J = 1;
// //   double Delta = 1.;
// //   std::cout << "done 1" << std::endl;
// //   auto energy_vec = define_xxz2d_sos(basis.total_refs_, basis.TI_map_, J, Delta, Ly, Lx);
// //   std::cout << "done 2" << std::endl;

// //   basis.set_energy_vec(energy_vec, E_upper, E_lower);
// //   std::cout << "done 3" << std::endl;
// //   std::pair<std::vector<int>, std::vector<int>> pos{{0, 0}, {2, 2}};
// //   std::pair<std::string, std::string> dirs{"x", "x"};
// //   auto corr_func = define_correlation_function_sos(basis.total_refs_, basis.TI_map_, dirs, pos, Ly, Lx);
// //   // define_xxz2d_sos( basis.total_refs_,basis.TI_map_, J, Delta, Ly, Lx);
// //   //
// //   basis.set_b(corr_func);
// //   std::cout << "done 3" << std::endl;
// //   basis.fix_constrains();

// //   auto h = basis.get_costfunction();

// //   //  basis.M_->objective(ObjectiveSense::Maximize, h);
// //   basis.M_->objective(ObjectiveSense::Minimize, h);
// //   basis.M_->dataReport();
// //   M->setLogHandler([=](const std::string &msg)
// //                    { std::cout << msg << std::flush; });
// //   basis.M_->solve();

// //   std::cout << "Solution : " << std::endl;
// //   std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;

// //   //  //double sol=M->primalObjValue(); //-0.720847431

// //   //   //	   if(std::abs(sol+0.44670126)>1e-06)
// //   //   // {std::cout<<"error, not converging properly"<<std::endl;}
// //   //     return;
// // }
template <typename K, typename V>
std::vector<K> keys_with_value(const std::map<K, V> &m, std::string value)
{
  std::vector<K> keys;
  for (const auto &[k, v] : m)
  {
    if (v.first == value)
    {
      keys.push_back(k);
    }
  }
  return keys;
}
void test_multiple_blocks_higher_order_2d_rdm()
{

  // std::cout << "xxx WARNING! Takes a lot of memory" << std::endl;
  // int Lx = 4;
  // int Ly = 4;
  // // op_vec v_p = {spin_op("x", {3, 0}, {Lx, Ly}), spin_op("x", {0, 0}, {Lx, Ly}), spin_op("y", {1, 2}, {Lx, Ly}), spin_op("y", {1, 1}, {Lx, Ly})};
  // // auto b = all_translations(v_p, Lx, Ly);
  // // std::cout << "start  " << print_op(v_p) << std::endl;
  // // for (auto g : b)
  // // {
  // //   std::cout << print_op(g) << std::endl;
  // // }
  // // // std::cout << print_op(v0) << std::endl;
  // // // auto [coeff, nf] = get_normal_form(v0);
  // // // std::cout << print_op(nf) << std::endl;
  // auto lattice = SquareLattice(Ly, Lx, true, false, "xyz", "xyz");
  // auto map_sec = get_sector_map();
  // basis_structure states = get_states();
  // get_order_one_monomials(states, map_sec, Ly, Lx, true);

  // int r = 1;
  // get_order_two_monomials(states, map_sec, Ly, Lx, r, r, -r, -r, true);
  // // get_order_three_monomials(states, map_sec, Ly, Lx, true);
  // basis_structure states_2;
  // states_2[0] = states[0];
  // states_2[1] = states[1];
  // lattice.states_ = states_2;

  // // for (auto s : states_2)
  // // {
  // //   std::cout << s.first << "with size " << s.second.size() << std::endl;
  // //   for (auto m : s.second)
  // //   {
  // //     std::cout << print_op(m) << std::endl;
  // //   }
  // // }
  // // // for (auto a : lattice.TI_map_)
  // // // {
  // // //   std::cout << a.first << "  " << a.second.first << std::endl;
  // // // }
  // // // get_order_four_monomials(states, map_sec, Ly, Lx, true);

  // // // //   for (auto a : states)
  // // // //   {

  // // // //     std::cout << "sec " << a.first << " and size " << a.second.size() << " and " << a.second.size() / (Lx * Lx) << std::endl;
  // // // //     for (auto l : a.second)
  // // // //     {
  // // // //       //	std::cout<< print_op(l)<<std::endl;
  // // // //     }
  // // // //   }
  // auto data = get_rdms(Lx, Lx);
  // data = {};
  // rdms_struct rdms(data);

  // Model::t M = new Model("sdo1");
  // auto _M = finally([&]()
  //                   { M->dispose(); });
  // std::cout << "start " << std::endl;
  // auto basis = momentum_symmetry_solver_dual<SquareLattice>(lattice, M, rdms);
  // // for (auto a : lattice.variable_map_)
  // // {
  // //   std::cout << " " << a.first << std::endl;
  // // }

  // // auto g = keys_with_value(lattice.TI_map_, "s_[x,(0,0)]s_[z,(1,0)]s_[z,(0,1)]s_[x,(3,3)]");
  // // std::cout << g.size() << std::endl;
  // // for (auto g_ : g)
  // // {
  // //   std::cout << g_ << std::endl;
  // // }
  // // auto it = lattice.TI_map_.find("s_[x,(0,0)]s_[z,(1,0)]s_[z,(0,1)]s_[x,(3,3)]");
  // // if (it != lattice.TI_map_.end())
  // // {
  // //   std::cout << "found and key was " << it->second.first << std::endl;
  // // }

  // // auto it2 = lattice.variable_map_.find("s_[x,(0,0)]s_[z,(1,0)]s_[z,(0,1)]s_[x,(3,3)]");
  // // if (it2 != lattice.variable_map_.end())
  // // {
  // //   std::cout << "found 2 " << std::endl;
  // // }
  // // // // // // for(auto a: basis.TI_map_))
  // double J = 1;
  // double Delta = 1.;

  // auto b = define_xxz2d_sos(lattice, J, Delta);

  // basis.set_b(b);
  // basis.fix_constrains();
  // auto h = basis.get_costfunction();

  // // // // //     {std::pair<int,int> a(0,0);
  // // // //   std::pair<int,intget_basis_2d(Lx, 3, -3, true);> b(0,1);
  // // // //    std::pair<int,int> c(0,2);
  // // // // //generate_rmds_primal({a, b, c},basis.total_refs_,basis.TI_map_ , basis.variables_, Lx, M);
  // // // //   }
  // basis.M_->objective(ObjectiveSense::Minimize, h);
  // basis.M_->dataReport();
  // M->setLogHandler([=](const std::string &msg)
  //                  { std::cout << msg << std::flush; });
  // basis.M_->solve();

  // std::cout << "Solution : " << std::endl;
  // std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;
  //-0.702827317
}
std::vector<int> get_sign(op_vec op)
{
  std::vector<int> res;
  std::vector<std::string> dirs = {"x", "y", "z"};
  for (auto dir_ : dirs)
  {
    int fac = 1;
    for (auto a : op)
    {

      if (a.get_dir() == dir_)
      {
        fac *= -1;
      }
    }    
    res.push_back(fac);
}

  return res;
}
void test_multiple_blocks_higher_order_2d_rdm_sos()
{
   int Lx = 4;
   int Ly = 4;
 
  // std::cout << "WARNING! Takes a lot of memory" << std::endl;

  // op_vec op = {spin_op("z", {3, 2}, {Lx, Ly}), spin_op("z", {3, 3}, {Lx, Ly}), spin_op("z", {1, 0}, {Lx, Ly}), spin_op("z", {0, 0}, {Lx, Ly})};
  // auto all_p = generate_all_permutations_xyz(op);
  // std::cout<<all_p.size()<<std::endl;
  // for(auto p: all_p)
  // {
  //   std::cout<<print_op(p)<<std::endl;
  // }
  // std::cout << "start " << print_op(op) << std::endl;
  // for (auto a : all_p)
  // {
  //   std::cout << print_op(a) << std::endl;
  // }
   auto lattice = SquareLattice(Ly, Lx, true, false, "xyz", "xyz");


  auto map_sec = get_sector_map();
  basis_structure states = get_states(1);

std::cout<<states.size()<<std::endl;
for(auto b: states)
{
  std::cout<<b.second.size()<<std::endl;
}
   get_order_one_monomials(states, map_sec, Ly, Lx, true);
  int r = 1;
  get_order_two_monomials(states, map_sec, Ly, Lx, r, r, -r, -r, true);
    // get_order_three_monomials(states, map_sec, Ly, Lx, true);

// for(auto a: states)
// {
//   std::cout<<"SS "<< a.first<<std::endl;
//   for(auto b: a.second)
//   {
//     std::cout<< "b "<< b.first<<std::endl;
//     for(auto c: b.second)
//     {
//       auto res=get_sign(c);
//       std::cout<<print_op(c)<< " signe "<<res[0]<< " "<< res[1]<< "  "<<  res[2]<<std::endl;
//     }
//   }
// }
//   //get_order_four_monomials(states, map_sec, Ly, Lx, true);
  basis_structure states_2;
  for(int i=0; i<2; i++)
  {for(int j=0; j<1; j++)
    //states[i].size(); j++)
    {
      states_2[i][j] = states[i][j];
    }
  }

for(auto a: states_2)
{
  std::cout<< "first sector "<<a.first<<std::endl;
  for(auto b: a.second)
  {
    std::cout<< "second sector "<<b.first<<std::endl;
    std::cout<<b.second.size()<<std::endl;
  }
}
std::cout<< "end sector analysis"<<std::endl;
  lattice.states_ = states_2;
  auto data = get_rdms(4, Lx);

  rdms_struct rdms={};
  //(data); //{}; // data);
  std::cout << rdms.size() << std::endl;

  Model::t M = new Model("sdo1");
  auto _M = finally([&]()
                    { M->dispose(); });
  auto basis = momentum_symmetry_solver_sos(lattice, M, rdms);
  for (auto a : lattice.TI_map_)
  {
    if (std::abs(a.second.second) < 1e-4)
    {
      std::cout << "XXXXX" << std::endl;
    }
  }
  //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
  //   for(auto k: basis.total_refs_)
  //   {

  // //	std::cout<<k.first<<std::endl;
  //   }
  double J = 1;
  double Delta = 1.;

  auto b = define_xxz2d_sos(lattice, J, Delta);

  basis.set_b(b);
  std::cout << "start fixing constraints " << std::endl;
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
  int i = 0;
  for (auto val : lattice.variable_map_)
  {
    if (val.first == "1" or val.first == "0")
    {
      std::cout << val.first << " " << -1. * (*(M->getConstraint(i)->dual()))[0] << std::endl;
      // np_vec(i, 0) = -1. * (*(M->getConstraint(i)->dual()))[0];

      i++;
    }
  }
 // thord order -0.702827317
  return;
}
void test_multiple_blocks_higher_order_2d_rdm_sos_with_linear_constraints()
{

//   // std::cout << "WARNING! Takes a lot of memory" << std::endl;
//   int Lx = 4;
//   int Ly = 4;
//   // op_vec op = {spin_op("y", {3, 2}, {Lx, Ly}), spin_op("z", {0, 2}, {Lx, Ly}), spin_op("x", {0, 3}, {Lx, Ly}), spin_op("x", {1, 0}, {Lx, Ly}), spin_op("z", {1, 3}, {Lx, Ly}), spin_op("y", {0, 3}, {Lx, Ly})};
//   // auto all_p = generate_all_permutations_xyz(op);
//   // std::cout << "start " << print_op(op) << std::endl;
//   // for (auto a : all_p)
//   // {
//   //   std::cout << print_op(a) << std::endl;
//   // }
//   auto lattice = SquareLattice(Ly, Lx, true, false, "xyz", "xyz");
//   auto map_sec = get_sector_map();
//   basis_structure states = get_states();
//   get_order_one_monomials(states, map_sec, Ly, Lx, true);
//   int r = 1;
//   get_order_two_monomials(states, map_sec, Ly, Lx, r, r, -r, -r, true);
//   // get_order_three_monomials(states, map_sec, Ly, Lx, true);

//   // get_order_four_monomials(states, map_sec, Ly, Lx, true);
//   basis_structure states_2;
//   states_2[0] = states[0];
//   states_2[1] = states[1];

//   lattice.states_ = states_2;
//   auto data = get_rdms(4, Lx);

//   rdms_struct rdms(data);
//   //(data); //{}; // data);
//   std::cout << rdms.size() << std::endl;

//   Model::t M = new Model("sdo1");
//   auto _M = finally([&]()
//                     { M->dispose(); });

//     // std::vector<SumOfOperators> linear_constratints;
//      std::cout<<lattice.states_[0].size()<<std::endl;
//      auto sz=spin_op("z", {0, 0}, states[0][0][0].offset_);
//     operator_and_coeff sz_with_coeff(1.0, {sz});
// //    operator_and_coeff szcons(1.0, {spin_op("z" {0, 0}, states[0][0][0].offset_)});
//   SumOfOperators sumofops;
//   sumofops.insert(sz_with_coeff);
//   std::vector<SumOfOperators> linear_constraint;
//   auto convs=convert_linear_constraints(lattice, linear_constraint);
//    auto basis = momentum_symmetry_solver_sos(lattice, M, rdms);
//    basis.set_linear_constraints_vec(convs);
//   // // for (auto a : lattice.TI_map_)
//   // // {
//   // //   if (std::abs(a.second.second) < 1e-4)
//   // //   {
//   // //     std::cout << "XXXXX" << std::endl;
//   // //   }
//   // // }
//   // //  {std::cout<< a.first << " -> "<<a.second.first<<std::endl;}
//   // //   for(auto k: basis.total_refs_)
//   // //   {

//   // // //	std::cout<<k.first<<std::endl;
//   // //   }
//   double J = 1;
//   double Delta = 1.;

//   auto b = define_xxz2d_sos(lattice, J, Delta);

//   basis.set_b(b);
//   std::cout << "start fixing constraints " << std::endl;
//   basis.fix_constrains();
//   // auto h = basis.get_costfunction();
//   // // //   std::cout<<C->toString()<<std::endl;

//   // std::cout << "starting solving SDP" << std::endl;
//   // basis.M_->objective(ObjectiveSense::Maximize, h);
//   // basis.M_->dataReport();
//   // M->setLogHandler([=](const std::string &msg)
//   //                  { std::cout << msg << std::flush; });
//   // basis.M_->solve();
//   // auto cons = M->getConstraint(0)->dual();
//   // std::cout << cons << std::endl;

//   // std::cout << "Solution : " << std::endl;
//   // std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;
//   // int i = 0;
//   // for (auto val : lattice.variable_map_)
//   // {
//   //   if (val.first == "1" or val.first == "0")
//   //   {
//   //     std::cout << val.first << " " << -1. * (*(M->getConstraint(i)->dual()))[0] << std::endl;
//   //     // np_vec(i, 0) = -1. * (*(M->getConstraint(i)->dual()))[0];

//   //     i++;
//   //   }
//   // }
//   // thord order -0.702827317
  return;
}
// // void test_y()
// // {
// //   int L=4;
// //  op_vec vec={spin_op("x", {2,2}, {L,L}),spin_op("x", {0,0}, {L,L})};

// // std::cout<< print_op(vec)<<std::endl;
// //  auto [fac, nf] =get_normal_form(vec);
// //  std::cout<< fac << " and "<< print_op(nf)<<std::endl;

// // 	    return;
// // }

// // void test_d8_symm()
// // {
// //   int L=6;
// // //  op_vec vec={spin_op("x", {0,0}, {L,L}),spin_op("x", {1,1}, {L,L}),spin_op("y", {2,1}, {L,L})};
// // // std::cout<< print_op(vec)<<std::endl;
// // // auto all_d8=generate_all_d8(vec,  L);
// // // std::cout<< "start "<<std::endl;
// // // for(auto b: all_d8)
// // // {
// // // std::cout<<print_op(b)<<std::endl;

// // // }
// // //s_[x,(1,0)]s_[y,(2,0)]s_[x,(4,3)]s_[y,(5,3)]
// // // s_[x,(1,0)]s_[y,(2,0)]s_[x,(4,3)]s_[y,(5,3)]
// // // s_[x,(0,5)]s_[y,(0,4)]s_[x,(3,2)]s_[y,(3,1)]
// // // s_[x,(5,0)]s_[y,(4,0)]s_[x,(2,3)]s_[y,(1,3)]
// // // s_[x,(0,1)]s_[y,(0,2)]s_[x,(3,4)]s_[y,(3,5)]
// // // s_[x,(5,0)]s_[y,(4,0)]s_[x,(2,3)]s_[y,(1,3)]
// // // s_[x,(0,1)]s_[y,(0,2)]s_[x,(3,4)]s_[y,(3,5)]
// // // s_[x,(1,0)]s_[y,(2,0)]s_[x,(4,3)]s_[y,(5,3)]
// // // s_[x,(0,5)]s_[y,(0,4)]s_[x,(3,2)]s_[y,(3,1)]
// // op_vec vec={spin_op("x", {1,0}, {L,L}),spin_op("y", {2,0}, {L,L}),spin_op("x", {4,3}, {L,L}),spin_op("y", {5,3}, {L,L})};
// // std::cout<<print_op(vec)<<std::endl;
// // auto all_d8=generate_all_d8(vec,  L);
// // for(auto b: all_d8)
// // {
// // std::cout<<print_op(b)<<std::endl;

// // }
// // 	    return;
// // }
// void test_primal_and_dual()
// {
//   int Lx = 4;
//   int Ly = 4;
//   auto lattice1 = SquareLattice(Ly, Lx, true, false, "xyz", "xyz");
//   auto lattice2 = SquareLattice(Ly, Lx, true, false, "xyz", "xyz");
//   auto map_sec = get_sector_map();
//   basis_structure states1 = get_states();
//   get_order_one_monomials(states1, map_sec, Ly, Lx, true);
//   get_order_two_monomials(states1, map_sec, Ly, Lx, 3, 3, -3, -3, true);

//   basis_structure states2 = get_states();
//   get_order_one_monomials(states2, map_sec, Ly, Lx, true);
//   get_order_two_monomials(states2, map_sec, Ly, Lx, 3, 3, -3, -3, true);

//   auto data = get_rdms(Lx, Lx);

//   rdms_struct rdms(data);
//   Model::t M1 = new Model("sdo1");
//   auto _M1 = finally([&]()
//                      { M1->dispose(); });

//   Model::t M2 = new Model("sdo1");
//   auto _M2 = finally([&]()
//                      { M2->dispose(); });

//   auto basis1 = momentum_symmetry_solver_sos(lattice1, states1, M1, rdms);
//   auto basis2 = momentum_symmetry_solver_dual(lattice2, states2, M2, rdms);
//   double J = 1;
//   double Delta = 1.;

//   auto b1 = define_xxz2d_sos(basis1.total_refs_, lattice1, J, Delta);

//   basis1.set_b(b1);

//   basis1.fix_constrains();
//   auto h1 = basis1.get_costfunction();

//   auto b2 = define_xxz2d_sos(basis2.total_refs_, lattice2, J, Delta);

//   basis2.set_b(b2);

//   basis2.fix_constrains();
//   auto h2 = basis2.get_costfunction();

//   basis1.M_->objective(ObjectiveSense::Maximize, h1);
//   basis1.M_->dataReport();
//   M1->setLogHandler([=](const std::string &msg)
//                     { std::cout << msg << std::flush; });
//   basis1.M_->solve();
//   // auto cons=M->getConstraint(0)->dual();
//   // std::cout<<cons<<std::endl;
//   basis2.M_->objective(ObjectiveSense::Minimize, h2);
//   basis2.M_->dataReport();
//   M2->setLogHandler([=](const std::string &msg)
//                     { std::cout << msg << std::flush; });
//   basis2.M_->solve();

//   std::cout << "Solution : " << std::endl;
//   std::cout << std::setprecision(9) << M1->primalObjValue() << std::endl;

//   std::cout << "Solution : " << std::endl;
//   std::cout << std::setprecision(9) << M2->primalObjValue() << std::endl;
//   auto y_fusion = Matrix::dense(basis2.total_refs_.size(), 1, basis2.y_->level());
//   int i = 0;
//   double max = 0.;
//   for (auto val : basis2.total_refs_)
//   {
//     if (val.first != "1" and val.first != "0")
//     {
//       // auto c_fusion=Matrix::dense(basis2.total_refs_.size(), 1, basis2.y_->level());
//       //  std::cout<< (*(basis2.y_->level()))[val.second]<< " and "<< -1.*(*(M1->getConstraint(i)->dual()))[0]<<std::endl;
//       max = std::max(max, std::abs((*(basis2.y_->level()))[val.second] - (-1. * (*(M1->getConstraint(i)->dual()))[0])));
//       i++;
//     }
//   }
//   // std::cout<< max << " and "<<std::abs(M1->primalObjValue()-M2->primalObjValue())<<std::endl;
//   assert(max < 1e-3);
//   assert(std::abs(M1->primalObjValue() - M2->primalObjValue()) < 1e-6);

//   return;
// }
void test_multiple_blocks_higher_order_J1J2_2d_rdm_sos()
{
  std::cout << "WARNING! Takes a lot of memory" << std::endl;
  int Lx = 4;
  int Ly = 4;

  auto lattice = SquareLattice(Ly, Lx, true, false, "xyz", "xyz");
  auto map_sec = get_sector_map();
  basis_structure states = get_states();
  get_order_one_monomials(states, map_sec, Ly, Lx, true);
  get_order_two_monomials(states, map_sec, Ly, Lx, 3, 3, -3, -3, true);
  get_order_three_monomials(states, map_sec, Ly, Lx, true);

  get_order_four_monomials(states, map_sec, Ly, Lx, true);
  auto data = get_rdms(Lx, Lx);
  basis_structure states2;
  states2[0] = states[0];
  states2[1] = states[1];
  states = states2;
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
  auto basis = momentum_symmetry_solver_sos(lattice, states, M, rdms);

  double J1 = 1;
  double J2 = 0.5;

  auto b = define_J1J22d_sos(basis.total_refs_, lattice, J1, J2);

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

// void test_multiple_blocks_higher_order_J1J2_1d_rdm_sos()
// {
//   std::cout << "WARNING! Takes a lot of memory" << std::endl;
//   int Lx = 16;
//   int Ly = 1;

//   auto lattice = SquareLattice(Ly, Lx, false, false, "xyz", "xyz");
//   auto map_sec = get_sector_map();
//   basis_structure states = get_states();
//   get_order_one_monomials(states, map_sec, Ly, Lx, true);
//   get_order_two_monomials(states, map_sec, Ly, Lx, 0, int(Lx / 2), 0, 0, true);
//   get_order_three_monomials_1d(states, map_sec, Ly, Lx, true);

//   get_order_four_monomials_1d(states, map_sec, Ly, Lx, true);
//   // for (auto a : states)
//   // {
//   //   std::cout << "sec " << a.first << std::endl;
//   //   for (auto b : a.second)
//   //   {
//   //     std::cout << print_op(b) << std::endl;
//   //   }
//   // }

//   auto data = get_rdms_1d(Lx, Lx);

//   rdms_struct rdms(data);

//   Model::t M = new Model("sdo1");
//   auto _M = finally([&]()
//                     { M->dispose(); });
//   auto basis = momentum_symmetry_solver_sos(lattice, states, M, rdms);

//   double J1 = 1;
//   double J2 = 0.2;

//   auto b = define_J1J2_1d_sos(basis.total_refs_, lattice, J1, J2);

//   basis.set_b(b);

//   basis.fix_constrains();
//   auto h = basis.get_costfunction();
//   // //   std::cout<<C->toString()<<std::endl;

//   std::cout << "starting solving SDP" << std::endl;
//   basis.M_->objective(ObjectiveSense::Maximize, h);
//   basis.M_->dataReport();
//   M->setLogHandler([=](const std::string &msg)
//                    { std::cout << msg << std::flush; });
//   basis.M_->solve();
//   auto cons = M->getConstraint(0)->dual();
//   std::cout << cons << std::endl;

//   std::cout << "Solution : " << std::endl;
//   std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;

//   // thord order -0.703680777

//   return;
// }

// void test_multiple_blocks_higher_order_xxz_1d_rdm_sos()
// {
//   std::cout << "WARNING! Takes a lot of memory" << std::endl;
//   int Lx = 16;
//   int Ly = 1;

//   auto lattice = SquareLattice(Ly, Lx, false, false, "xyz", "xyz");
//   auto map_sec = get_sector_map();
//   basis_structure states = get_states();
//   get_order_one_monomials(states, map_sec, Ly, Lx, true);
//   get_order_two_monomials(states, map_sec, Ly, Lx, 0, int(Lx / 2), 0, 0, true);
//   get_order_three_monomials_1d(states, map_sec, Ly, Lx, true);

//   get_order_four_monomials_1d(states, map_sec, Ly, Lx, true);

//   auto data = get_rdms_1d(Lx, Lx);

//   rdms_struct rdms(data);

//   Model::t M = new Model("sdo1");
//   auto _M = finally([&]()
//                     { M->dispose(); });
//   auto basis = momentum_symmetry_solver_sos(lattice, states, M, rdms);

//   double J = 1;
//   double Delta = 1.;

//   auto b = define_xxz_1d_sos(basis.total_refs_, lattice, J, Delta);

//   basis.set_b(b);

//   basis.fix_constrains();
//   auto h = basis.get_costfunction();
//   // //   std::cout<<C->toString()<<std::endl;

//   std::cout << "starting solving SDP" << std::endl;
//   basis.M_->objective(ObjectiveSense::Maximize, h);
//   basis.M_->dataReport();
//   M->setLogHandler([=](const std::string &msg)
//                    { std::cout << msg << std::flush; });
//   basis.M_->solve();
//   auto cons = M->getConstraint(0)->dual();
//   std::cout << cons << std::endl;

//   std::cout << "Solution : " << std::endl;
//   std::cout << std::setprecision(9) << M->primalObjValue() << std::endl;

//   return;
// }