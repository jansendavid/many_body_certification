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
int main()
{ int Lx = 4;
    int Ly = 4;

   
 
 
   auto map_sec = get_sector_map();
   basis_structure states = get_states();
 
//  std::cout<<states.size()<<std::endl;
//  for(auto b: states)
//  {
//    std::cout<<b.second.size()<<std::endl;
//  }
// adding one

    get_order_one_monomials(states, map_sec, Ly, Lx, true);
   int r = 3;
  get_order_two_monomials(states, map_sec, Ly, Lx, r, r, -r, -r, true);
 get_order_three_monomials(states, map_sec, Ly, Lx, true);
 
  get_order_four_monomials(states, map_sec, Ly, Lx, true);
std::cout<<states.size()<<std::endl;
for(auto b: states)
{

   std::cout<< "b first "<<b.first<<std::endl;
  std::cout<<b.second.size()<<std::endl;
 
}
   basis_structure states_2;
   for(int i=0; i<2; i++)
   {
     {
       states_2[i]= states[i];
     }
   }
   states=states_2;

 std::cout<< "end sector analysis"<<std::endl;
  //lattice.states_ = states_2;
   auto lattice = SquareLattice(states,Ly, Lx, true, false, "xyz", "xyz");
   auto data = get_rdms(Lx ,7);
   std::cout<< "data size "<<data.size()<<std::endl;
   rdms_struct rdms(data); //{}; // data);
  
   std::cout << "rdm size "<<rdms.size() << std::endl;
 
   Model::t M = new Model("sdo1");
   auto _M = finally([&]()
                     { M->dispose(); });
                     bool U1=true;
                     bool maximize=true;
   auto basis = momentum_symmetry_solver_sos(lattice, M, rdms, maximize, U1);


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
   auto con = M->getConstraint(0);

    auto dual = con->dual();
   for (auto val : lattice.variable_map_)
   {
    if (val.first == "1" || val.first == "0")
     {
      std::cout <<val.first << " " << -1*(*dual)[i] << std::endl;
       //std::cout << val.first << " " << -1. * (*(M->getConstraint(i)->dual()))[0] << std::endl;
       // np_vec(i, 0) = -1. * (*(M->getConstraint(i)->dual()))[0];
 
       i++;
     }
   }
    // thord order -0.702827317
    return 0;
}