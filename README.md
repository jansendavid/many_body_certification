# many_body_certification
This library implements semidefinite programs to compute lower bounds for ground-state energies and upper and lower bounds for observables for quantum spin systems.
The code uses different symmetries, including translation invariance, sign symmetry of the Hamiltonian, D8 symmetry and more. The supported models are the XXZ chain in one and two dimension, XY chain in one and two dimensions, the transversefield ising model in one and two dimensions, the frustrated bilayer Heisenber model (with and without next nearest neighbor interaction). Other Hamiltonians should be easy to implement. This code was used to prepare the manuscript [Mapping phase diagrams of quantum spin systems through semidefinite-programming relaxations](https://arxiv.org/abs/2507.03137).
## Installation
This code is written in C++ depends on [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [xtensor](https://xtensor.readthedocs.io/en/latest/) and [Mosek fusion](https://docs.mosek.com/11.0/cxxfusion/index.html) (license required, but students and researchers should be able to get an academic one). Then, just clone this directory and set the paths in the Make file.
## Running code
First, innclude the necissary libraries:
```C++
#include "sos/complex_momentum_dual.hpp"
#include "spin_hamiltonians_TIsym.hpp"
#include "functions.hpp"
#include "reduced_dms.hpp"
using namespace mosek::fusion;
using namespace monty;
```
Formulate the problem the 2d J1-J1 model:
```C++
  # system size
  int Lx = 4;
  int Ly = 4;
  # define the lattice with symmetries to include
  auto lattice = SquareLattice(Ly, Lx, true, false, "xyz", "xyz");
  # define a set of monomials of different order, using sign symmetry sectors
  auto map_sec = get_sector_map();
  basis_structure states = get_states();
  get_order_one_monomials(states, map_sec, Ly, Lx, true);
  get_order_two_monomials(states, map_sec, Ly, Lx, 3, 3, -3, -3, true);
  get_order_three_monomials(states, map_sec, Ly, Lx, true);

  get_order_four_monomials(states, map_sec, Ly, Lx, true);
  auto data = get_rdms(Lx, Lx);
  // for this problem, only the first two sectors are needed
  basis_structure states2;
  states2[0] = states[0];
  states2[1] = states[1];
  states = states2;
  // include reduced density matrices for semidefinite positivity
  rdms_struct rdms(data);
```
Formulate problem:
```C++
  Model::t M = new Model("sdo1");
  auto _M = finally([&]()
                    { M->dispose(); });
  // prepare SDP
  auto basis = momentum_symmetry_solver_sos(lattice, states, M, rdms);

  double J1 = 1;
  double J2 = 0.5;
  // define Hamiltonian
  auto b = define_J1J22d_sos(basis.total_refs_, lattice, J1, J2);

  basis.set_b(b);
  
  basis.fix_constrains();
```
Solve SDP:
```C++
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
```
Many more problem examples can be found in cpp/include/test_functions.hpp
