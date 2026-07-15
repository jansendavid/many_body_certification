#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include "fusion.h"
#include "functions.hpp"
#include "lattices.hpp"
#include "operator_operations.hpp"
#include "sos/complex_momentum_dual_double.hpp"
#include "spin_hamiltonians_TIsym.hpp"
#include "symmetries.hpp"

using namespace mosek::fusion;
using namespace monty;

SumOfOperators define_periodic_xxz2d_operator(int Lx, int Ly, double J,
                                               double Delta)
{
  assert(Lx > 2);
  assert(Ly > 2);

  SumOfOperators hamiltonian;
  const std::vector<int> offset{Lx, Ly};
  const std::vector<std::string> directions{"x", "y", "z"};

  for (int x = 0; x < Lx; ++x)
  {
    for (int y = 0; y < Ly; ++y)
    {
      const std::vector<int> site{x, y};
      const std::vector<int> right{(x + 1) % Lx, y};
      const std::vector<int> up{x, (y + 1) % Ly};

      for (const auto &direction : directions)
      {
        const double coupling = direction == "z" ? Delta : J;
        const double coefficient = coupling / 4.;

        hamiltonian.insert(operator_and_coeff(
            coefficient,
            {spin_op(direction, site, offset),
             spin_op(direction, right, offset)}));
        hamiltonian.insert(operator_and_coeff(
            coefficient,
            {spin_op(direction, site, offset),
             spin_op(direction, up, offset)}));
      }
    }
  }

  return hamiltonian;
}

int main()
{
  const int Lx = 4;
  const int Ly = 4;
	const double J = 1.;
	const double Delta = 1.;

  auto map_sec = get_sector_map();
  basis_structure_with_sub states = get_states_with_sub();

  get_order_one_monomials_double(states, map_sec, Ly, Lx, true);
  const int r = 3;
  get_order_two_monomials_double(states, map_sec, Ly, Lx, r, r, -r, -r,
                                 true);
  get_order_three_monomials_double(states, map_sec, Ly, Lx, true);
  get_order_four_monomials_double(states, map_sec, Ly, Lx, true);

  basis_structure_with_sub retained_states;
  for (int sector = 0; sector < 2; ++sector)
    retained_states[sector] = states[sector];
  states = retained_states;

	auto hamiltonian = define_periodic_xxz2d_operator(Lx, Ly, J, Delta);
	const std::size_t expected_hamiltonian_terms =
		static_cast<std::size_t>(6 * Lx * Ly);
	assert(hamiltonian.get_terms().size() == expected_hamiltonian_terms);
	std::cout << "Full periodic Hamiltonian terms: "
			  << hamiltonian.get_terms().size() << std::endl;

	auto lattice = SquareLattice(states, Ly, Lx, true, false, "xyz", "xyz",
								 {}, hamiltonian, 3);
  auto data = get_rdms(Lx, 7);
  rdms_struct rdms(data);

  Model::t M = new Model("state_optimality");
  auto dispose_model = finally([&]() { M->dispose(); });
  const bool U1 = true;
  const bool maximize = true;
	const bool enable_state_optimality_conditions = true;
  auto basis = momentum_symmetry_solver_sos_double(lattice, M, rdms,
													maximize, U1,
													enable_state_optimality_conditions);

  auto objective_coefficients = define_xxz2d_sos(lattice, J, Delta);
  basis.set_b(objective_coefficients);
  basis.fix_constrains();
  auto objective = basis.get_costfunction();

  basis.M_->objective(ObjectiveSense::Maximize, objective);
  basis.M_->dataReport();
  M->setLogHandler([=](const std::string &message) {
    std::cout << message << std::flush;
  });
  basis.M_->solve();

  std::cout << "Solution: " << std::setprecision(9)
            << M->primalObjValue() << std::endl;
  return 0;
}
