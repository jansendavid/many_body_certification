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
#include "symmetries.hpp"
using namespace mosek::fusion;
using namespace monty;
using string_pair = std::pair<std::string, std::string>;

// file with spin hamiltonians for translation invariant system

// HAMILTONIANS WHEN USING TRANSLATION SYMMETRY
// NOTE: all are of the form vec{S_i}vec{S_j}

std::vector<double> define_correlation_function_sos(std::map<std::string, int> refs, SquareLattice lattice, std::pair<std::string, std::string> dirs, std::pair<std::vector<int>, std::vector<int>> pos)
{

      std::vector<double> vals(refs.size(), 0);
      auto offset_vector = lattice.get_offset_vec();
      op_vec v_p = {spin_op(dirs.first, pos.second, offset_vector), spin_op(dirs.first, pos.first, offset_vector)};
      auto [fac_p, nf_p] = get_normal_form(v_p);
      assert(abs(fac_p.imag()) < 1e-9);
      auto [key_p, coeff_map_p] = lattice.TI_map_.at(print_op(nf_p));
      assert(abs(coeff_map_p.imag()) < 1e-9);
      auto el_p = refs.at(key_p);

      vals[el_p] += fac_p.real() * coeff_map_p.real() / 4.;

      return vals;
}

std::vector<double> define_bilayer_correlation_sos(std::map<std::string, int> refs, SquareLattice lattice, std::pair<std::string, std::string> dirs, std::pair<std::vector<int>, std::vector<int>> pos)
{

      std::vector<double> vals(refs.size(), 0);
      auto offset_vector = lattice.get_offset_vec();
      op_vec v_p = {spin_op(dirs.first, pos.first, offset_vector), spin_op(dirs.second, pos.second, offset_vector)};

      auto [fac_p, nf_p] = get_normal_form(v_p);
      assert(abs(fac_p.imag()) < 1e-9);
      auto [key_p, coeff_map_p] = lattice.TI_map_.at(print_op(nf_p));
      assert(abs(coeff_map_p.imag()) < 1e-9);
      auto el_p = refs.at(key_p);

      vals[el_p] += fac_p.real() * coeff_map_p.real() / 4.;

      return vals;
}

std::vector<double> define_xxz2d_sos(std::map<std::string, int> refs, SquareLattice lattice, double J, double Delta)
{
      std::vector<string_pair> dirs{string_pair("x", "x"), string_pair("z", "z"), string_pair("y", "y")};
      std::vector<double> vals(refs.size(), 0);
      auto offset_vector = lattice.get_offset_vec();
      for (auto term : dirs)
      {

            op_vec v_p = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {1, 0}, offset_vector)};
            auto [fac_p, nf_p] = get_normal_form(v_p);

            auto [key_p, coeff_map_p] = lattice.TI_map_.at(print_op(nf_p));

            auto el_p = refs.at(key_p);

            double coeff = J;
            if (term == string_pair("z", "z"))
            {
                  coeff = Delta;
            }

            if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }
            {

                  vals[el_p] += coeff * fac_p.real() * coeff_map_p.real() / 4.;
            }

            op_vec v_t = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {0, 1}, offset_vector)};
            auto [fac_t, nf_t] = get_normal_form(v_t);
            auto [key_t, coeff_map_t] = lattice.TI_map_.at(print_op(nf_t));
            auto el_t = refs.at(key_t);
            coeff = J;
            if (term == string_pair("z", "z"))
            {
                  coeff = Delta;
            }

            if (std::abs(fac_t.imag()) > 1e-9 or std::abs(coeff_map_t.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }

            {

                  vals[el_t] += coeff * fac_t.real() * coeff_map_t.real() / 4.;
            }
      }

      return vals;
}
std::vector<double> define_J1J22d_sos(std::map<std::string, int> refs, SquareLattice lattice, double J1, double J2)
{
      std::vector<string_pair> dirs{string_pair("x", "x"), string_pair("z", "z"), string_pair("y", "y")};
      std::vector<double> vals(refs.size(), 0);
      auto offset_vector = lattice.get_offset_vec();
      // J1 nearest neighbour interaction
      for (auto term : dirs)
      {

            op_vec v_p = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {1, 0}, offset_vector)};
            auto [fac_p, nf_p] = get_normal_form(v_p);

            auto [key_p, coeff_map_p] = lattice.TI_map_.at(print_op(nf_p));

            auto el_p = refs.at(key_p);

            double coeff = J1;

            if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }
            {

                  vals[el_p] += coeff * fac_p.real() * coeff_map_p.real() / 4.;
            }

            op_vec v_t = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {0, 1}, offset_vector)};
            auto [fac_t, nf_t] = get_normal_form(v_t);
            auto [key_t, coeff_map_t] = lattice.TI_map_.at(print_op(nf_t));
            auto el_t = refs.at(key_t);

            if (std::abs(fac_t.imag()) > 1e-9 or std::abs(coeff_map_t.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }

            {

                  vals[el_t] += coeff * fac_t.real() * coeff_map_t.real() / 4.;
            }
      }
      // J2 next nearest neighbour interaction
      for (auto term : dirs)
      {

            op_vec v_p = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {1, 1}, offset_vector)};
            auto [fac_p, nf_p] = get_normal_form(v_p);

            auto [key_p, coeff_map_p] = lattice.TI_map_.at(print_op(nf_p));

            auto el_p = refs.at(key_p);

            double coeff = J2;

            if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }
            {

                  vals[el_p] += coeff * fac_p.real() * coeff_map_p.real() / 4.;
            }

            op_vec v_t = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {1, lattice.Lx_ - 1}, offset_vector)};
            auto [fac_t, nf_t] = get_normal_form(v_t);
            auto [key_t, coeff_map_t] = lattice.TI_map_.at(print_op(nf_t));
            auto el_t = refs.at(key_t);

            if (std::abs(fac_t.imag()) > 1e-9 or std::abs(coeff_map_t.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }

            {

                  vals[el_t] += coeff * fac_t.real() * coeff_map_t.real() / 4.;
            }
      }

      return vals;
}
std::vector<double> define_heisenberg_bilayer_sos(std::map<std::string, int> refs, SquareLattice lattice, double J_perpendicular, double J_parallel, double J_x)
{
      std::vector<string_pair> dirs{string_pair("x", "x"), string_pair("z", "z"), string_pair("y", "y")};
      std::vector<double> vals(refs.size(), 0);
      // J_parallel terms
      int layers = 2;
      auto offset_vector = lattice.get_offset_vec();
      for (auto term : dirs)
      { // iterate over layers
            for (int i = 0; i < layers; i++)
            {
                  {
                        op_vec v_p = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {i, 1, 0}, offset_vector)};
                        auto [fac_p, nf_p] = get_normal_form(v_p);

                        auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);
                        // lattice.find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);
                        assert(found);
                        auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                        auto el_p = refs.at(new_key_p);

                        if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                        {
                              std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                        }
                        {

                              vals[el_p] += J_parallel * fac_p.real() * coeff_map_p.real() / 4.;
                        }
                  }
                  {
                        op_vec v_p = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {i, 0, 1}, offset_vector)};
                        auto [fac_p, nf_p] = get_normal_form(v_p);

                        // auto [key_p, coeff_map_p] = find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);
                        auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                        assert(found);
                        auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                        auto el_p = refs.at(new_key_p);

                        if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                        {
                              std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                        }
                        {

                              vals[el_p] += J_parallel * fac_p.real() * coeff_map_p.real() / 4.;
                        }
                  }
            }
      }
      // J perpendicular
      for (auto term : dirs)
      {

            {
                  op_vec v_p = {spin_op(term.first, {0, 0, 0}, offset_vector), spin_op(term.second, {1, 0, 0}, offset_vector)};
                  auto [fac_p, nf_p] = get_normal_form(v_p);

                  auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                  assert(found);
                  auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                  //                auto [key_p, coeff_map_p] = find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);

                  auto el_p = refs.at(new_key_p);

                  if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                  {
                        std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                  }
                  {

                        vals[el_p] += J_perpendicular * fac_p.real() * coeff_map_p.real() / 4.;
                  }
            }
      }
      // J_x
      // iterate over layers
      for (int i = 0; i < layers; i++)
      {
            for (auto term : dirs)
            {

                  {
                        op_vec v_p = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {(i + 1) % layers, 0, 1}, offset_vector)};
                        auto [fac_p, nf_p] = get_normal_form(v_p);
                        auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                        assert(found);
                        auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                        // auto [key_p, coeff_map_p] = find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);

                        auto el_p = refs.at(new_key_p);

                        if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                        {
                              std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                        }
                        {

                              vals[el_p] += J_x * fac_p.real() * coeff_map_p.real() / 4.;
                        }
                  }

                  {
                        op_vec v_p = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {(i + 1) % layers, 1, 0}, offset_vector)};
                        auto [fac_p, nf_p] = get_normal_form(v_p);

                        // auto [key_p, coeff_map_p] = find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);
                        auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                        assert(found);
                        auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                        auto el_p = refs.at(new_key_p);

                        if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                        {
                              std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                        }
                        {

                              vals[el_p] += J_x * fac_p.real() * coeff_map_p.real() / 4.;
                        }
                  }
            }
      }

      return vals;
}

std::vector<double> define_heisenberg_bilayer_J2_sos(std::map<std::string, int> refs, SquareLattice lattice, double J_perpendicular, double J_parallel, double J_x, double J2)
{
      std::vector<string_pair> dirs{string_pair("x", "x"), string_pair("z", "z"), string_pair("y", "y")};
      std::vector<double> vals(refs.size(), 0);
      // J_parallel terms
      int layers = 2;
      auto offset_vector = lattice.get_offset_vec();
      for (auto term : dirs)
      {
            for (int i = 0; i < layers; i++)
            {
                  {
                        op_vec v_p = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {i, 1, 0}, offset_vector)};
                        auto [fac_p, nf_p] = get_normal_form(v_p);

                        // auto [key_p, coeff_map_p] = find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);
                        auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                        assert(found);
                        auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);

                        auto el_p = refs.at(new_key_p);

                        if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                        {
                              std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                        }
                        {

                              vals[el_p] += J_parallel * fac_p.real() * coeff_map_p.real() / 4.;
                        }
                  }
                  {
                        op_vec v_p = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {i, 0, 1}, offset_vector)};
                        auto [fac_p, nf_p] = get_normal_form(v_p);

                        // auto [key_p, coeff_map_p] = find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);
                        auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                        assert(found);
                        auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                        auto el_p = refs.at(new_key_p);

                        if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                        {
                              std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                        }
                        {

                              vals[el_p] += J_parallel * fac_p.real() * coeff_map_p.real() / 4.;
                        }
                  }
            }
      }
      // J perpendicular
      for (auto term : dirs)
      {

            {
                  op_vec v_p = {spin_op(term.first, {0, 0, 0}, offset_vector), spin_op(term.second, {1, 0, 0}, offset_vector)};
                  auto [fac_p, nf_p] = get_normal_form(v_p);

                  // auto [key_p, coeff_map_p] = find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);
                  auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                  assert(found);
                  auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                  auto el_p = refs.at(new_key_p);

                  if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                  {
                        std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                  }
                  {

                        vals[el_p] += J_perpendicular * fac_p.real() * coeff_map_p.real() / 4.;
                  }
            }
      }
      // J_x
      for (int i = 0; i < layers; i++)
      {
            for (auto term : dirs)
            {

                  {
                        op_vec v_p = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {(i + 1) % layers, 0, 1}, offset_vector)};
                        auto [fac_p, nf_p] = get_normal_form(v_p);

                        // auto [key_p, coeff_map_p] = find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);
                        auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                        assert(found);
                        auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                        auto el_p = refs.at(new_key_p);

                        if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                        {
                              std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                        }
                        {

                              vals[el_p] += J_x * fac_p.real() * coeff_map_p.real() / 4.;
                        }
                  }

                  {
                        op_vec v_p = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {(i + 1) % layers, 1, 0}, offset_vector)};
                        auto [fac_p, nf_p] = get_normal_form(v_p);

                        /// auto [key_p, coeff_map_p] = find_index_of_operator(nf_p, map, Lx, Ly, permuts, bilayer);
                        auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                        assert(found);
                        auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                        auto el_p = refs.at(new_key_p);

                        if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                        {
                              std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                        }
                        {

                              vals[el_p] += J_x * fac_p.real() * coeff_map_p.real() / 4.;
                        }
                  }
            }
      }
      // J2
      // J2 next nearest neighbour interaction
      for (int i = 0; i < layers; i++)
      {
            for (auto term : dirs)
            {

                  op_vec v_p = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {i, 1, 1}, offset_vector)};
                  auto [fac_p, nf_p] = get_normal_form(v_p);

                  // auto [key_p, coeff_map_p] = map.at(print_op(nf_p));
                  auto [found, key_p] = lattice.check_if_operator_exists(nf_p, lattice.TI_map_, false);

                  assert(found);
                  auto [new_key_p, coeff_map_p] = lattice.TI_map_.at(key_p);
                  auto el_p = refs.at(new_key_p);

                  double coeff = J2;

                  if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
                  {
                        std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                  }
                  {

                        vals[el_p] += coeff * fac_p.real() * coeff_map_p.real() / 4.;
                  }

                  op_vec v_t = {spin_op(term.first, {i, 0, 0}, offset_vector), spin_op(term.second, {i, 1, lattice.Lx_ - 1}, offset_vector)};
                  auto [fac_t, nf_t] = get_normal_form(v_t);
                  // auto [key_t, coeff_map_t] = map.at(print_op(nf_t));
                  auto [found_t, key_t] = lattice.check_if_operator_exists(nf_t, lattice.TI_map_, false);

                  assert(found_t);
                  auto [new_key_t, coeff_map_t] = lattice.TI_map_.at(key_t);
                  auto el_t = refs.at(new_key_t);

                  if (std::abs(fac_t.imag()) > 1e-9 or std::abs(coeff_map_t.imag()) > 1e-9)
                  {
                        std::cout << "error: Hamiltonian contains complex elements " << std::endl;
                  }

                  {

                        vals[el_t] += coeff * fac_t.real() * coeff_map_t.real() / 4.;
                  }
            }
      }
      return vals;
}
// 1d models

std::vector<double> define_xxz_1d_sos(std::map<std::string, int> refs, SquareLattice lattice, double J, double Delta)
{
      std::vector<string_pair> dirs{string_pair("x", "x"), string_pair("z", "z"), string_pair("y", "y")};
      std::vector<double> vals(refs.size(), 0);
      auto offset_vector = lattice.get_offset_vec();
      for (auto term : dirs)
      {

            op_vec v_t = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {0, 1}, offset_vector)};
            auto [fac_t, nf_t] = get_normal_form(v_t);
            auto [key_t, coeff_map_t] = lattice.TI_map_.at(print_op(nf_t));
            auto el_t = refs.at(key_t);
            double coeff = J;
            if (term == string_pair("z", "z"))
            {
                  coeff = Delta;
            }

            if (std::abs(fac_t.imag()) > 1e-9 or std::abs(coeff_map_t.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }

            {

                  vals[el_t] += coeff * fac_t.real() * coeff_map_t.real() / 4.;
            }
      }

      return vals;
}
std::vector<double> define_J1J2_1d_sos(std::map<std::string, int> refs, SquareLattice lattice, double J1, double J2)
{
      std::vector<string_pair> dirs{string_pair("x", "x"), string_pair("z", "z"), string_pair("y", "y")};
      std::vector<double> vals(refs.size(), 0);
      auto offset_vector = lattice.get_offset_vec();
      // J1 nearest neighbour interaction
      for (auto term : dirs)
      {

            op_vec v_t = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {0, 1}, offset_vector)};
            auto [fac_t, nf_t] = get_normal_form(v_t);
            auto [key_t, coeff_map_t] = lattice.TI_map_.at(print_op(nf_t));
            auto el_t = refs.at(key_t);
            double coeff = J1;
            if (std::abs(fac_t.imag()) > 1e-9 or std::abs(coeff_map_t.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }

            {

                  vals[el_t] += coeff * fac_t.real() * coeff_map_t.real() / 4.;
            }
      }
      // J2 next nearest neighbour interaction
      for (auto term : dirs)
      {

            op_vec v_p = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {0, 2}, offset_vector)};
            auto [fac_p, nf_p] = get_normal_form(v_p);

            auto [key_p, coeff_map_p] = lattice.TI_map_.at(print_op(nf_p));

            auto el_p = refs.at(key_p);

            double coeff = J2;

            if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }
            {

                  vals[el_p] += coeff * fac_p.real() * coeff_map_p.real() / 4.;
            }
      }

      return vals;
}

std::vector<double> define_TFI_1d_sos(std::map<std::string, int> refs, SquareLattice lattice, double J, double h)
{
      std::vector<string_pair> dirs{string_pair("z", "z")};
      std::vector<double> vals(refs.size(), 0);
      auto offset_vector = lattice.get_offset_vec();
      // J1 nearest neighbour interaction
      for (auto term : dirs)
      {

            op_vec v_t = {spin_op(term.first, {0, 0}, offset_vector), spin_op(term.second, {0, 1}, offset_vector)};
            auto [fac_t, nf_t] = get_normal_form(v_t);
            auto [key_t, coeff_map_t] = lattice.TI_map_.at(print_op(nf_t));
            auto el_t = refs.at(key_t);

            if (std::abs(fac_t.imag()) > 1e-9 or std::abs(coeff_map_t.imag()) > 1e-9)
            {
                  std::cout << "error: Hamiltonian contains complex elements " << std::endl;
            }

            {

                  vals[el_t] += J * fac_t.real() * coeff_map_t.real();
            }
      }
      // J2 next nearest neighbour interaction
      std::string term = "x";
      op_vec v_p = {spin_op(term, {0, 0}, offset_vector)};
      auto [fac_p, nf_p] = get_normal_form(v_p);

      auto [key_p, coeff_map_p] = lattice.TI_map_.at(print_op(nf_p));

      auto el_p = refs.at(key_p);

      if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
      {
            std::cout << "error: Hamiltonian contains complex elements " << std::endl;
      }
      {

            vals[el_p] += h * fac_p.real() * coeff_map_p.real();
      }

      return vals;
}

std::vector<double> define_magnetization_sos(std::map<std::string, int> refs, SquareLattice lattice, std::string term)
{

      auto offset_vector = lattice.get_offset_vec();
      std::vector<double> vals(refs.size(), 0);
      op_vec v_p = {spin_op(term, {0, 0}, offset_vector)};

      auto [fac_p, nf_p] = get_normal_form(v_p);

      auto [key_p, coeff_map_p] = lattice.TI_map_.at(print_op(nf_p));

      auto el_p = refs.at(key_p);

      if (std::abs(fac_p.imag()) > 1e-9 or std::abs(coeff_map_p.imag()) > 1e-9)
      {
            std::cout << "error: Hamiltonian contains complex elements " << std::endl;
      }
      {

            vals[el_p] += fac_p.real() * coeff_map_p.real();
      }

      return vals;
}
