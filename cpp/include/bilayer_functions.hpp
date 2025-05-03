#pragma once
#include "fusion.h"
#include "spins.hpp"
#include <unordered_map>
#include <memory>
#include "symmetries.hpp"
using namespace mosek::fusion;
using namespace monty;

void get_order_one_monomials_bilayer(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int layers, int Ly, int Lx, bool use_symm)
{

  assert(Lx == Ly);
  std::vector<std::string> dirs = {"x", "y", "z"};

  for (int i = 0; i < layers; i++)
  {
    for (auto s : dirs)
    {

      op_vec v0 = {spin_op(s, {i, 0, 0}, {layers, Ly, Lx})};

      auto sign = get_sec(v0);

      if (use_symm)
      {
        add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
      }
      else
      {
        add_state(states, v0, map_sec);
      }
    }
  }
}

void get_order_two_monomials_bilayer(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int layers, int Ly, int Lx, int r, int start, bool use_symm)
{
  std::vector<std::string> dirs = {"x", "y", "z"};
  assert(Lx == Ly);

  for (int layer = 0; layer < layers; layer++)
  {
    for (int i = start; i <= r; i++)
    {

      for (int j = start; j <= r; j++)
      {
        for (auto s1 : dirs)
        {

          for (auto s2 : dirs)
          {

            if (i != 0 or j != 0)
            {
              int ind1 = (Ly + i) % Ly;
              int ind2 = (Lx + j) % Lx;
              op_vec v0 = {spin_op(s1, {0, 0, 0}, {layers, Ly, Lx}), spin_op(s2, {layer, ind1, ind2}, {layers, Ly, Lx})};

              auto [fac, vec] = get_normal_form(v0);

              if (use_symm)
              {
                add_state_with_symmetries(states, vec, map_sec, Ly, Lx);
              }
              else
              {
                add_state(states, vec, map_sec);
              }
            }
          }
        }
      }
    }
  }
}
void get_order_three_monomials_bilayer(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int layers, int Ly, int Lx, bool use_symm)
{

  std::vector<std::string> dirs = {"x", "y", "z"};
  std::vector<int> offsets_ = {layers, Ly, Lx};
  for (auto s1 : dirs)
  {
    for (auto s2 : dirs)
    {
      for (auto s3 : dirs)
      {

        {
          op_vec v0 = {spin_op(s1, {0, 0, 0}, offsets_), spin_op(s2, {0, 0, 1}, offsets_), spin_op(s3, {0, 0, 2}, offsets_)};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
        {
          op_vec v0 = {spin_op(s1, {0, 0, 0}, offsets_), spin_op(s2, {0, 0, 1}, offsets_), spin_op(s3, {1, 1, 1}, offsets_)};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
        {
          op_vec v0 = {spin_op(s1, {1, 0, 0}, offsets_), spin_op(s2, {0, 1, 0}, offsets_), spin_op(s3, {1, 2, 1}, offsets_)};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
      }
    }
  }
}
void get_order_four_monomials_bilayer(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int layers, int Ly, int Lx, bool use_symm)
{
  std::vector<std::string> dirs = {"x", "y", "z"};
  std::vector<int> offsets_ = {layers, Ly, Lx};
  for (auto s1 : dirs)
  {
    for (auto s2 : dirs)
    {
      for (auto s3 : dirs)
      {

        {
          for (auto s4 : dirs)
          {

            {
              op_vec v0 = {spin_op(s1, {0, 0, 0}, offsets_), spin_op(s2, {1, 0, 0}, offsets_), spin_op(s3, {0, 1, 1}, offsets_), spin_op(s4, {1, 1, 1}, {offsets_})};
              auto [fac, vec] = get_normal_form(v0);

              if (use_symm)
              {
                add_state_with_symmetries(states, vec, map_sec, Ly, Lx);
              }
              else
              {
                add_state(states, vec, map_sec);
              }
            }
          }
        }
      }
    }
  }
}

bool see_if_rdm_is_included(std::vector<std::string> included_rdms, op_vec monomial, int Ly, int Lx, bool bilayer)
{

  assert(Ly == Lx);
  auto all_ty = generate_all_translations_y(monomial, Ly, 1);

  int cnt{0};

  for (auto op_ty : all_ty)
  {
    auto all_t = generate_all_translations(op_ty, Lx);

    for (auto op_t : all_t)
    {
      cnt = count(included_rdms.begin(), included_rdms.end(), print_op(op_t));

      if (cnt > 0)
      {

        return true;
      }

      auto all_d8sym = generate_all_d8(op_t, Lx);

      for (auto d8s : all_d8sym)
      {
        cnt = count(included_rdms.begin(), included_rdms.end(), print_op(d8s));

        if (cnt > 0)
        {

          return true;
        }

        auto op_mirror = mirror(d8s);

        cnt = count(included_rdms.begin(), included_rdms.end(), print_op(op_mirror));

        if (cnt > 0)
        {

          return true;
        }
        if (bilayer)
        {
          auto op_flip_layer = flip_layer(d8s);
          cnt = count(included_rdms.begin(), included_rdms.end(), print_op(op_flip_layer));

          if (cnt > 0)
          {

            return true;
          }
          op_flip_layer = flip_layer(op_mirror);
          cnt = count(included_rdms.begin(), included_rdms.end(), print_op(op_flip_layer));

          if (cnt > 0)
          {

            return true;
          }
        }
      }
    }
  }
  return false;
}
rdms_struct get_rdms_bilayer(int Ly, int Lx, int dim, bool bilayer)
{
  rdms_struct data;
  int layers = 2;
  std::vector<std::string> rdm_ops = {};

  std::string s = "x";
  if (dim >= 2)
  {
    for (int n = 0; n < layers; n++)
    {

      for (int m = 0; m < layers; m++)
      {
        int i = 0;
        for (int j = 0; j < Lx; j++)
        {
          for (int a = 0; a < Lx; a++)
          {
            for (int b = 0; b < Lx; b++)
            {
              std::set<std::vector<int>> set_with_vector;
              set_with_vector.insert({n, i, j});
              set_with_vector.insert({m, a, b});
              if (set_with_vector.size() == 2)
              {
                op_vec v0 = {spin_op(s, {n, i, j}, {layers, Ly, Lx}), spin_op(s, {m, a, b}, {layers, Ly, Lx})};
                auto [coeff, nf] = get_normal_form(v0);
                auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
                if (!found)
                {
                  rdm_ops.push_back(print_op(nf));
                  rdm_operator newstate;
                  for (auto site_op : nf)
                  {
                    newstate.op_.push_back(site_op.site_);
                  }
                  data.add_operator(newstate);
                }
              }
            }
          }
        }
      }
    }
  }
  std::cout << "number of two site denisty matrices " << data.size() << std::endl;
  if (dim >= 3)
  {
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({0, 0, 1});
      set_with_vector.insert({0, 0, 1});
      if (set_with_vector.size() == 3)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({0, 1, 0});
      set_with_vector.insert({0, 1, 1});
      if (set_with_vector.size() == 3)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({1, 0, 0});
      set_with_vector.insert({0, 1, 0});
      if (set_with_vector.size() == 3)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({1, 1, 1});
      set_with_vector.insert({0, 1, 0});
      if (set_with_vector.size() == 3)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({1, 0, 1});
      set_with_vector.insert({1, 2, 0});
      if (set_with_vector.size() == 3)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
  }
  std::cout << "after three " << data.size() << std::endl;
  if (dim >= 4)
  {
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({0, 1, 0});
      set_with_vector.insert({0, 2, 0});
      set_with_vector.insert({0, 3, 0});
      if (set_with_vector.size() == 4)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({1, 0, 0});
      set_with_vector.insert({0, 1, 1});
      set_with_vector.insert({1, 1, 1});
      if (set_with_vector.size() == 4)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({0, 1, 1});
      set_with_vector.insert({0, 2, 0});
      set_with_vector.insert({1, 3, 1});
      if (set_with_vector.size() == 4)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
  }
  std::cout << "after four " << data.size() << std::endl;

  if (dim >= 5)
  {
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({1, 0, 0});
      set_with_vector.insert({0, 1, 1});
      set_with_vector.insert({1, 1, 1});
      set_with_vector.insert({0, 2, 2});

      if (set_with_vector.size() == 5)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 4), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({1, 0, 0});
      set_with_vector.insert({0, 1, 1});
      set_with_vector.insert({1, 1, 1});
      set_with_vector.insert({1, 2, 2});
      if (set_with_vector.size() == 5)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 4), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({1, 0, 0});
      set_with_vector.insert({0, 1, 1});
      set_with_vector.insert({1, 1, 1});
      set_with_vector.insert({1, 3, 2});
      if (set_with_vector.size() == 4)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 4), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
  }

  std::cout << "after five " << data.size() << std::endl;

  if (dim >= 6)
  {

    if (Lx > 5)
    {
      std::set<std::vector<int>> set_with_vector;
      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({0, 1, 0});
      set_with_vector.insert({0, 2, 0});
      set_with_vector.insert({0, 3, 0});
      set_with_vector.insert({0, 4, 0});
      set_with_vector.insert({0, 5, 0});

      if (set_with_vector.size() == 6)
      {
        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 4), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 5), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }

          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;

      set_with_vector.insert({0, 0, 0});
      set_with_vector.insert({0, 0, 1});
      set_with_vector.insert({1, 1, 1});
      set_with_vector.insert({0, 2, 2});
      set_with_vector.insert({0, 2, 3});
      set_with_vector.insert({1, 3, 3});
      if (set_with_vector.size() == 6)
      {

        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 4), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 5), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
    {
      std::set<std::vector<int>> set_with_vector;

      set_with_vector.insert({1, 0, 0});
      set_with_vector.insert({0, 1, 0});
      set_with_vector.insert({1, 2, 1});
      set_with_vector.insert({1, 1, 1});
      set_with_vector.insert({0, 2, 1});
      set_with_vector.insert({1, 3, 2});
      if (set_with_vector.size() == 6)
      {

        op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 4), {layers, Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 5), {layers, Ly, Lx})};
        auto [coeff, nf] = get_normal_form(v0);
        auto found = see_if_rdm_is_included(rdm_ops, nf, Ly, Lx, bilayer);
        if (!found)
        {
          rdm_ops.push_back(print_op(nf));
          rdm_operator newstate;
          for (auto site_op : nf)
          {
            newstate.op_.push_back(site_op.site_);
          }
          data.add_operator(newstate);
        }
      }
    }
  }

  std::cout << "total number of density operators= " << data.size() << std::endl;
  return data;
}