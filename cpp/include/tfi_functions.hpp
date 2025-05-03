#pragma once
#include "fusion.h"
#include "spins.hpp"
#include <unordered_map>
#include <memory>
#include "symmetries.hpp"
using namespace mosek::fusion;
using namespace monty;
void add_state_tfi(basis_structure &states, op_vec op, int Ly, int Lx)
{
    // adds a state to a basis

    auto [fac, nf] = get_normal_form(op);
    bool print = false;

    auto sign = 0;
    bool found = false;
    if (nf.size() > 0)
    {
        auto all_ty = generate_all_translations_y(nf, Ly, 1);

        bool found = false;

        for (auto op_ty : all_ty)
        {
            auto all_t = generate_all_translations(op_ty, Lx);
            for (auto op_t : all_t)
            {
                if (print)
                {
                    std::cout << print_op(op_t) << std::endl;
                }
                auto it = std::find(states.at(0).begin(), states.at(0).end(), op_t);
                if (it != states.at(0).end())
                {

                    found = true;
                    break;
                }
            }
        }

        if (!found)
        {
            states.at(0).push_back(nf);
        }
    }

    return;
}
void get_order_one_monomials_tfi(basis_structure &states, int Ly, int Lx)
{

    std::vector<std::string> dirs = {"x", "y", "z"};
    assert(Ly == 1);
    for (auto s : dirs)
    {

        op_vec v0 = {spin_op(s, {0, 0}, {Ly, Lx})};

        auto sign = get_sec(v0);
        add_state_tfi(states, v0, Ly, Lx);
    }
    return;
}

void get_order_two_monomials_tfi(basis_structure &states, int Ly, int Lx, int r)
{
    std::vector<std::string> dirs = {"x", "y", "z"};
    assert(Ly == 1);
    for (int i = 1; i <= r; i++)
    {

        for (auto s1 : dirs)
        {

            for (auto s2 : dirs)
            {

                int ind = (Lx + i) % Lx;
                op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {0, ind}, {Ly, Lx})};

                auto [fac, vec] = get_normal_form(v0);

                add_state_tfi(states, vec, Ly, Lx);
            }
        }
    }
    return;
}

void get_order_three_monomials_tfi(basis_structure &states, int Ly, int Lx)
{

    std::vector<std::string> dirs = {"x", "y", "z"};
    std::vector<int> offsets_ = {Ly, Lx};
    assert(Ly == 1);
    for (auto s1 : dirs)
    {
        for (auto s2 : dirs)
        {
            for (auto s3 : dirs)
            {

                {
                    op_vec v0 = {spin_op(s1, {0, 0}, offsets_), spin_op(s2, {0, 1}, offsets_), spin_op(s3, {0, 2}, offsets_)};
                    auto [fac, vec] = get_normal_form(v0);

                    add_state_tfi(states, vec, Ly, Lx);
                }
            }
        }
    }
    return;
}
void get_order_four_monomials_tfi(basis_structure &states, int Ly, int Lx)
{
    std::vector<std::string> dirs = {"x", "y", "z"};
    std::vector<int> offsets_ = {Ly, Lx};
    assert(Ly == 1);
    for (auto s1 : dirs)
    {
        for (auto s2 : dirs)
        {
            for (auto s3 : dirs)
            {

                {
                    for (auto s4 : dirs)
                    {

                        op_vec v0 = {spin_op(s1, {0, 0}, offsets_), spin_op(s2, {0, 1}, offsets_), spin_op(s3, {0, 2}, offsets_), spin_op(s4, {0, 3}, {offsets_})};
                        auto [fac, vec] = get_normal_form(v0);

                        add_state_tfi(states, vec, Ly, Lx);
                    }
                }
            }
        }
    }
    return;
}

rdms_struct get_rdms_tfi(int Ly, int Lx, int dim, bool bilayer)
{
    rdms_struct data;

    std::vector<std::string> rdm_ops = {};

    std::string s = "x";
    if (dim >= 2)
    {

        for (int b = 0; b < int(Lx / 2); b++)
        {
            std::set<std::vector<int>> set_with_vector;
            set_with_vector.insert({0, 0});
            set_with_vector.insert({0, b});
            if (set_with_vector.size() == 2)
            {
                op_vec v0 = {spin_op(s, {0, 0}, {Ly, Lx}), spin_op(s, {0, b}, {Ly, Lx})};
                auto [coeff, nf] = get_normal_form(v0);

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

    std::cout << "number of two site denisty matrices " << data.size() << std::endl;
    if (dim >= 3)
    {
        {
            std::set<std::vector<int>> set_with_vector;
            set_with_vector.insert({0, 0});
            set_with_vector.insert({0, 1});
            set_with_vector.insert({0, 2});
            if (set_with_vector.size() == 3)
            {
                op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {Ly, Lx})};
                auto [coeff, nf] = get_normal_form(v0);

                {

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
            set_with_vector.insert({0, 0});
            set_with_vector.insert({0, 1});
            set_with_vector.insert({0, 2});
            set_with_vector.insert({0, 3});
            if (set_with_vector.size() == 4)
            {
                op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {Ly, Lx})};
                auto [coeff, nf] = get_normal_form(v0);

                {

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
            set_with_vector.insert({0, 0});
            set_with_vector.insert({0, 1});
            set_with_vector.insert({0, 2});
            set_with_vector.insert({0, 3});
            set_with_vector.insert({0, 4});

            if (set_with_vector.size() == 5)
            {
                op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 4), {Ly, Lx})};
                auto [coeff, nf] = get_normal_form(v0);

                {

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

        std::set<std::vector<int>> set_with_vector;
        set_with_vector.insert({0, 0});
        set_with_vector.insert({0, 1});
        set_with_vector.insert({0, 2});
        set_with_vector.insert({0, 3});
        set_with_vector.insert({0, 4});
        set_with_vector.insert({0, 5});

        if (set_with_vector.size() == 6)
        {
            op_vec v0 = {spin_op(s, *next(set_with_vector.begin(), 0), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 1), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 2), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 3), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 4), {Ly, Lx}), spin_op(s, *next(set_with_vector.begin(), 5), {Ly, Lx})};
            auto [coeff, nf] = get_normal_form(v0);

            {
                rdm_operator newstate;
                for (auto site_op : nf)
                {
                    newstate.op_.push_back(site_op.site_);
                }

                data.add_operator(newstate);
            }
        }
    }

    std::cout << "total number of density operators= " << data.size() << std::endl;
    return data;
}