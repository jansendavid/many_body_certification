#pragma once
// functions to acton on operators, like building sum, commutator etc.
#include <string>
#include <iostream>
#include <iomanip>
#include "spins.hpp"
class operator_and_coeff
{
    std::complex<double> coeff_;
    op_vec op_;

public:
    op_vec get_op()
    {
        return op_;
    }

    op_vec get_op() const
    {
        return op_;
    }
    void add_coeff(std::complex<double> to_add)
    {
        coeff_ += to_add;
        return;
    }
    void multiply_coeff(std::complex<double> to_multiply)
    {
        coeff_ *= to_multiply;
        return;
    }

    std::complex<double> get_coeff()
    {
        return coeff_;
    }

public:
    operator_and_coeff(std::complex<double> coeff, op_vec op) : coeff_(coeff), op_(op) {}

    friend inline bool operator==(const operator_and_coeff &lhs, const operator_and_coeff &rhs)
    {

        return (print_op(lhs.op_) == print_op(rhs.op_)) and (std::abs(lhs.coeff_ - rhs.coeff_) < 1e-9);
    }
    friend std::ostream &operator<<(std::ostream &os, const operator_and_coeff &dt);
    bool operator<(const operator_and_coeff &other) const
    {
        return print_op(op_) < print_op(other.get_op());
    }
};
std::vector<operator_and_coeff> multiply_two_ops(std::vector<operator_and_coeff> &A, std::vector<operator_and_coeff> &B)
{
    std::map<std::string, operator_and_coeff> res;

    for (auto a : A)
    {

        for (auto b : B)
        {

            auto v_x = a.get_op();
            auto vec_insert = b.get_op();
            std::complex<double> fac{1.};
            op_vec vec;
            if (print_op(v_x) == "1" and print_op(vec_insert) == "1")
            {
            }
            else
            {
                v_x.insert(v_x.end(), vec_insert.begin(), vec_insert.end());
                auto [fac_t, vec_t] = get_normal_form(v_x);
                fac = fac_t;
                vec = vec_t;
            }
            // if (print_op(vec) == "1")
            // {
            //     std::cout << print_op(a.get_op()) << " " << print_op(b.get_op()) << std::endl;
            // // }
            operator_and_coeff op_with_coeff(a.get_coeff() * b.get_coeff() * fac, vec);
            auto it = res.find(print_op(vec));

            if (it != res.end())
            {
                // if (print_op(vec) == "1")
                // {
                //     std::cout << "found " << a.get_coeff() << " " << b.get_coeff() << " " << fac << " " << op_with_coeff.get_coeff() << std::endl;
                // }
                auto m = op_with_coeff.get_coeff();
                it->second.add_coeff(op_with_coeff.get_coeff()); // op_with_coeff.get_coeff());
            }
            else
            {
                res.insert({print_op(vec), op_with_coeff});
            }
        }
    }
    std::vector<operator_and_coeff> res_vec;
    for (auto it = res.begin(); it != res.end(); ++it)
    {
        if (std::abs(it->second.get_coeff()) > 1e-9)
        {
            res_vec.push_back(it->second);
        }
    }
    return res_vec;
}
std::vector<operator_and_coeff> compute_anticommutator(std::vector<operator_and_coeff> &A, std::vector<operator_and_coeff> &B)
{
    // return [A,B]
    std::vector<operator_and_coeff> result;
    auto A_B = multiply_two_ops(A, B);

    auto B_A = multiply_two_ops(B, A);

    for (auto b_a : B_A)
    {
        auto it = std::find_if(A_B.begin(), A_B.end(), [&b_a](operator_and_coeff op_temp)
                               { return print_op(op_temp.get_op()) == print_op(b_a.get_op()); });
        if (it != A_B.end())
        {
            it->add_coeff(b_a.get_coeff());
        }
        else
        {
            A_B.push_back(operator_and_coeff(b_a.get_coeff(), b_a.get_op()));
        }
    }

    return A_B;
}

std::vector<operator_and_coeff> compute_commutator(std::vector<operator_and_coeff> &A, std::vector<operator_and_coeff> &B)
{
    // return [A,B]
    std::vector<operator_and_coeff> result;
    auto A_B = multiply_two_ops(A, B);

    auto B_A = multiply_two_ops(B, A);

    std::vector<operator_and_coeff> operators_to_delete;
    // for (auto)
    for (auto a_b : A_B)
    {
        auto it = find(B_A.begin(), B_A.end(), a_b);
        if (it != B_A.end())
        {

            operators_to_delete.push_back(operator_and_coeff(a_b.get_coeff(), a_b.get_op()));
        }

        //   // Check if the target value was found
    }

    // remove the values that cancel

    for (auto op_to_delete : operators_to_delete)
    {
        {
            auto it = find(A_B.begin(), A_B.end(), op_to_delete);
            A_B.erase(it);
        }
        {
            auto it = find(B_A.begin(), B_A.end(), op_to_delete);
            B_A.erase(it);
        }
    }
    for (auto b_a : B_A)
    {
        auto it = std::find_if(A_B.begin(), A_B.end(), [&b_a](operator_and_coeff op_temp)
                               { return print_op(op_temp.get_op()) == print_op(b_a.get_op()); });
        if (it != A_B.end())
        {
            it->add_coeff(-1. * b_a.get_coeff());
        }
        else
        {
            A_B.push_back(operator_and_coeff(-1. * b_a.get_coeff(), b_a.get_op()));
        }
    }

    return A_B;
}

struct SumOfOperators
{ // L=\sum_i l_i
    int N_;
    std::map<int, std::vector<std::pair<std::complex<double>, op_vec>>> op_;
    SumOfOperators(int N) : N_(N)
    { // Constructor
    }
};
std::ostream &operator<<(std::ostream &os, const operator_and_coeff &dt)
{
    os << "coeff: " << dt.coeff_ << " op: " << print_op(dt.op_);
    return os;
}