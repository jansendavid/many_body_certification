#pragma once
// functions to acton on operators, like building sum, commutator etc.
#include <string>
#include <iostream>
#include <iomanip>
#include <set>
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

    std::complex<double> get_coeff() const
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
class SumOfOperators{
    std::unordered_map<std::string, operator_and_coeff> terms_;

    public:
    
    std::unordered_map<std::string, operator_and_coeff>& get_terms()
    {
        return terms_;
    }
    const std::unordered_map<std::string, operator_and_coeff>& get_terms() const
    {
        return terms_;
    }
    void merge(std::unordered_map<std::string, operator_and_coeff>& to_merge)
    {
        terms_.merge(to_merge);
    }
    void erase(std::string key)
    {terms_.erase(key);
    return;}
     
     void insert(operator_and_coeff term)
     {
        auto key=print_op(term.get_op());
         auto it=terms_.find(key);
         if(it!=terms_.end())
         {
            it->second.add_coeff(term.get_coeff());
         }
         else{
            terms_.insert({key, term});
         }
    }
};

bool pauli_strings_anticommute(const op_vec &left, const op_vec &right)
{
    int different_overlaps = 0;
    for (const auto &left_factor : left)
    {
        for (const auto &right_factor : right)
        {
            if (left_factor.pos() == right_factor.pos() &&
                left_factor.get_dir() != right_factor.get_dir())
            {
                ++different_overlaps;
            }
        }
    }
    return (different_overlaps % 2) != 0;
}

SumOfOperators build_state_optimality_entry(
    const op_vec &v, const op_vec &w_dagger,
    const SumOfOperators &hamiltonian)
{
    SumOfOperators result;
    for (const auto &[label, hamiltonian_term] : hamiltonian.get_terms())
    {
        (void)label;
        const auto h = hamiltonian_term.get_op();
        const int multiplicity =
            static_cast<int>(pauli_strings_anticommute(v, h)) +
            static_cast<int>(pauli_strings_anticommute(w_dagger, h));
        if (multiplicity == 0)
            continue;

        op_vec product = v;
        product.insert(product.end(), h.begin(), h.end());
        product.insert(product.end(), w_dagger.begin(), w_dagger.end());
        auto [normal_phase, normal_form] = get_normal_form(std::move(product));
        const auto coefficient = static_cast<double>(multiplicity) *
                                 hamiltonian_term.get_coeff() * normal_phase;
        if (std::abs(coefficient) > 1e-12)
            result.insert(operator_and_coeff(coefficient, std::move(normal_form)));
    }

    std::vector<std::string> cancelled_terms;
    for (const auto &[label, term] : result.get_terms())
        if (std::abs(term.get_coeff()) <= 1e-12)
            cancelled_terms.push_back(label);
    for (const auto &label : cancelled_terms)
        result.erase(label);

    return result;
}

SumOfOperators build_state_optimality_entry_from_anticommuting_terms(
    const op_vec &v, const op_vec &w_dagger,
    const SumOfOperators &hamiltonian,
    const std::vector<std::string> &anticommuting_with_v,
    const std::vector<std::string> &anticommuting_with_w)
{
    std::unordered_map<std::string, int> multiplicities;
    for (const auto &label : anticommuting_with_v)
        ++multiplicities[label];
    for (const auto &label : anticommuting_with_w)
        ++multiplicities[label];

    SumOfOperators result;
    for (const auto &[label, multiplicity] : multiplicities)
    {
        const auto &hamiltonian_term = hamiltonian.get_terms().at(label);
        const auto h = hamiltonian_term.get_op();
        op_vec product = v;
        product.insert(product.end(), h.begin(), h.end());
        product.insert(product.end(), w_dagger.begin(), w_dagger.end());
        auto [normal_phase, normal_form] = get_normal_form(std::move(product));
        const auto coefficient = static_cast<double>(multiplicity) *
                                 hamiltonian_term.get_coeff() * normal_phase;
        if (std::abs(coefficient) > 1e-12)
            result.insert(operator_and_coeff(coefficient, std::move(normal_form)));
    }

    std::vector<std::string> cancelled_terms;
    for (const auto &[label, term] : result.get_terms())
        if (std::abs(term.get_coeff()) <= 1e-12)
            cancelled_terms.push_back(label);
    for (const auto &label : cancelled_terms)
        result.erase(label);
    return result;
}

template<typename LattceType>
std::vector<std::vector<double>> convert_linear_constraints(LattceType& lattice, std::vector<SumOfOperators>& elements)
{
    std::vector<std::vector<double>> results;
    for(int i=0; i<elements.size(); i++)
    {
        std::vector<double> real_part(lattice.variable_map_.size(), 0);
            std::vector<double> imag_part(lattice.variable_map_.size(), 0);
            bool include_imag=false;
            bool include_real=false;
        for(auto& term :elements[i].get_terms() )
        {
            
            auto [key, coeff_map] = lattice.TI_map_.at(key_dir_pos(term.second.get_op()));
            auto el = lattice.variable_map_.at(op_key_label(key));
            std::complex<double> total_coeff=coeff_map*term.second.get_coeff();
            if(std::abs(total_coeff.real())>1e-9)
            {
                real_part[el]+=total_coeff.real();
                include_real=true;
            }
            if(std::abs(total_coeff.imag())>1e-9)
            {

                imag_part[el]+=total_coeff.imag();
                include_imag=true;
            }
         

        }
        if(include_imag)
        {
            results.push_back(imag_part);
        }
        if(include_real)
        {
            results.push_back(real_part);
        }
    }
return results;
}
SumOfOperators multiply_two_ops(SumOfOperators &A, SumOfOperators &B)
{
    SumOfOperators res;

    for (auto a : A.get_terms())
    {

        for (auto b : B.get_terms())
        {

            auto v_x = a.second.get_op();
            auto vec_insert = b.second.get_op();
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
            auto fac_tot=a.second.get_coeff() * b.second.get_coeff() * fac;
            if(std::abs(fac_tot)>1e-9){
            operator_and_coeff op_with_coeff(fac_tot, vec);
            res.insert(op_with_coeff);
            }
        }
    
  
    }
    return res;
}
// std::vector<operator_and_coeff> compute_anticommutator(std::vector<operator_and_coeff> &A, std::vector<operator_and_coeff> &B)
// {
//     // return [A,B]
//     std::vector<operator_and_coeff> result;
//     auto A_B = multiply_two_ops(A, B);

//     auto B_A = multiply_two_ops(B, A);

//     for (auto b_a : B_A)
//     {
//         auto it = std::find_if(A_B.begin(), A_B.end(), [&b_a](operator_and_coeff op_temp)
//                                { return print_op(op_temp.get_op()) == print_op(b_a.get_op()); });
//         if (it != A_B.end())
//         {
//             it->add_coeff(b_a.get_coeff());
//         }
//         else
//         {
//             A_B.push_back(operator_and_coeff(b_a.get_coeff(), b_a.get_op()));
//         }
//     }

//     return A_B;
// }

SumOfOperators compute_commutator(SumOfOperators &A, SumOfOperators &B)
{
    // return [A,B]
    SumOfOperators result;
    auto A_B = multiply_two_ops(A, B);

    auto B_A = multiply_two_ops(B, A);
    for(auto term : B_A.get_terms())
    {
        operator_and_coeff op({-1.*term.second.get_coeff(), term.second.get_op()});
        A_B.insert(op);
    }
    std::set<std::string> to_delete_vec;
    for(auto term : A_B.get_terms())
    {
        if(std::abs(term.second.get_coeff())<1e-9)
        {
            to_delete_vec.insert(term.first);
    
        }

    }
    for(auto key_to_delete:to_delete_vec )
    {
        A_B.erase(key_to_delete);
    }


    return A_B;
}


std::ostream &operator<<(std::ostream &os, const operator_and_coeff &dt)
{
    os << "coeff: " << dt.coeff_ << " op: " << print_op(dt.op_);
    return os;
}
