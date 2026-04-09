
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include "sdp.hpp"

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
std::pair<std::string, std::complex<double>> get_key(op_vec spin_op)
{

    auto [fac, nf] = get_normal_form(spin_op);

    std::string key = print_op(nf);

   
        if (is_zero_signsym_xyz(nf))
        {
            key = "0";
        }
   
    return std::pair<std::string, std::complex<double>>(key, fac);
}
std::vector<double> get_hamiltonian(std::unordered_map<std::string, int> terms_mapping, std::unordered_map<std::string, std::string> &elements_mapping, int L)
{
    // defines the vector b_
    std::vector<double> expressions(terms_mapping.size(), 0);
std::vector<std::string> dir={"x","y","z"};
for(auto s1: dir)
{
    for (int i = 0; i < L; i++)
    {
   
        op_vec v0 = {spin_op(s1, {i}, {1}), spin_op(s1, {(i+1)%L}, {1})};

        auto [fac0, nf0] = get_normal_form(v0);
        auto key = elements_mapping.at(print_op(nf0));
        auto el = terms_mapping.at(key);
        assert(abs(fac0.imag()) < 1e-9);
        // std::cout << key << " fac x " << fac0 << std::endl;
        expressions[el] += fac0.real()/4;
       
    }
}
    return expressions;
}
int main()
{
    int L = 6;
    std::map<std::string,op_vec> map_string_op;
    std::vector<op_vec> block_0;
    std::vector<op_vec> block_1;
    int substructure=2;
    basis_structure_with_sub states=get_states_with_sub();
    auto map_sec = get_sector_map();
    std::vector<std::string> sym = {"x", "y", "z"};
    states[0][0].push_back({});
    for (int i = 0; i < L; i++)
    {
        for (auto s1 : sym)
        {
            op_vec v1 = {spin_op(s1, {i}, {1})};
            auto [fac1, nf1] = get_normal_form(v1);
            auto sign = get_sec(nf1);
            states[map_sec.at(sign)][nf1.size()%substructure].push_back(nf1);
          
        }
    }
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
        for (auto s1 : sym)
        {
            for (auto s2 : sym)
            {
            if(i!=j)
            {
            op_vec v1 = {spin_op(s1, {i}, {1}),spin_op(s2, {j}, {1})};
            auto [fac1, nf1] = get_normal_form(v1);
            auto sign = get_sec(nf1);
           states[map_sec.at(sign)][nf1.size()%substructure].push_back(nf1);
            }
       }
        }
    }
}
for (int i = 0; i < L; i++)
{
    // for (int j = 0; j < L; j++)
    // {
    for (auto s1 : sym)
    {
        for (auto s2 : sym)
        {
            for (auto s3 : sym)
            {
       // if(i!=j)
        // {
        op_vec v1 = {spin_op(s1, {i}, {1}),spin_op(s2, {(i+1)%L}, {1}),spin_op(s3, {(i+2)%L}, {1})};
        auto [fac1, nf1] = get_normal_form(v1);
        auto sign = get_sec(nf1);
        states[map_sec.at(sign)][nf1.size()%substructure].push_back(nf1);
        // }
   }
    }
}
    }
    basis_structure_with_sub states2;
    states2[0]=states[0];
    states2[1]=states[1];
     states=states2;
     for(auto se: states)
     {
        std::cout<< se.first<<std::endl;
        for(auto be: se.second)
        {
            std::cout<< be.first <<" and "<<be.second.size()<<std::endl;
        }
     }
    Model::t M = new Model("sdo1");
    auto _M = finally([&]()
                      { M->dispose(); });
    std::unordered_map<std::string, int> terms_mapping;
    std::unordered_map<int,std::unordered_map<std::string, matrix_organizer>> matrix_mapping;
    auto mat_terms = get_mat_terms_double(states);
    int index_old = 0;
    std::unordered_map<std::string, std::string> elements_mapping;
    for (auto symm_sec : mat_terms)
    {
        matrix_mapping.insert({symm_sec.first, {}});
    
        for (auto op : symm_sec.second)
        {

            bool found = false;
            {
                auto it = elements_mapping.find(print_op(op));
                if (it != elements_mapping.end())
                {
                    found = true;
                    elements_mapping.insert({print_op(op), it->second});
                }
            }

            if (not found)
            {
                elements_mapping.insert({print_op(op), print_op(op)});
                matrix_mapping[symm_sec.first].insert({print_op(op), matrix_organizer()});
                terms_mapping.insert({print_op(op), index_old});
                index_old++;
            }
        }
     
        }
        auto y = M->variable("T", terms_mapping.size());
        std::cout << "y shape " << terms_mapping.size() << std::endl;
    
        auto el = terms_mapping.at("1");
        M->constraint(y->index(el), Domain::equalsTo(1.0));
        M->constraint(y, Domain::lessThan(1.0));
        M->constraint(y, Domain::greaterThan(-1.0));
 
        for (auto & sec : states)
        {
            get_sdp_block_general_double(M, states[sec.first], matrix_mapping[sec.first], elements_mapping, terms_mapping, y);
 

        }
    // }
    auto obs = get_hamiltonian(terms_mapping, elements_mapping, L);
    auto obs_arr = monty::new_array_ptr<double>(obs);
    auto objective = Expr::dot(obs_arr, y);
    // auto objective = Expr::dot(real_part_arr, y);
    //M->objective(ObjectiveSense::Maximize, objective);
    M->objective(ObjectiveSense::Minimize, objective);
    M->dataReport();
    M->setLogHandler([=](const std::string &msg)
                     { std::cout << msg << std::flush; });
    M->solve();
    std::cout << "Solution : " << std::endl;
    std::cout << std::setprecision(9) << M->primalObjValue()/L << std::endl;
    return 0;
}