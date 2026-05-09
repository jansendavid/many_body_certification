#pragma once

#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include "spins.hpp"
#include <unordered_map>
#include <Eigen/Dense>
#include "util.hpp"
#include "spins.hpp"
#include "symmetries.hpp"

using namespace mosek::fusion;
using namespace monty;
using string_pair = std::pair<std::string, std::string>;

template <typename T>
 std::unordered_map<int,std::vector<op_vec>> get_mat_terms(T &states)
{
    std::unordered_map<int,std::vector<op_vec>> mat_terms;
    for (auto sec : states)
    {
        mat_terms.insert({sec.first, {}});
     
        //     mat_terms[sec.first].insert({subsec.first, {}});
        auto operator_sector = states[sec.first];

        for (int i = 0; i < operator_sector.size(); i++)
        {
         
            auto op_dagger = dagger_operator(operator_sector[i]);
            for (int j = 0; j < operator_sector.size(); j++)
            {

                auto v_x = op_dagger;

                v_x.insert(v_x.end(), operator_sector[j].begin(), operator_sector[j].end());
                auto [coeff_, nf] = sdp_get_form(v_x);
                mat_terms[sec.first].push_back(nf);
            }
        
    }
}
    return mat_terms;
}

template <typename T>
void make_block(std::vector<T> &v_tot, std::unordered_map<std::string, matrix_organizer> &matrix_mapping, std::unordered_map<std::string, std::string> &elements_mapping, int shift, int unit)
{
    // make induvidual matrix block
    // generating all matrix elements
    for (int i = 0; i < v_tot.size(); i++)
    {

        for (int j = 0; j < v_tot.size(); j++)
        {

            auto v_x = v_tot[i];
            auto op_dagger = dagger_operator(v_x);

            auto [coeff_x, nf_x] = sdp_get_form(op_dagger);

            auto [coeff_y, nf_y] = sdp_get_form(v_tot[j]);

            nf_x.insert(nf_x.end(), nf_y.begin(), nf_y.end());

            auto [coeff_, nf] = get_normal_form(nf_x);
            auto coeff_tot = coeff_x * coeff_y * coeff_;
            auto found_string = elements_mapping.at(print_op(nf));

            if (std::abs(coeff_tot.real()) > 1e-9)
            {
                // using 1/4 or 1/2 prefacor?, since H=X_1+X2, but I add A+A^T

                matrix_mapping[found_string]
                    .add_values({i + unit, j + unit}, 1. / 2 * coeff_tot.real());
                matrix_mapping[found_string].add_values({i + unit + shift, j + unit + shift}, 1. / 2 * coeff_tot.real());
            }
            if (std::abs(coeff_tot.imag()) > 1e-9)
            {

                matrix_mapping[found_string].add_values({i + unit, shift + j + unit}, -1. / 2 * coeff_tot.imag());
                matrix_mapping[found_string].add_values({shift + i + unit, j + unit}, 1. / 2 * coeff_tot.imag());

                // matrix_mapping[found_string].add_values({i + unit, shift + j + unit}, -1. / 2 * coeff.imag());
                // matrix_mapping[print_op(op)].add_values({shift + j + 1, i + 1}, -1. / 2 * coeff.imag());
            }
        }
    }

    return;
}
template <typename T>
void get_sdp_block_zero(Model::t &M, std::vector<T> &v_tot, std::unordered_map<std::string, matrix_organizer> &matrix_mapping, std::unordered_map<std::string, std::string> &elements_mapping,std::unordered_map<std::string, int> &terms_mapping, Variable::t &y)
{

    int dim = 2 * (v_tot.size() + 1);

    int shift = v_tot.size() + 1;
    int unit = 1;
    matrix_mapping["1"].add_values({0, 0}, 1. / 2);
    matrix_mapping["1"].add_values({shift, shift}, 1. / 2);
    make_block(v_tot, matrix_mapping, elements_mapping, shift, unit);

    for (int i = 0; i < v_tot.size(); i++)
    {

        auto [coeff, op] = sdp_get_form(v_tot[i]);

        auto found_string = elements_mapping.at(print_op(op));

        auto el = terms_mapping.at(found_string);

        if (std::abs(coeff.real()) > 1e-9)
        {
            // std::cout << matrix_mapping[print_op(op)].make_matrix(dim, dim)->toString() << std::endl;
            matrix_mapping[found_string].add_values({0, i + 1}, 1. / 2 * coeff.real());
            matrix_mapping[found_string].add_values({i + 1, 0}, 1. / 2 * coeff.real());

            matrix_mapping[found_string].add_values({shift, i + 1 + shift}, 1. / 2 * coeff.real());
            matrix_mapping[found_string].add_values({i + 1 + shift, shift}, 1. / 2 * coeff.real());

            // std::cout << matrix_mapping[print_op(op)].make_matrix(dim, dim)->toString() << std::endl;
        }

        if (std::abs(coeff.imag()) > 1e-9)
        {
            // fccc std::cout<<i + 1<< ","<< shift<<std::endl
            matrix_mapping[found_string].add_values({i + 1, shift}, 1. / 2 * coeff.imag());
            matrix_mapping[found_string].add_values({0, i + 1 + shift}, -1. / 2 * coeff.imag());
            matrix_mapping[found_string].add_values({i + 1 + shift, 0}, -1. / 2 * coeff.imag());
            matrix_mapping[found_string].add_values({shift, i + 1}, 1. / 2 * coeff.imag());
            // //   std::cout << "error imaginare non zero" << std::endl;
        }
    }
    // // std::cout << "finished fisrt loop" << std::endl;
    // std::cout << "start, size= " << dim << std::endl;
    auto it = matrix_mapping.begin();
    // std::cout << it->first << std::endl;
    // matrix_mapping[it->first].print();
    auto expression = Expr::mul(y->index(terms_mapping.at(it->first)), matrix_mapping[it->first].make_matrix(dim, dim));

    it++;
    int l = 1;
    while (it != matrix_mapping.end())
    {

        expression = Expr::add(expression, Expr::mul(y->index(terms_mapping.at(it->first)), matrix_mapping[it->first].make_matrix(dim, dim)));
        it++;
        l++;
    }

    M->constraint(expression, Domain::inPSDCone());
    return;
}
template <typename T>
void get_sdp_block_general(Model::t &M, std::vector<T> &v_tot, std::unordered_map<std::string, matrix_organizer> &matrix_mapping, std::unordered_map<std::string, std::string> &elements_mapping, std::unordered_map<std::string, int> &terms_mapping, Variable::t &y)
{

    // correct this should be but must remove in make_block
    int dim = 2 * (v_tot.size());

    int unit = 0;

    int shift = v_tot.size();

    make_block(v_tot, matrix_mapping, elements_mapping, shift, unit);

    auto it = matrix_mapping.begin();
    auto expression = Expr::mul(y->index(terms_mapping.at(it->first)), matrix_mapping[it->first].make_matrix(dim, dim));

    it++;
    int l = 1;
    while (it != matrix_mapping.end())
    {

        expression = Expr::add(expression, Expr::mul(y->index(terms_mapping.at(it->first)), matrix_mapping[it->first].make_matrix(dim, dim)));
        it++;
        l++;
    }

    M->constraint(expression, Domain::inPSDCone());
    return;
}

/// using double symmetry
template <typename T>
void run_sdp_loop_double(std::vector<T> &operator_sector_1, std::vector<T> &operator_sector_2, std::unordered_map<std::string, matrix_organizer> &matrix_mapping, std::unordered_map<std::string, std::string> &elements_mapping, std::pair<int,int> shift, std::complex<double> prefac)
{
   std::cout<<"start "<<shift.first << " and "<<shift.second<<std::endl;
    int i=0;
  for (auto it1=operator_sector_1.begin(); it1!=operator_sector_1.end(); ++it1)
    { auto op1=*it1;
        auto op_dagger = dagger_operator(op1);
int j=0;
        for (auto it2=operator_sector_2.begin(); it2!=operator_sector_2.end(); ++it2)
        {
           
            auto op2=*it2;
            auto vx=op_dagger;
            //std::cout<< " op1 "<<print_op(op1)<< " op2 "<<print_op(op2)<<std::endl;
            vx.insert(vx.end(), op2.begin(), op2.end());
            //std::cout<<print_op(vx)<<std::endl;
            auto [coeff_, nf] = get_normal_form(vx);
            //std::cout<<print_op(nf)<< " with suze "<< nf.size()<< " and factor "<<coeff_<<std::endl;
          //  auto coeff_tot = coeff_x * coeff_y * coeff_;
          auto coeff_tot_=coeff_*prefac;
          if(nf.size()%2==0)
          {
         
            auto found_string = elements_mapping.at(print_op(nf));
     
            //         
        
                    if (std::abs(coeff_tot_.real()) > 1e-9)
                    {
                        // using 1/4 or 1/2 prefacor?, since H=X_1+X2, but I add A+A^T
       // std::cout<< "indices "<<shift.first+i<< " "<<shift.second+j<<std::endl;
                        matrix_mapping[found_string]
                            .add_values({shift.first+i, shift.second+j}, 1. / 2 * coeff_tot_.real());

                    }
                    if (std::abs(coeff_tot_.imag()) > 1e-9)
                    {
                        std::cout<< "error"<<std::endl;
        }
    }
   
     j++;}
     i++;
    }

}
template <typename T>
void make_block_double(std::map<int,std::vector<T>> &v_tot, std::unordered_map<std::string, matrix_organizer> &matrix_mapping, std::unordered_map<std::string, std::string> &elements_mapping)
{
    std::pair<int, int> shift={0,0};
    std::complex<double> prefac(1.,0);
    run_sdp_loop_double(v_tot[0], v_tot[0], matrix_mapping, elements_mapping, shift, prefac);

    shift={v_tot[0].size(),v_tot[0].size()};
    run_sdp_loop_double(v_tot[1], v_tot[1], matrix_mapping, elements_mapping, shift, prefac);
    prefac={0.,-1.};
    shift={v_tot[0].size(),0};
    run_sdp_loop_double(v_tot[1], v_tot[0], matrix_mapping, elements_mapping, shift, prefac);
    prefac={0.,-1.};
    shift={0,v_tot[0].size()};
    run_sdp_loop_double(v_tot[0], v_tot[1], matrix_mapping, elements_mapping, shift, prefac);


    return;
}
template <typename T>
void get_sdp_block_general_double(Model::t &M, std::map<int,std::vector<T>> &v_tot, std::unordered_map<std::string, matrix_organizer> &matrix_mapping, std::unordered_map<std::string, std::string> &elements_mapping, std::unordered_map<std::string, int> &terms_mapping, Variable::t &y)
{

    // correct this should be but must remove in make_block
    int dim = v_tot[0].size()+v_tot[1].size();
    std::cout<< "dim "<<dim <<std::endl;
    // int unit = 0;

    // int shift = v_tot.size();

    make_block_double(v_tot, matrix_mapping, elements_mapping);

    auto it = matrix_mapping.begin();
    auto expression = Expr::mul(y->index(terms_mapping.at(it->first)), matrix_mapping[it->first].make_matrix(dim, dim));

    it++;
    int l = 1;
    while (it != matrix_mapping.end())
    {

        expression = Expr::add(expression, Expr::mul(y->index(terms_mapping.at(it->first)), matrix_mapping[it->first].make_matrix(dim, dim)));
        it++;
        l++;
    }

    M->constraint(expression, Domain::inPSDCone());
    return;
}

template <typename T>
void run_loop(std::vector<T>& mat_terms,std::vector<T>& operator_sector_1, std::vector<T>& operator_sector_2, int sec)
{
    for(auto it1=operator_sector_1.begin(); it1!=operator_sector_1.end(); ++it1)
    {
        auto op1=*it1;
        auto op_dagger = dagger_operator(op1);
       for(auto it2=operator_sector_2.begin(); it2!=operator_sector_2.end(); ++it2)
    {
     
        auto op2=*it2;
    
       auto v_x=op_dagger;
           v_x.insert(v_x.end(), op2.begin(), op2.end());
               auto [coeff_, nf] = sdp_get_form(v_x);
            
               if(nf.size()%2==0)
               {
               mat_terms.push_back(nf);
               }
    }
    }
return;
}
template <typename T>
 std::unordered_map<int,std::vector<op_vec>> get_mat_terms_double(T &states)
{
    std::unordered_map<int,std::vector<op_vec>> mat_terms;
    for (auto sec : states)
    {

         run_loop(mat_terms[sec.first],states[sec.first][0], states[sec.first][0], sec.first);
     
         run_loop(mat_terms[sec.first],states[sec.first][0], states[sec.first][1], sec.first);
         run_loop(mat_terms[sec.first],states[sec.first][1], states[sec.first][1], sec.first);
     
          
        }

    return mat_terms;
}