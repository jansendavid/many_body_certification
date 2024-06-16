#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>
#include"electrons.hpp"
#include"test.hpp"
using namespace mosek::fusion;
using namespace monty;
struct mat_el{
  mat_el(int ind_1, int ind_2, double coeff):i1_(ind_1), i2_(ind_2), coeff_(coeff){}
  int i1_;
  int i2_;
  double coeff_;
  mat_el()=default;
};
std::ostream& operator<<(std::ostream& os, const mat_el& op)
{
  os <<"("<<op.i1_ << ","<<op.i2_<< ") coeff: "<< op.coeff_<<std::endl;
    return os;
}


Matrix::t get_ham(  std::unordered_map<std::string, mat_el> mat_map, int L, int dim)
{
  std::vector<int> rows={};
  std::vector<int> cols={};
  std::vector<double> values={};
 
  for(int i=0; i<L; i++)
    {
       op_vec v1={electron_op("up", {i}, true),electron_op("up", {(i+1)%L}, false)};
       auto coeff= bubbleSort(v1, v1.size());
       //std::cout<< <<std::endl;
       auto el=mat_map[print_op(v1)];
       rows.push_back(el.i1_);
       cols.push_back(el.i2_);
       values.push_back(-1*coeff*el.coeff_);
       //std::cout<< el << print_op(v1)<< " tot cof"<< -1*coeff*el.coeff_<<std::endl;
       v1={electron_op("up", {(i+1)%L}, true),electron_op("up", {i}, false)};
        coeff= bubbleSort(v1, v1.size());
       
        el=mat_map[print_op(v1)];
	//	std::cout<< el << print_op(v1)<< " tot cof"<< -1*coeff*el.coeff_<<std::endl;
       rows.push_back(el.i1_);
       cols.push_back(el.i2_);
       values.push_back(-1*coeff*el.coeff_);
    }
     

       
  auto rows_pt=monty::new_array_ptr<int>(rows);
  auto cols_pt=monty::new_array_ptr<int>(cols);
  auto values_pt=monty::new_array_ptr<double>(values);
   auto ms = Matrix::sparse(dim, dim, rows_pt, cols_pt, values_pt);
   std::cout<< "el numer "<< values.size()<<std::endl;
   return ms;
}

int main()
{


  
   int L=3;
   
  std::vector<op_vec> v_tot;
  std::unordered_map<std::string, int> vec_map;
  int indi=0;
  for(int i=0; i<L; i++)
    {
      op_vec v1={electron_op("up", {i}, true)};
      op_vec v2={electron_op("up", {i}, false)};
      v_tot.push_back(v1);
      vec_map.insert({print_op(v1), indi});
      indi+=1;
      
      v_tot.push_back(v2);
      vec_map.insert({print_op(v2), indi});
      indi+=1;
    }
  // for(auto a: v_tot)
  //   {
  //     std::cout<<print_op(a)<<std::endl;
  //   }
     
     std::unordered_map<std::string, int> terms;
     std::unordered_map<std::string, mat_el> mat_terms;
     
     indi=0;
     for(int i=0; i<v_tot.size(); i++)
       {
	 mat_terms.insert({print_op(v_tot[i]), mat_el(0,i+1, 1)});
   for(int j=0; j<v_tot.size(); j++)
     {
       auto op_dagger=dagger_operator(v_tot[j]);
       int coeff= bubbleSort(op_dagger, op_dagger.size());
       
       mat_terms.insert({print_op(op_dagger), mat_el(j+1,0, coeff)});

    
       auto v_x=op_dagger;
       v_x.insert(v_x.end(), v_tot[i].begin(),v_tot[i].end());
            coeff*= bubbleSort(v_x, v_x.size());
              bool dummy_variable=false;
       auto all_terms=generate_all_terms(v_x,dummy_variable);
       if(all_terms.size()==1 and !dummy_variable)
	 {
	   mat_terms.insert({print_op(all_terms[0].second), mat_el(j+1,i+1, all_terms[0].first*coeff)});
	 }
       
     }
     }



     
     for(auto a: mat_terms)
       {
	 // std::cout<< a.first<< " and "<< a.second<<std::endl;
       }

  Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });	
      auto X=M->variable("X", Domain::inPSDCone(v_tot.size()+1));
     M->constraint( X->index(0,0), Domain::equalsTo(1.0));


  for(int i=0; i<v_tot.size(); i++)
       {
	 
   for(int j=i; j<v_tot.size(); j++)
     {
       auto op_dagger=dagger_operator(v_tot[j]);
       int coeff= bubbleSort(op_dagger, op_dagger.size());
       if(print_op(v_tot[i])==print_op(op_dagger))
	 {
	   std::cout<< "i,j "<< i+1 << "  "<< j+1<<   " coeff "<< coeff<<std::endl;
	   M->constraint( Expr::add(X->index(0,i+1),Expr::mul(-1.,X->index(j+1,0))), Domain::equalsTo(0.0));
	   M->constraint( Expr::add(X->index(0,j+1),Expr::mul(-1.,X->index(i+1,0))), Domain::equalsTo(0.0));
	 }
     }}
  
    for(int i=0; i<v_tot.size(); i++)
       {
	 
   for(int j=0; j<v_tot.size(); j++)
     {
       auto op_dagger=dagger_operator(v_tot[j]);
       int coeff= bubbleSort(op_dagger, op_dagger.size());
       
       

    
       auto v_x=op_dagger;
       v_x.insert(v_x.end(), v_tot[i].begin(),v_tot[i].end());
       coeff*= bubbleSort(v_x, v_x.size());
       bool dummy_variable=false;
       auto all_terms=generate_all_terms(v_x,dummy_variable);
       
       auto el=mat_terms[print_op(all_terms[0].second)];
       if(el.i1_!=j+1 and el.i2_!=i+1 and all_terms.size()>0 and !dummy_variable)
	 {
	   //  std::cout<< j+1<< " " << i+1 <<" "<< print_op(all_terms[0].second)<< " " <<el<<std::endl;
	    M->constraint( Expr::add(X->index(j+1,i+1),Expr::mul(1.*coeff,X->index(el.i1_,el.i2_))), Domain::equalsTo(0.0));
	 }
       
       
     }
     }


    
  
   for(int i=1; i<v_tot.size(); i+=2)
       {
	 //std::cout<< "pairs "<< i << " "<< i+1<<std::endl;
	 M->constraint( Expr::add(X->index(i+1,i+1),Expr::add(Expr::mul(1.,X->index(i,i)), -1)), Domain::equalsTo(0.0));
	 
       }
  //for (auto& [key, value]:mat_terms) {  std::cout << key << "  "<< mat_terms[key]<< std::endl; }
    
   auto Ms=get_ham(mat_terms, L,v_tot.size()+1);

   M->objective( ObjectiveSense::Minimize,  Expr::dot(Ms,X) );
   //M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
      M->solve();
  std::cout<< M->primalObjValue()  <<std::endl;
std::cout << "  X = " << *(X->level()) << std::endl;
 

     
 

      
}
