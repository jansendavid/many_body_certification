#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>
#include"spins.hpp"
#include<unordered_map>
#include <Eigen/Dense>
#include"complex_momentum_symm.hpp"
using namespace mosek::fusion;
using namespace monty;


Expression::t define_ham( std::map<std::string, Variable::t> terms_mapping, double J, double Delta, int L)
{
  std::vector<Expression::t> expressions={};
  std::vector<std::string> sym={"x","y","z"};
  auto ham=Expr::constTerm(0.);
 
    
  for(auto s: sym)
      {
    for(int i=0; i<L; i++)
      {
	op_vec v0={spin_op(s, {i}),spin_op(s, {(i+1)%L})};
     auto [fac0, nf0] =get_normal_form(v0);
     

     if(s=="z")
       {
	 ham=Expr::add(ham,Expr::mul(Delta/4,terms_mapping.at(print_op(nf0))));

       }
     else{
	 ham=Expr::add(ham,Expr::mul(J/4,terms_mapping.at(print_op(nf0))));
      }
      if((std::abs(fac0.real()-1.)>1e-8 or std::abs(fac0.imag())>1e-8 ))
       {std::cout<< "error"<<std::endl;}
      }
}
   
   
  

  return ham;
}


int main()
{
  //test_translation();
    int L=5;
    basis_structure states;
    std::vector<op_vec> v_tot;
      std::vector<std::string> sym={"z", "x", "y"};
      for(int i=0; i<L; i++)
	{
	  for(auto s: sym)
	    {
	      op_vec v1={spin_op(s, {i}), spin_op(s, {(i+1)%L})};
	      auto [fac1, nf1] =get_normal_form(v1);
	      v_tot.push_back(nf1);
	    }
 

	}
             sym={"x","y","z"};
      
      for(int i=0; i<L; i++)
	{
	  for(auto s: sym)
	    {
	      op_vec v1={spin_op(s, {i})};
	      v_tot.push_back(v1);
	    }
 

	}
      std::vector<std::pair<std::string, std::string>> occ={std::pair<std::string, std::string>("x","y"),std::pair<std::string, std::string>("y","x"),std::pair<std::string, std::string>("y","z"),std::pair<std::string, std::string>("z","y"),std::pair<std::string, std::string>("x","z"),std::pair<std::string, std::string>("z","x")};
   for(int i=0; i<L; i++)
	{
	  for(auto& a: occ)
	    {
	      op_vec v1={spin_op(a.first, {i}), spin_op(a.second, {(i+1)%L})};
	      auto [fac1, nf1] =get_normal_form(v1);
	      v_tot.push_back(nf1);
	    }
	      

	}



       std::map<std::string, std::complex<double>> mat_terms;
      
       
        for(int i=0; i<v_tot.size(); i++)
       {
	 
   for(int j=i; j<v_tot.size(); j++)
     {

 
 auto op_dagger=dagger_operator(v_tot[j]);
 
        auto v_x=op_dagger;

       v_x.insert(v_x.end(), v_tot[i].begin(),v_tot[i].end());
       
       	  auto [coeff, op]=get_normal_form(v_x);
        
	  mat_terms.insert({print_op(op), coeff});

     }
       }
	std::cout<< mat_terms.size()<<std::endl;
	Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });	
	   auto X=M->variable("X", Domain::inPSDCone(2*(v_tot.size()+1)));
	   int shift=v_tot.size()+1;
	   M->constraint( Expr::add(X->index(0,0), X->index(shift,shift)), Domain::equalsTo(1.0));
	  
	 auto XT=M->variable("T", mat_terms.size(),Domain::inRange(-1., 1));
	 
	 std::map<std::string, Variable::t> terms_mapping;
	 int index=0;
	 for(auto a :mat_terms)
	   {
	     
	     terms_mapping.insert({a.first, XT->index(index)});
	     index++;
	   }
	 
	 auto el=terms_mapping.at(print_op({}));
	 M->constraint( el, Domain::equalsTo(1.0));
        for(int i=0; i<v_tot.size(); i++)
       {
	 
   for(int j=i; j<v_tot.size(); j++)
     {
       
       auto op_dagger=dagger_operator(v_tot[j]);
        auto v_x=op_dagger;
       v_x.insert(v_x.end(), v_tot[i].begin(),v_tot[i].end());
       auto [coeff, op]=get_normal_form(v_x);
       
       auto el=terms_mapping.at(print_op(op));
       
	M->constraint( Expr::add(Expr::add(X->index(i+1,j+1),X->index(i+1+shift,j+1+shift)),Expr::mul(-1.*coeff.real(),el)), Domain::equalsTo(0.0));
	

       
		  M->constraint( Expr::add(Expr::add(X->index(shift+i+1,j+1),Expr::mul(-1.,X->index(shift+j+1,i+1))),Expr::mul(-1.*coeff.imag(),el)), Domain::equalsTo(0.0)); 
	  
	  if(j!=i)
	    {
	       M->constraint( Expr::add(Expr::add(X->index(j+1,i+1),X->index(j+1+shift,i+1+shift)),Expr::mul(-1.*coeff.real(),el)), Domain::equalsTo(0.0));
	      M->constraint( Expr::add(Expr::add(X->index(shift+j+1,i+1),Expr::mul(-1.,X->index(shift+i+1,j+1))),Expr::mul(1.*coeff.imag(),el)), Domain::equalsTo(0.0));

	    }

       }}

	 	 
	       for(int i=0; i<v_tot.size(); i++)
       {
	 auto [coeff, op]=get_normal_form(v_tot[i]);
	 auto el=terms_mapping.at(print_op(op));
	 
	 M->constraint( Expr::add(Expr::add(X->index(0,i+1),X->index(shift,i+1+shift)),Expr::mul(-1.*coeff.real(),el)), Domain::equalsTo(0.0));

		  M->constraint( Expr::add(Expr::add(X->index(shift,i+1),X->index(shift+i+1,0)),Expr::mul(-1.*coeff.imag(),el)), Domain::equalsTo(0.0));

	 

       }

	       
	 double J=1;
   	double Delta=1;
   	auto h=define_ham( terms_mapping, J, Delta,L);
	
	for(auto a : terms_mapping)
	  {
	    //std::cout<< a.first<<std::endl;
	  }	
  // 	// double numer_of_part=0.;
	// M->constraint( Expr::mul(L,part_nr_op), Domain::equalsTo(numer_of_part));
        M->objective(ObjectiveSense::Minimize, h);
	 M->dataReport();
	  M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
	   M->solve();
	  
	  
	  std::cout << "Solution : " << std::endl;
	  std::cout<<std::setprecision(9)<<M->primalObjValue()/L  <<std::endl;
 
	  //auto x=Matrix::dense(block.total_refs_.size(), 1, block.XT_->level());
     return 0;
}
