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

void test_translation()
{
  int L=3;
   
    op_vec v1={spin_op("x", {0}),spin_op("y", {1})};

    //std::cout<<print_op(v1)<<std::endl;
    // for(int i=0; i<=L;i++)
    //   {
    // 	auto O=translation(v1, i, L);
    // 	//	std::cout<<print_op(O)<<std::endl;
    //   }
    auto all_T=generate_all_translations(v1, L);
    for(auto& a: all_T)
      {
	std::cout<<print_op(a)<<std::endl;
      }
    
}
Expression::t define_ham( std::map<std::string, Variable::t> terms_mapping, double J, double Delta, int L)
{
  std::vector<Expression::t> expressions={};
  std::vector<std::string> sym={"x","y","z"};
  auto ham=Expr::constTerm(0.);
    for(auto s: sym)
      {
	op_vec v0={spin_op(s, {0}),spin_op(s, {1})};
     auto [fac0, nf0] =get_normal_form(v0);
     
    
     if(s=="z")
       {
	 ham=Expr::add(ham,Expr::mul(Delta/4,terms_mapping.at(print_op(nf0))));

       }
     else{
	 ham=Expr::add(ham,Expr::mul(J/4,terms_mapping.at(print_op(nf0))));

     }
     if((std::abs(fac0.real()-1.)>1e-8 or std::abs(fac0.imag())>1e-8 )  )
       {std::cout<< "error"<<std::endl;}
      }
    
  for(auto s: sym)
      {
    for(int i=1; i<L; i++)
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
    int L=4;
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
             sym={"z"};
      
      for(int i=0; i<L; i++)
	{
	  for(auto s: sym)
	    {
	      op_vec v1={spin_op(s, {i})};
	      v_tot.push_back(v1);
	    }
 

	}
      std::vector<std::pair<std::string, std::string>> occ={std::pair<std::string, std::string>("x","y")};//,std::pair<std::string, std::string>("y","x")};//,std::pair<std::string, std::string>("y","z"),std::pair<std::string, std::string>("z","y"),std::pair<std::string, std::string>("x","z"),std::pair<std::string, std::string>("z","x")};
   for(int i=0; i<L; i++)
	{
	  for(auto& a: occ)
	    {
	      op_vec v1={spin_op(a.first, {i}), spin_op(a.second, {(i+1)%L})};
	      auto [fac1, nf1] =get_normal_form(v1);
	      v_tot.push_back(nf1);
	    }
	      

	}

   // occ={std::pair<std::string, std::string>("x","y"),std::pair<std::string, std::string>("y","x")};
   // for(int i=0; i<L; i++)
   // 	{
   // 	  for(auto& a: occ)
   // 	    {
   // 	      op_vec v1={spin_op(a.first, {i}), spin_op(a.second, {(i+2)%L})};
   // 	      auto [fac1, nf1] =get_normal_form(v1);
   // 	      v_tot.push_back(nf1);
   // 	    }
	      

   //	}
   
   //std::cout<< "size "<< v_tot.size()<< " and "<< occ.size()*L+2*L*sym.size()<<std::endl;

       std::map<std::string, std::complex<double>> mat_terms;
      
       
        for(int i=0; i<v_tot.size(); i++)
       {
	 
   for(int j=i; j<v_tot.size(); j++)
     {

       //std::cout<< "j "<< print_op(v_tot[j])<<std::endl;
 auto op_dagger=dagger_operator(v_tot[j]);
 // std::cout<< "dag "<< print_op(op_dagger)<<std::endl;
        auto v_x=op_dagger;

       v_x.insert(v_x.end(), v_tot[i].begin(),v_tot[i].end());
       
       	  auto [coeff, op]=get_normal_form(v_x);
	  //	  std::cout<< "tot "<< print_op(v_x)<< "  nf-> "<< print_op(op)<< "  coeff "<< coeff<< std::endl;
	  //if(op.size()<1)
	  //{std::cout<<"dd "<< coeff<<std::endl;}
	  mat_terms.insert({print_op(op), coeff});

     }
       }

	Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });	
	   auto X=M->variable("X", Domain::inPSDCone(2*(v_tot.size()+1)));
	   int shift=v_tot.size()+1;
	   M->constraint( X->index(0,0), Domain::equalsTo(1.0));
	   M->constraint( X->index(shift,shift), Domain::equalsTo(1.0));
	 auto XT=M->variable("T", mat_terms.size(),Domain::inRange(-1., 1));
	 
	 std::map<std::string, Variable::t> terms_mapping;
	 int index=0;
	 // for(auto a :mat_terms)
	 //   {

	 //     std::cout<<"x"<<std::endl;
	 //     std::cout<<a.first<<std::endl;
	 //   }
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
	  // if(op.size()==0){std::cout<<"h'"<<std::endl;}
	  auto el=terms_mapping.at(print_op(op));
	  M->constraint( Expr::add(X->index(i+1,j+1),Expr::mul(-1.*coeff.real(),el)), Domain::equalsTo(0.0));
	  M->constraint( Expr::add(X->index(shift+i+1,j+1),Expr::mul(-1.*coeff.imag(),el)), Domain::equalsTo(0.0));
	  if(j!=i)
	    {
	      M->constraint( Expr::add(X->index(j+1,i+1),Expr::mul(-1.*coeff.real(),el)), Domain::equalsTo(0.0));
	  M->constraint( Expr::add(X->index(shift+j+1,i+1),Expr::mul(1.*coeff.imag(),el)), Domain::equalsTo(0.0));
	    }

     }}

	       for(int i=0; i<v_tot.size(); i++)
       {
	 auto [coeff, op]=get_normal_form(v_tot[i]);
	 auto el=terms_mapping.at(print_op(op));
	 
	  M->constraint( Expr::add(X->index(0,i+1),Expr::mul(-1.*coeff.real(),el)), Domain::equalsTo(0.0));
	  M->constraint( Expr::add(X->index(i+1,0),Expr::mul(-1.*coeff.real(),el)), Domain::equalsTo(0.0));
	  M->constraint( Expr::add(X->index(shift,i+1),Expr::mul(-1.*coeff.imag(),el)), Domain::equalsTo(0.0));
	  M->constraint( Expr::add(X->index(shift+i+1,0),Expr::mul(-1.*coeff.imag(),el)), Domain::equalsTo(0.0));
	 

       }
	 
	   for(int l=0; l<shift;l++)
	     {
	       for(int m=l;m<shift; m++)
		 {
		   
		   M->constraint( Expr::add(X->index(l,m),Expr::mul(-1.,X->index(l+shift,m+shift))), Domain::equalsTo(0.0));
		    M->constraint( Expr::add(X->index(l,m+shift),Expr::mul(1.,X->index(m,l+shift))), Domain::equalsTo(0.0));

		 }
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
	//   // 	  block.M_->dataReport();
	//   // M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
	   M->solve();
	  
	  
	  std::cout << "Solution : " << std::endl;
	  std::cout<<std::setprecision(9)<<M->primalObjValue()/L  <<std::endl;
 
	  //auto x=Matrix::dense(block.total_refs_.size(), 1, block.XT_->level());
     return 0;
}
