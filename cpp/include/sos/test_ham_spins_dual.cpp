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
#include <Eigen/Sparse>
#include"util.hpp"
using namespace mosek::fusion;
using namespace monty;

/*
solved the dual problem max(b^t y)
\sum_i A_i y_i <=c
*/ 
 std::vector<double> define_ham( std::map<std::string, int> terms_mapping, double J, double Delta, int L)
{
	// defines the vector b_
  std::vector<double> expressions(terms_mapping.size(), 0);
  std::vector<std::string> sym={"x", "y", "z"};
 

  for(auto s: sym)
      {
    for(int i=0; i<L; i++)
      {
	op_vec v0={spin_op(s, {i}),spin_op(s, {(i+1)%L})};


     auto [fac0, nf0] =get_normal_form(v0);
	
     auto el=terms_mapping.at(print_op(nf0));
	
     if(s=="z")
       {
		expressions[el]+=Delta/4.;

       }
     else{
   expressions[el]+=J/4.;
	
	 
      }
      if((std::abs(fac0.real()-1.)>1e-8 or std::abs(fac0.imag())>1e-8 ))
       {std::cout<< "error"<<std::endl;}
      }
}
   

  return expressions;
}
 std::vector<double> define_ham_2( std::map<std::string, int> terms_mapping, double J, double Delta, int L)
{
	// defines the vector b_
  std::vector<double> expressions(terms_mapping.size(), 0);
  std::vector<std::string> sym={"x"};
 

  for(auto s: sym)
      {
    for(int i=0; i<L; i++)
      {
	op_vec v0={spin_op(s, {i})};


     auto [fac0, nf0] =get_normal_form(v0);
	
     auto el=terms_mapping.at(print_op(nf0));
	
     if(s=="z")
       {
		expressions[el]+=Delta/4.;

       }
     else{
   expressions[el]+=J/4.;
	
	 
      }
      if((std::abs(fac0.real()-1.)>1e-8 or std::abs(fac0.imag())>1e-8 ))
       {std::cout<< "error"<<std::endl;}
      }
}
   
   
  

  return expressions;
}



int main()
{
  //test_translation();
    int L=5;
    basis_structure states;
    std::vector<op_vec> v_tot;
      std::vector<std::string> sym={"x", "y", "z"};
      for(int i=0; i<L; i++)
	{
	  for(auto s: sym)
	    {
	      op_vec v1={spin_op(s, {i}), spin_op(s, {(i+1)%L})};
	      auto [fac1, nf1] =get_normal_form(v1);
	     v_tot.push_back(nf1);
	    }
 

	}
             sym={"x", "y", "z"};

      
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
	  mat_terms.insert({print_op(v_tot[i]), 1.});
   for(int j=i; j<v_tot.size(); j++)
     {


       auto op_dagger=dagger_operator(v_tot[j]);

        auto v_x=op_dagger;

       v_x.insert(v_x.end(), v_tot[i].begin(),v_tot[i].end());
       
       	  auto [coeff, op]=get_normal_form(v_x);

	  mat_terms.insert({print_op(op), coeff});

     }
       }

	
	
	// strategy: write two blocks, H1, and H2 (where H2 is diagonal with all variables and 1-H2>>0)
	std::cout<<"mat term size "<< mat_terms.size()<<std::endl;
	Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });
	int dim=2*(v_tot.size()+1);

	int shift=v_tot.size()+1;
	// these are the diagonal elements in H2
		 
	 std::map<std::string, int> terms_mapping;
	 std::map<std::string, matrix_organizer> matrix_mapping;
	 //maps string to the number in the diagonal matrix
;

	 int index_old=0;
	 for(auto a :mat_terms)
	   {
	     
	     terms_mapping.insert({a.first, index_old});
		 matrix_mapping.insert({a.first, matrix_organizer()});
	     index_old++;
	   }
 	auto y=M->variable("T", terms_mapping.size());
	auto el=terms_mapping.at("1");
	  M->constraint( y->index(el), Domain::equalsTo(1.0));
	
	// for(int i=0; i<dim; i++)
	// {
	// matrix_mapping["1"].add_values({i,i},1./2);

	// }
	auto KK_0=matrix_organizer();

	matrix_mapping["1"].add_values({0,0},1./2);
	matrix_mapping["1"].add_values({shift,shift},1./2);

		KK_0.add_values({0,0},1./2);
	KK_0.add_values({shift,shift},1./2);
	
// 	// seocnd part
// 	auto el=terms_mapping.at(print_op({}));
// 	 b[variable_index]+=2;
//    	 second_block[el]=Expr::add(second_block[el],X->index(variable_index));
// 	 rows_second_block.push_back(el);
// 	 cols_second_block.push_back(variable_index);
// 	 vals_second_block.push_back(1.);
	 
// 	 variable_index+=1;
// 	 int element=0;

        for(int i=0; i<v_tot.size(); i++)
       {
	 
   for(int j=i; j<v_tot.size(); j++)
     {
       
       auto op_dagger=dagger_operator(v_tot[j]);
       auto v_x=op_dagger;
       v_x.insert(v_x.end(), v_tot[i].begin(),v_tot[i].end());
       auto [coeff, op]=get_normal_form(v_x);
      
       
        if(std::abs(coeff.real())>1e-9){
			// using 1/4 or 1/2 prefacor?, since H=X_1+X2, but I add A+A^T
		matrix_mapping[print_op(op)].add_values({j+1,i+1},1./2*coeff.real());
		matrix_mapping[print_op(op)].add_values({j+1+shift,i+1+shift},1./2*coeff.real());

		
				KK_0.add_values({j+1,i+1},1./2*coeff.real());
		KK_0.add_values({j+1+shift,i+1+shift},1./2*coeff.real());
if(i!=j)
{
		matrix_mapping[print_op(op)].add_values({i+1,j+1},1./2*coeff.real());
		matrix_mapping[print_op(op)].add_values({i+1+shift,j+1+shift},1./2*coeff.real());

				KK_0.add_values({i+1,j+1},1./2*coeff.real());
		KK_0.add_values({i+1+shift,j+1+shift},1./2*coeff.real());
}

		// change i and j

       }
	   if(std::abs(coeff.imag())>1e-9)
 	   {
			matrix_mapping[print_op(op)].add_values({j+1,shift+i+1},1./2*coeff.imag());
			matrix_mapping[print_op(op)].add_values({i+1,shift+j+1},-1./2*coeff.imag());
		//matrix_mapping[print_op(op)].add_values({shift+j+1,i+1},-1./2*coeff.imag());

		

					KK_0.add_values({j+1,shift+i+1},1./2*coeff.imag());
		KK_0.add_values({i+1,shift+j+1},-1./2*coeff.imag());
		if(i!=j)
		{
	//matrix_mapping[print_op(op)].add_values({i+1,shift+j+1},-1./2*coeff.imag());
		matrix_mapping[print_op(op)].add_values({shift+i+1,j+1},1./2*coeff.imag());
		matrix_mapping[print_op(op)].add_values({shift+j+1,i+1},-1./2*coeff.imag());

			KK_0.add_values({shift+i+1,j+1},1./2*coeff.imag());
		KK_0.add_values({shift+j+1,i+1},-1./2*coeff.imag());
		}

       }
//        // if(std::abs(coeff.real())>1e-9 and std::abs(coeff.imag())>1e-9){std::cout<< "error"<<std::endl;}
// 	   if(j!=i)
// 	   {
// 	     std::vector<double> val={1./2,1./2,1./2,1./2};
// 	     std::vector<int> row={j+1, i+1,j+1+shift,i+1+shift };
// 	     std::vector<int> col={i+1,j+1,i+1+shift,j+1+shift};
//         Matrix::t Beta_i  = Matrix::sparse(dim,dim, nint(row), nint(col), ndou(val));
// 	second_block[el]=Expr::add(second_block[el], Expr::mul(-1.*coeff.real(),X->index(variable_index)));
// 	first_block=Expr::add(first_block, Expr::mul(Beta_i,X->index(variable_index)));
	
// 	rows_second_block.push_back(el);
// 	cols_second_block.push_back(variable_index);
// 	vals_second_block.push_back(-1.*coeff.real());
// 	b[variable_index]+=-1.*coeff.real();
	
// 	variable_index+=1;

    	
// 	   }

	 
 	     }
      

            }

	 for(int i=0; i<v_tot.size(); i++)
       {
	 auto [coeff, op]=get_normal_form(v_tot[i]);
	 auto el=terms_mapping.at(print_op(op));
	      std::cout<< print_op(op)<< " coeff "<<coeff<<std::endl;
 	 if(std::abs(coeff.real())>1e-9){
		std::cout<<matrix_mapping[print_op(op)].make_matrix(dim, dim)->toString()<<std::endl;
		matrix_mapping[print_op(op)].add_values({0,i+1},1./2*coeff.real());
		matrix_mapping[print_op(op)].add_values({i+1,0},1./2*coeff.real());

		KK_0.add_values({0,i+1},1./2*coeff.real());
		KK_0.add_values({i+1,0},1./2*coeff.real());


		 	matrix_mapping[print_op(op)].add_values({shift, i+1+shift },1./2*coeff.real());
		 matrix_mapping[print_op(op)].add_values({i+1+shift,shift},1./2*coeff.real());

		  	KK_0.add_values({shift, i+1+shift },1./2*coeff.real());
		 KK_0.add_values({i+1+shift,shift},1./2*coeff.real());
		 std::cout<<matrix_mapping[print_op(op)].make_matrix(dim, dim)->toString()<<std::endl;
	


         }
	
	 if(std::abs(coeff.imag())>1e-9){std::cout<< "error imaginare non zero"<<std::endl;}

  }

	 
// std::cout<< "################################"<<std::endl;
// 	std::cout<< "dim "<<dim<< " shift "<< shift<< " vsize "<< v_tot.size()<<std::endl;






auto it=matrix_mapping.begin();
auto expression=Expr::mul(y->index(terms_mapping.at(it->first)),matrix_mapping[it->first].make_matrix(dim, dim));

it++;
 int l=1;
while(it!=matrix_mapping.end())
{
	std::cout<<l<< "  and "<<it->first<<"  and "<< terms_mapping.at(it->first)<<std::endl;
	std::cout<<matrix_mapping[it->first].make_matrix(dim, dim)->toString()<<std::endl;
	expression=Expr::add(expression, Expr::mul(y->index(terms_mapping.at(it->first)), matrix_mapping[it->first].make_matrix(dim, dim)));
it++;
l++;

}
std::cout<< "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"<<std::endl;
std::cout<<KK_0.make_matrix(dim, dim)->toString()<<std::endl;


 M->constraint(expression, Domain::inPSDCone());
	
		 double J=1;
  	 double Delta=1;

        
// // //  // 	 // add Hamiltonian
         auto b=define_ham( terms_mapping, J, Delta,L);

 auto b_arr = monty::new_array_ptr<double>(b);
  auto h=Expr::dot(b_arr, y);

        M->objective(ObjectiveSense::Maximize, (h));
	M->dataReport();
	M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
	M->solve();
	  
	  
 	    std::cout << "Solution : " << std::endl;
  //  std::cout<<std::setprecision(9)<<M->primalObjValue() <<std::endl;
		 	    std::cout<<std::setprecision(9)<<M->primalObjValue()/L <<std::endl;
	   
     return 0;
}
/*auto KK_0=matrix_organizer();
int n=4;
for(int j=0;j<2*n; j++)
{

	KK_0.add_values({j,j},1);
}
auto KK_1=matrix_organizer();

KK_1.add_values({0,1},1);
KK_1.add_values({1,0},1);
KK_1.add_values({n,n+1},1);
KK_1.add_values({n+1,n},1);

auto KK_2=matrix_organizer();

KK_2.add_values({0,2},1);
KK_2.add_values({2,0},1);
KK_2.add_values({n,n+2},1);
KK_2.add_values({n+2,n},1);


std::cout<<KK_1.make_matrix(2*n, 2*n)->toString()<<std::endl;
auto expression=Expr::mul(y->index(1), KK_1.make_matrix(2*n, 2*n));
expression=Expr::add(expression,Expr::mul(y->index(0), KK_0.make_matrix(2*n, 2*n)));
expression=Expr::add(expression,Expr::mul(y->index(2), KK_2.make_matrix(2*n, 2*n)));
  M->objective(ObjectiveSense::Minimize, Expr::add(y->index(1),y->index(2)));*/