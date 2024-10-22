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
using namespace mosek::fusion;
using namespace monty;

std::shared_ptr<ndarray<int,1>>    nint(const std::vector<int> &X)    { return new_array_ptr<int>(X); }
std::shared_ptr<ndarray<double,1>> ndou(const std::vector<double> &X) { return new_array_ptr<double>(X); }

// void define_ham( std::map<std::string, int> terms_mapping, double J, double Delta, int L, std::vector<Expression::t>& block)
// {
//   std::vector<Expression::t> expressions={};
//   std::vector<std::string> sym={"x","y","z"};
//   // since sum_i y_i A_i<=C
 
    
//   for(auto s: sym)
//       {
//     for(int i=0; i<L; i++)
//       {
// 	op_vec v0={spin_op(s, {i}),spin_op(s, {(i+1)%L})};

//      auto [fac0, nf0] =get_normal_form(v0);
	
//      auto el=terms_mapping.at(print_op(nf0));
//      if(s=="z")
//        {

	 
// 	 block[el]=Expr::add(block[el],-1*Delta/4 );

//        }
//      else{
// 	 block[el]=Expr::add(block[el],-1*J/4 );
	 
	 
//       }
//       if((std::abs(fac0.real()-1.)>1e-8 or std::abs(fac0.imag())>1e-8 ))
//        {std::cout<< "error"<<std::endl;}
//       }
// }
   
   
  

//   return;
// }
 Matrix::t define_ham( std::map<std::string, int> terms_mapping, double J, double Delta, int L)
{
  std::vector<Expression::t> expressions={};
  std::vector<std::string> sym={"x","y","z"};
  // since sum_i y_i A_i<=C
 
    std::vector<int> rows;
    std::vector<int> cols;
  std::vector<double> vals;
  for(auto s: sym)
      {
    for(int i=0; i<L; i++)
      {
	op_vec v0={spin_op(s, {i}),spin_op(s, {(i+1)%L})};

     auto [fac0, nf0] =get_normal_form(v0);
	
     auto el=terms_mapping.at(print_op(nf0));
     if(s=="z")
       {

	 
	 //	 block[el]=Expr::add(block[el],-1*Delta/4 );
	 rows.push_back(el);
	 cols.push_back(0);
	  vals.push_back(Delta/4.);
       }
     else{
       //	 block[el]=Expr::add(block[el],-1*J/4 );
	 	 rows.push_back(el);
		 cols.push_back(0);
	  vals.push_back(J/4.);
	 
      }
      if((std::abs(fac0.real()-1.)>1e-8 or std::abs(fac0.imag())>1e-8 ))
       {std::cout<< "error"<<std::endl;}
      }
}
   
   
  
 Matrix::t Alpha_e_1  = Matrix::sparse(terms_mapping.size(),1, nint(rows), nint(cols), ndou(vals)); 

  return Alpha_e_1;
}
// struct matrix_organizer{
//   std::vector<int_pair> matrix_positions;
//   std::vector<double> matrix_values;
  
//   std::vector<double> b;
//   Variable::t variables;
//   int variable_index{0};
  
//   void add_values(int_pair position,double value)
//   {

//     auto index=getIndex(matrix_positions, position);
//     // Check if the target value was found
//     if(index<0)
//       {
// 	matrix_positions.push_back(position);
// 	matrix_values.push_back(value);
	
//       }
//     else
//       {
// 	matrix_values[index]+=value;
//       }
    
//   }
//   Matrix::t make_matrix(int dim1, int dim2)
//   {
//     std::vector<int> rows;
//     std::vector<int> cols;
//     std::vector<double> T;
//     int i=0;
//     for(auto& p:matrix_positions )
//       {
// 	rows.push_back(p.first);
// 	cols.push_back(p.second);
// 	T.push_back(matrix_values[i]);
// 	break;
// 	i++;
//       }
    
  
//   return Matrix::sparse(dim1, dim2, nint(rows), nint(cols), ndou(T));
//   }
// };


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
	 //maps string to the number in the diagonal matrix
;

	 int index_old=0;
	 for(auto a :mat_terms)
	   {
	     
	     terms_mapping.insert({a.first, index_old});
	     index_old++;
	   }

	 
	 // Convention: X= two blocks, X_1, X_2, lower is diagonal. with entries 1+sigma
	 // constratis:
	 //Tr[X_1\gamma]=1
	 // Tr[X_2 \alpha]=2
	 // Tr[\beta X_1]-cTr[delta C_2]=-c
	 

	 std::vector<int> rows_second_block;
	 std::vector<int> cols_second_block;
	 std::vector<double> vals_second_block;
	 
	 std::vector<Expression::t> second_block(mat_terms.size(),Expr::constTerm(0.) );
	 

	 
	 
	 int variable_dimension=2*v_tot.size()*v_tot.size()-v_tot.size()+2*v_tot.size()+2;
	 auto X=M->variable("X", variable_dimension);
	 int variable_index=0;
	 std::vector<double> b(variable_dimension, 0);
	 // one
	 std::vector<double> val={1.,1.};
	std::vector<int> row={0, shift};
	std::vector<int> col={0, shift};
  	Matrix::t Alpha  = Matrix::dense(Matrix::sparse(dim,dim, nint(row), nint(col), ndou(val)));
	auto first_block=Expr::mul(X->index(variable_index),Alpha);
	b[variable_index]+=1;
	variable_index+=1;
	
	// seocnd part
	auto el=terms_mapping.at(print_op({}));
	 b[variable_index]+=2;
   	 second_block[el]=Expr::add(second_block[el],X->index(variable_index));
	 rows_second_block.push_back(el);
	 cols_second_block.push_back(variable_index);
	 vals_second_block.push_back(1.);
	 
	 variable_index+=1;
	 int element=0;
        for(int i=0; i<v_tot.size(); i++)
       {
	 
   for(int j=i; j<v_tot.size(); j++)
     {
       
       auto op_dagger=dagger_operator(v_tot[j]);
       auto v_x=op_dagger;
       v_x.insert(v_x.end(), v_tot[i].begin(),v_tot[i].end());
       auto [coeff, op]=get_normal_form(v_x);
       auto el=terms_mapping.at(print_op(op));
       
       {
   std::vector<double> val={1./2,1./2};
   std::vector<int> row={j+1, j+1+shift };
   std::vector<int> col={i+1,i+1+shift};
   
   Matrix::t Beta_i  = Matrix::sparse(dim,dim, nint(row), nint(col), ndou(val));
   Matrix::t Beta_i_t=Matrix::sparse(dim,dim, nint(col), nint(row), ndou(val));

   second_block[el]=Expr::add(second_block[el], Expr::mul(-1.*coeff.real(),X->index(variable_index)));

   
   first_block=Expr::add(first_block, Expr::mul(Beta_i,X->index(variable_index)));
	first_block=Expr::add(first_block, Expr::mul(Beta_i_t,X->index(variable_index)));
       	
	
	rows_second_block.push_back(el);
	cols_second_block.push_back(variable_index);
	vals_second_block.push_back(-1.*coeff.real());
	b[variable_index]+=-1.*coeff.real();
	
	variable_index+=1;

       }
	   {
	           std::vector<double> val={1./2,-1./2};
	 std::vector<int> row={j+1,shift+j+1 };
	 std::vector<int> col={shift+i+1,i+1};
	 Matrix::t Beta_i  = Matrix::sparse(dim,dim, nint(row), nint(col), ndou(val));
	 Matrix::t Beta_i_t  = Matrix::sparse(dim,dim, nint(col), nint(row), ndou(val));

	second_block[el]=Expr::add(second_block[el], Expr::mul(-1.*coeff.imag(),X->index(variable_index)));
	first_block=Expr::add(first_block, Expr::mul(Beta_i,X->index(variable_index)));
	first_block=Expr::add(first_block, Expr::mul(Beta_i_t,X->index(variable_index)));

	rows_second_block.push_back(el);
	cols_second_block.push_back(variable_index);
	vals_second_block.push_back(-1.*coeff.imag());
	b[variable_index]+=-1.*coeff.imag();
	
	variable_index+=1;
       }
       // if(std::abs(coeff.real())>1e-9 and std::abs(coeff.imag())>1e-9){std::cout<< "error"<<std::endl;}
	   if(j!=i)
	   {
	     std::vector<double> val={1./2,1./2,1./2,1./2};
	     std::vector<int> row={j+1, i+1,j+1+shift,i+1+shift };
	     std::vector<int> col={i+1,j+1,i+1+shift,j+1+shift};
        Matrix::t Beta_i  = Matrix::sparse(dim,dim, nint(row), nint(col), ndou(val));
	second_block[el]=Expr::add(second_block[el], Expr::mul(-1.*coeff.real(),X->index(variable_index)));
	first_block=Expr::add(first_block, Expr::mul(Beta_i,X->index(variable_index)));
	
	rows_second_block.push_back(el);
	cols_second_block.push_back(variable_index);
	vals_second_block.push_back(-1.*coeff.real());
	b[variable_index]+=-1.*coeff.real();
	
	variable_index+=1;

    	
	   }

	 
 	     }
      

            }

	 for(int i=0; i<v_tot.size(); i++)
       {
	 auto [coeff, op]=get_normal_form(v_tot[i]);
	 auto el=terms_mapping.at(print_op(op));
	      
	 {
  std::vector<double> val={1./2,1./2,1./2,1./2};
	     std::vector<int> row={0, i+1,shift,i+1+shift };
	     std::vector<int> col={i+1,0,i+1+shift,shift};
        Matrix::t Beta_i  = Matrix::sparse(dim,dim, nint(row), nint(col), ndou(val));
	 
    	second_block[el]=Expr::add(second_block[el], Expr::mul(-1.*coeff.real(),X->index(variable_index)));
    	first_block=Expr::add(first_block, Expr::mul(Beta_i,X->index(variable_index)));
		rows_second_block.push_back(el);
	cols_second_block.push_back(variable_index);
	vals_second_block.push_back(-1.*coeff.real());
	b[variable_index]=-1.*coeff.real();
	
	variable_index+=1;
        }
	 
   {
         std::vector<double> val={1./2,1./2,-1./2, -1./2};
	 std::vector<int> row={0, i+1,shift,i+1+shift };
	 std::vector<int> col={i+1,0,i+1+shift,shift};
	 Matrix::t Beta_i  = Matrix::sparse(dim,dim, nint(row), nint(col), ndou(val));
	 second_block[el]=Expr::add(second_block[el], Expr::mul(-1.*coeff.imag(),X->index(variable_index)));
	 first_block=Expr::add(first_block, Expr::mul(Beta_i,X->index(variable_index)));
	 rows_second_block.push_back(el);
	cols_second_block.push_back(variable_index);
	vals_second_block.push_back(-1.*coeff.imag());
	b[variable_index]=-1.*coeff.imag();
	
	variable_index+=1;
	 }
 }

	 
	std::cout<< "################################"<<std::endl;
	 double J=1;
  	 double Delta=1;

        
// //  // 	 // add Hamiltonian
        auto C=define_ham( terms_mapping, J, Delta,L);

	Matrix::t variable_mat = Matrix::sparse(mat_terms.size(),variable_dimension, nint(rows_second_block), nint(cols_second_block), ndou(vals_second_block));
	auto tot=Expr::constTerm(mat_terms.size(),0.);
	tot=Expr::add(tot,Expr::mul(X, variable_mat->transpose()));
   //Expr::mul(variable_mat,X);
 // std::cout<<variable_mat->toString()<<std::endl;
	for(int i=0; i<mat_terms.size(); i++)
	  {
	    //M->constraint(second_block[i], Domain::lessThan(C->get(i,0)) );
	    M->constraint(tot->index(i), Domain::lessThan(C->get(i,0)) );
	  }
	//M->constraint(tot, Domain::lessThan(C) );
 		M->constraint(Expr::neg(first_block), Domain::inPSDCone());


	 


	
auto a = monty::new_array_ptr<double>(b);
 auto h=Expr::dot(a, X);
	 
        M->objective(ObjectiveSense::Maximize, h);
	M->dataReport();
	M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
	M->solve();
	  
	  
	    std::cout << "Solution : " << std::endl;
	    std::cout<<std::setprecision(9)<<M->primalObjValue()/L-Delta/4-2*J/4  <<std::endl;
	   
     return 0;
}
