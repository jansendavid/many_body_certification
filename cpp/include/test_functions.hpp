#pragma once
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xnpy.hpp>
using namespace mosek::fusion;
using namespace monty;
#include"functions.hpp"
void test_write_to_numpy()
{
  /*
    Writing the matrix
1, 3+2j, 2.3 -1j
3-2j, 3, 6.3-4.4j,
2.3 +1j, 6.3+4.4j, 3.3 

   */
  //Matrix::t C  = Matrix::dense ( new_array_ptr<double, 2>({{1., 3., 2.3, 0, 2, -1.}, {3., 3., 6.3, -2, 0., -4.4}, {2.3, 6.3, 3.3,1.,4.4,0}, {0., -2., 1.,1.,3.,2.3}, {2., 0., 4.4,3.,3.,6.3}, {-1., -4.4, 0.,2.3,6.3,3.3}}));
  int dim=3;

  Model::t M_ = new Model("sdo1"); auto _M = finally([&]() { M_->dispose(); });
  auto X=M_->variable("X", Domain::inPSDCone(2*dim));
  M_->constraint( X->index(0,0), Domain::equalsTo(1.48201916));
  M_->constraint( X->index(0,1), Domain::equalsTo(0.95919097));
  M_->constraint( X->index(0,2), Domain::equalsTo(1.04858386));
  M_->constraint( X->index(0,3), Domain::equalsTo(0.));
  M_->constraint( X->index(0,4), Domain::equalsTo(1.57363419e-01));
  M_->constraint( X->index(0,5), Domain::equalsTo(8.70121947e-01));
   
  M_->constraint( X->index(1,0), Domain::equalsTo(0.95919097));
  M_->constraint( X->index(1,1), Domain::equalsTo(0.96956815));
  M_->constraint( X->index(1,2), Domain::equalsTo(0.99608958));
  M_->constraint( X->index(1,3), Domain::equalsTo(-1.57363419e-01));
  M_->constraint( X->index(1,4), Domain::equalsTo(0.));
  M_->constraint( X->index(1,5), Domain::equalsTo(+7.52863846e-01));
  
  M_->constraint( X->index(2,0), Domain::equalsTo(1.04858386));
  M_->constraint( X->index(2,1), Domain::equalsTo(0.99608958));
  M_->constraint( X->index(2,2), Domain::equalsTo(1.69414421));
  M_->constraint( X->index(2,3), Domain::equalsTo(-8.70121947e-01));
  M_->constraint( X->index(2,4), Domain::equalsTo(-7.52863846e-01));
  M_->constraint( X->index(2,5), Domain::equalsTo(0.));
  	   for(int l=0; l<dim;l++)
	     {
	       for(int m=l;m<dim; m++)
		 {
		   
		   M_->constraint( Expr::add(X->index(l,m),Expr::mul(-1.,X->index(l+dim,m+dim))), Domain::equalsTo(0.0));
		    M_->constraint( Expr::add(X->index(l,m+dim),Expr::mul(1.,X->index(m,l+dim))), Domain::equalsTo(0.0));
		    

		 }
	     }
	   M_->objective( ObjectiveSense::Minimize,  X->index(0,0) );
   
   //M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
      M_->solve();
  std::cout<< M_->primalObjValue()  <<std::endl;
    store_matrix_to_numpy("to_numpy_test.npy", 3, X);
  return; }

void test_extract_basis_from_file()
{
  int L=4;

 auto basis=load_basis_from_file("test_basis_bilayer.txt", L, true);
 for(auto b: basis)
 {
std::cout<<print_op(b)<<std::endl;

 }
  return;
}