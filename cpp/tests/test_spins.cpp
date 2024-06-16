#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>
#include"spins.hpp"
//#include"test.hpp"
using namespace mosek::fusion;
using namespace monty;




int main()
{
       op_vec v1={spin_op("x", {0})};
       op_vec v2={spin_op("y", {0})};
       //op_vec v3={spin_op("y", {2}),spin_op("z", {2}),spin_op("x", {1}),spin_op("y", {0}),spin_op("z", {4}),spin_op("z", {4}),spin_op("z", {4}),spin_op("z", {5})};
       // op_vec v3={spin_op("z", {4}),spin_op("z", {4})};
       op_vec v3={spin_op("x", {1}),spin_op("x", {0}),spin_op("x", {4}),spin_op("x", {0})};
     std::vector<op_vec> v_tot;
     v_tot.push_back(v1); 
    v_tot.push_back(v2);
     
    // v_tot.push_back(v3);
    //  v_tot.push_back(v4);
     
     std::cout<< print_op(v3)<<std::endl;
     //sort(v3.begin(), v3.end());
          std::cout<< print_op(v3)<<std::endl;
	  auto [coeff, op]=get_normal_form(v3);
	  
	  std::cout<< "start "<< op.size()<<std::endl;
	  std::cout<<print_op(op)<<std::endl;
	  std::cout<<coeff<<std::endl;
	  
    //  auto coeff= bubbleSort(v4, v4.size());
    //       std::cout<< print_op(v4)<<std::endl;
    // 	  std::cout<< coeff<<std::endl;
    // 	  std::cout<< "start"<<std::endl;
    //  bool note=false;
    //  auto n=generate_all_terms(v4, note);
    //  for(int l=0; l<n.size(); l++)
    //    {
    // 	 std::cout<<print_op(n[l].second)<< " coeff "<< n[l].first <<std::endl;
    //    }
    //  std::cout<< note<<std::endl;
  //test_commutator();
    // std::string a="up";
    //    std::string b="dn";
    //    bool t=c1>c0;
    //    std::cout<<t<<std::endl;
  // int d=2;
  //  Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });

  // // Setting up the variables
  //  Variable::t rho  = M->variable("X", Domain::inPSDCone(std::pow(d,2)));
  //   Matrix::t C  = Matrix::dense ( new_array_ptr<double, 2>({{2., 1., 0.}, {1., 2., 1.}, {0., 1., 2.}}));
  
}
