#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>
#include"electrons.hpp"
#include"test.hpp"
#include <iomanip>
#include"functions.hpp"
using namespace mosek::fusion;
using namespace monty;




Matrix::t get_ham(  std::unordered_map<std::string, mat_el_ref> mat_map, int L, int dim, double t0, double U)
{
  std::vector<int> rows={};
  std::vector<int> cols={};
  std::vector<double> values={};
 
  for(int i=0; i<L; i++)
    {

       op_vec v1={electron_op("up", {i}, true),electron_op("up", {(i+1)%L}, false)};
       auto coeff= bubbleSort(v1, v1.size());

       auto el=mat_map[print_op(v1)];
       rows.push_back(el.i1_);
       cols.push_back(el.i2_);
       values.push_back(-1*coeff*el.coeff_*t0);

       v1={electron_op("up", {(i+1)%L}, true),electron_op("up", {i}, false)};
        coeff= bubbleSort(v1, v1.size());
       
        el=mat_map[print_op(v1)];

       rows.push_back(el.i1_);
       cols.push_back(el.i2_);
       values.push_back(-1*coeff*el.coeff_*t0);
    }
     
  // for(int i=0; i<L; i++)
  //   {

  //     op_vec v1={electron_op("up", {i}, true),electron_op("up", {i}, false),electron_op("up", {(i+1)%L}, true),electron_op("up", {(i+1)%L}, false)};
  //      auto coeff= bubbleSort(v1, v1.size());
       
  //      auto el=mat_map[print_op(v1)];
  //      rows.push_back(el.i1_);
  //      cols.push_back(el.i2_);
  //      values.push_back(U*coeff*el.coeff_);
  //   }

  
       
  auto rows_pt=monty::new_array_ptr<int>(rows);
  auto cols_pt=monty::new_array_ptr<int>(cols);
  auto values_pt=monty::new_array_ptr<double>(values);
   auto ms = Matrix::sparse(dim, dim, rows_pt, cols_pt, values_pt);

   return ms;
}
void fix_particle_number(Model::t& M,int L, int part_nr,std::unordered_map<std::string, mat_el_ref> mat_map)
{
  op_vec v1={electron_op("up", {0}, true),electron_op("up", {0}, false)};
  auto el=mat_map[print_op(v1)];
  auto expr=Expr::mul(1*el.coeff_,el.var_);
  std::cout<< print_op(v1)<< " "<< el<<std::endl;
  for(int i=1; i<L; i++)
    {
       v1={electron_op("up", {i}, true),electron_op("up", {i}, false)};
       el=mat_map[print_op(v1)];
       //std::cout<< print_op(v1)<< " "<< el<<std::endl;
       
       expr=Expr::add(Expr::mul(1*el.coeff_,el.var_), expr);
    }
  
  M->constraint( expr, Domain::equalsTo(part_nr));
  return;
}
std::vector<op_vec> make_basis(int L)
{
  std::vector<op_vec> v_tot;

  for(int i=0; i<L; i++)
    {
      op_vec v1={electron_op("up", {i}, true)};
      op_vec v2={electron_op("up", {i}, false)};
       v_tot.push_back(v1);      
        v_tot.push_back(v2);
      

    }
 
//     for(int i=0; i<L; i++)
//     {
//       op_vec v1={electron_op("up", {i}, true),electron_op("up", {i}, false)};
//       v_tot.push_back(v1);      

      

//     }
//     for(int i=0; i<L; i++)
//     {
//       op_vec v1={electron_op("up", {i}, true),electron_op("up", {(i+1)%L}, false)};
//       auto coeff= bubbleSort(v1, v1.size());
//       v_tot.push_back(v1);
//     v1={electron_op("up", {i}, false),electron_op("up", {(i+1)%L}, true)};
//       coeff= bubbleSort(v1, v1.size());
//       v_tot.push_back(v1); 
// }

//         for(int i=0; i<L; i++)std::vector
//     {
//       op_vec v1={electron_op("up", {i}, true),electron_op("up", {(i+2)%L}, false)};
//       auto coeff= bubbleSort(v1, v1.size());
//       v_tot.push_back(v1);
//       v1={electron_op("up", {i}, false),electron_op("up", {(i+2)%L}, true)};
//       coeff= bubbleSort(v1, v1.size());
//       v_tot.push_back(v1);      
// }


//        for(int i=0; i<L; i++)
//     {
//       op_vec v1={electron_op("up", {i}, true),electron_op("up", {(i+3)%L}, false)};
//       auto coeff= bubbleSort(v1, v1.size());
//       v_tot.push_back(v1);
//       v1={electron_op("up", {i}, false),electron_op("up", {(i+3)%L}, true)};
//       coeff= bubbleSort(v1, v1.size());
//       v_tot.push_back(v1);      
// }	
 
//      for(int i=0; i<L; i++)
//     {
//       op_vec v1={electron_op("up", {i}, true),electron_op("up", {(i+4)%L}, false)};
//       auto coeff= bubbleSort(v1, v1.size());
//       v_tot.push_back(v1);
//       v1={electron_op("up", {i}, false),electron_op("up", {(i+4)%L}, true)};
//       coeff= bubbleSort(v1, v1.size());
//       v_tot.push_back(v1);      
// }    

	
           

	return v_tot;}
int main()
{


  
   int L=3;
   
   auto v_tot=make_basis(L);
    
  // for(auto a: v_tot)
  //   {
  //     std::cout<<print_op(a)<<std::endl;
  //   }
     
     std::unordered_map<std::string, int> terms;
     std::unordered_map<std::string, mat_el> mat_terms;
     std::unordered_map<std::string, mat_el> extra_terms;
     std::unordered_map<std::string, mat_el_ref> total_refs;
     std::vector<op_vec> extra_terms_vec;
     
     int indi=0;
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
       else{
	 for(int l=0; l<all_terms.size(); l++)
	   {
	     
	     if (mat_terms.find(print_op(all_terms[l].second)) == mat_terms.end()){
	    
	     extra_terms_vec.push_back(all_terms[l].second);

	     }
	   }
       }
       
     }
     }

     int new_ind=0;
     for(auto a: extra_terms_vec)
       {
	 if (mat_terms.find(print_op(a)) == mat_terms.end() and extra_terms.find(print_op(a)) == extra_terms.end() ){
	   extra_terms.insert({print_op(a), mat_el(new_ind,v_tot.size()+1, 1)});
	   std::cout<< "ex "<< print_op(a)<< "  " << new_ind<< std::endl;
	   new_ind++;
	 }
       }
 
     
     std::cout<<extra_terms.size() << " and "<< mat_terms.size()<<std::endl;
        for(auto a: mat_terms)
        {
	  //std::cout<< a.first<< " and "<< a.second<<std::endl;
        }
     std::cout<< "start"<<std::endl;  
   Model::t M = new Model("sdo1"); auto _M = finally([&]() { M->dispose(); });	
       auto X=M->variable("X", Domain::inPSDCone(v_tot.size()+1));
       M->constraint( X->index(0,0), Domain::equalsTo(1.0));
       //       auto XT=M->variable("XT", extra_terms.size(),Domain::unbounded());
       auto XT=M->variable("XT", extra_terms.size(),Domain::inRange(-1., 1));
       
         for(auto a: mat_terms)
       {
	 auto el=a.second;
	 total_refs.insert({a.first,mat_el_ref(X->index(el.i1_,el.i2_),el.i1_,el.i2_,el.coeff_) });
	 

       }

	 for(auto a: extra_terms)
       {
	 auto el=a.second;
	 //std::cout<< "extra "<<a.first << " "<< el<<std::endl; 
	 total_refs.insert({a.first,mat_el_ref(XT->index(el.i1_),el.i1_,el.i2_,el.coeff_) });
	 

       }  
     

  for(int i=0; i<v_tot.size(); i++)
       {
	 
   for(int j=i; j<v_tot.size(); j++)
     {
       auto op_dagger=dagger_operator(v_tot[j]);
       int coeff= bubbleSort(op_dagger, op_dagger.size());
       if(print_op(v_tot[i])==print_op(op_dagger))
	 {
	   //  std::cout<< "i,j "<< i+1 << "  "<< j+1<<   " coeff "<< coeff<<std::endl;
	   M->constraint( Expr::add(X->index(0,i+1),Expr::mul(-1.*coeff,X->index(j+1,0))), Domain::equalsTo(0.0));
	   M->constraint( Expr::add(X->index(0,j+1),Expr::mul(-1.*coeff,X->index(i+1,0))), Domain::equalsTo(0.0));
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
       auto v_x_copy=v_x;
       coeff*= bubbleSort(v_x, v_x.size());
       bool dummy_variable=false;
       auto all_terms=generate_all_terms(v_x,dummy_variable);
       
       
       
 //       //       std::cout<< print_op(v_x)<< " bool "<< dummy_variable<< " size "<< all_terms.size()<<std::endl;
       
       if( all_terms.size()>0)
	 {


  	   std::vector<Expression::t> expressions={};

  	   auto el=total_refs[print_op(all_terms[0].second)];//mat_terms[print_op(all_terms[0].second)];
	   //std::cout<< j+1 << "  " << i+1 << " ro be "<<  el<<std::endl;
 // 	     //
 // 	   //	    expression;
	 
	 
	 
  	       if(j+1!=el.i1_ or i+1!=el.i2_)		 
  		    {
		      // std::cout<< "enter "<<std::endl;
		      expressions.push_back(Expr::mul(-1.*el.coeff_*all_terms[0].first,el.var_));
		        if(check_zero(all_terms[0].second))
			{
			  
			  M->constraint( el.var_, Domain::equalsTo(0.0));
			}

  		      	//	      std::cout<<" coeff "<<print_op(v_x_copy)<<"  "<< -1.*el.coeff_*all_terms[0].first << "  "<< print_op(all_terms[0].second)<< " "<<"dummy "<<dummy_variable<< std::endl;
 // // std::cout<< el.coeff_<< "  "<<all_terms[0].first <<std::endl;
  		      	   if(dummy_variable){std::cout<< "-1"<<std::endl;}
			   // std::cout<< j+1 << "  " << i+1 << " ro be "<<  el<< " and  "<<all_terms[0].first <<std::endl;
 		    }
	            
	       //		     std::cout<< j+1 << "  " << i+1 << " ro be "<<  el<< " "<<std::endl;
	       for(int l=1; l<all_terms.size(); l++)
		 {
  	   auto el=total_refs[print_op(all_terms[l].second)];

  	       if(j+1!=el.i1_ or i+1!=el.i2_)
  		    {
		        if(check_zero(all_terms[l].second))
			{
			  
			  M->constraint( el.var_, Domain::equalsTo(0.0));
			}
		     
  		      expressions.push_back(Expr::mul(-1.*el.coeff_*all_terms[l].first,el.var_));
   // 	 	      std::cout << " ex fac "<<std::endl;
   // std::cout<<" coeff "<<print_op(v_x_copy)<<"  "<< -1.*el.coeff_*all_terms[l].first << "  "<< print_op(all_terms[l].second)<< " "<<std::endl;
   // std::cout<< j+1 << "  " << i+1 << " ro be "<<  el<< "  "<<std::endl;
	   
		    }
  		 }
		   if(expressions.size()>0)
		     {
		        auto expressions_pt=monty::new_array_ptr<Expression::t>(expressions);
			Expression::t exp=Expr::add(expressions_pt);
				   if(dummy_variable)
	   {
	     exp=Expr::add(-1., exp);
	   }
				   
				   M->constraint( Expr::add(Expr::mul(coeff,X->index(j+1,i+1)),exp), Domain::equalsTo(0.0));
				   
		     }
	
	       
	   }
       
       
     }
      }


    //    op_vec v4={electron_op("up", {0}, true),electron_op("up", {1}, true),electron_op("up", {0}, true)};
    int part_nr=1;
    double t0=1;
    double U=1;

    auto Ms=get_ham(total_refs, L,v_tot.size()+1, t0, U);
    fix_particle_number(M, L, part_nr,total_refs);
   M->objective( ObjectiveSense::Minimize,  Expr::dot(Ms,X) );
   
   //M->setLogHandler([=](const std::string & msg) { std::cout << msg << std::flush; } );
      M->solve();
  std::cout<< M->primalObjValue()  <<std::endl;
 // // //  //std::cout << "  X = " << *(X->level()) << std::endl; 
    auto a=Matrix::dense(v_tot.size()+1, v_tot.size()+1, X->level());
 // // // // std::cout<< a->toString()<<std::endl;


 // // // //  std::cout<< Ms->toString()<<std::endl;
 
  // for(int i=0; i<v_tot.size()+1; i++)
  //      {
  //      for(int j=0; j<v_tot.size()+1;j++)
  //   {
  //     std::cout<< Ms->get(i,j)<<"  ";
     
  //   }
  //      std::cout<< std::endl;
  //   }



 for(int i=0; i<v_tot.size()+1;i++)
   {
      for(int j=0; j<v_tot.size()+1;j++)
   {
     std::cout<<  std::fixed << std::setprecision(2)<<a->get(i,j)<<"  ";
     
   }
      std::cout<< std::endl;
   } 

      
  return 0;
}
