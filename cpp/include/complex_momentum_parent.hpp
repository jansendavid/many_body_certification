#pragma once
#include "fusion.h"
#include"spins.hpp"
#include<unordered_map>
#include <memory>
#include"symmetries.hpp"
using namespace mosek::fusion;
using namespace monty;
using basis_structure =std::map<int, std::vector<op_vec>>;
using int_pair=std::pair<int,int>;
using TI_map_type =std::map<std::string, std::pair<std::string, std::complex<double>>>;

class momentum_block_base{
public:
  int L_;
  std::vector<op_vec> operators_;
  Eigen::MatrixXcd& FT_;
  Model::t M_;
  std::string sector_label_;
    std::string permuts_;
  std::vector<std::vector<Variable::t>> blocks_variables_;
  // constants are stored
  std::vector<std::vector<Expression::t>> first_blocks_; // for each symmetri sector stores matrices and corresponding variables
  TI_map_type& TI_map_;
  std::map<std::string, int>& total_refs_;

  
  momentum_block_base(int L,std::vector<op_vec> operators, Model::t M, TI_map_type& TI_map,std::map<std::string, int>& total_refs,Eigen::MatrixXcd& FT,  std::string sector_label="",   std::string permuts="xyz"): L_(L), operators_(operators), M_(M),TI_map_(TI_map), total_refs_(total_refs), FT_(FT), sector_label_(sector_label), permuts_(permuts){
    if(permuts!="xyz" and permuts!="xy")
      {std::cout<< "permutation error"<<std::endl;}

      

  }
struct G_el{
    double prefac_;
    int pos_;
    G_el(double prefac, int pos ) : prefac_(prefac),pos_(pos) {};
            
  };
  // generate G element
    // convention, total_refs contains all elements appearing with prefact 1
  std::pair<G_el,G_el> generate_single_G_element(op_vec op1, op_vec op2, std::string key, int i)
  {


    auto op_dagger=dagger_operator(op1);
    op_vec new_op;

    if(i>0)
      {
	new_op=translation(op2, i, L_);
      }
    else{
      new_op=op2;
    }
    auto v_x=op_dagger;
    
    v_x.insert(v_x.end(), new_op.begin(),new_op.end());
     auto [fac, vec] =get_normal_form(v_x);
     auto [ti_key,ti_val]=TI_map_.at(print_op(vec));
     
       auto el=total_refs_.at(ti_key);
       cpx total_fac=fac;//fac*ti_val;

       auto real_val=G_el(total_fac.real(), el);
       
       //        M_->constraint( Expr::add(G_variables_imag.at(key+"_imag")->index(i),Expr::mul(-1.*total_fac.imag(),el.var_)), Domain::equalsTo(0.0));
	auto imag_val=G_el(total_fac.imag(), el);
	std::pair<G_el,G_el> construct(real_val, imag_val);
     return construct;
    
  }

  





    void check_if_operator_exists(op_vec op,std::map<std::string, op_vec>& mat_terms){
      // check if the operator is contained in functions
       if(op.size()<1)
	{
	  
       	  auto it=mat_terms.find(print_op(op));
       	  if(it != mat_terms.end())
       	       {
     		
       	       }
	  else{ mat_terms.insert({print_op(op), op });
       		 TI_map_.insert({print_op(op), { print_op(op),1}});}
      // 	  return;
       	}
       	  else
	    {

      auto all_t=generate_all_translations(op, L_);
      bool found=false;
	 for(auto op_t: all_t)
	   {
	     
	     
	     // generate all y translations, inc=1
	     auto all_ty=generate_all_translations_y(op_t, L_,1);
	     for(auto op_ty: all_ty)
	       {
	     // generate all permutations
		 std::vector<op_vec> all_p;	 
		 if(permuts_=="xyz")
		   {
	      all_p=generate_all_permutations_xyz(op_ty);
		   }
		 else if(permuts_=="xy")
		   {
		     all_p=generate_all_permutations_xy(op_ty);
		   }
		   
	    	 for(auto op_p: all_p)
	   { 
	     
	       auto it=mat_terms.find(print_op(op_p));
	       
	     if(it != mat_terms.end())
	       {
		 bool iszero_signsym_var=is_zero_signsym(op_p);
		 if(iszero_signsym_var)
		   {
		 TI_map_.insert({print_op(op), {"0",1.}});
		 found=true;
		 break;
		   }
		   else{
		 TI_map_.insert({print_op(op), { it->first,1.}});
		 found=true;
		 break;
		  }
	       }
	      auto op_mirror=mirror(op_p);
	       it=mat_terms.find(print_op(op_mirror));
	       if(it != mat_terms.end())
	       {
		 bool iszero_signsym_var=is_zero_signsym(op_mirror);
		 if(iszero_signsym_var)
		   {
		 TI_map_.insert({print_op(op), {"0",1.}});
		 found=true;
		 break;
		   }
		 else{
		 
		 TI_map_.insert({print_op(op), { it->first,1.}});
		 found=true;
		 break;
		  }
	       }
	   }
	   }
	   }
	 if(!found){

		 bool iszero_signsym_var=is_zero_signsym(op);
	    //	    std::cout<< print_op(op) << " and is zero "<< iszero<<std::endl;
	   if(iszero_signsym_var)
	     {
	       
	   mat_terms.insert({"0", op });
	   TI_map_.insert({print_op(op), { "0",1}});
	     }
	    else{
	     mat_terms.insert({print_op(op), op });
	   TI_map_.insert({print_op(op), { print_op(op),1}});
	   }

	 }}
	    

	 return;
  }



};

class momentum_block_child: public momentum_block_base{
public:

  bool is_zero_{0};
  std::vector<int> block_shifts;
  
  momentum_block_child(int L,std::vector<op_vec> operators, Model::t M, bool is_zero, TI_map_type& TI_map,std::map<std::string, int>& total_refs, Eigen::MatrixXcd& FT, std::string sector_label="",  std::string permuts="xyz"):is_zero_(is_zero),  momentum_block_base(L, operators, M,TI_map, total_refs,FT, sector_label, permuts){
    
  }
    void generate_TI_map(std::map<std::string, op_vec>& mat_terms){
    // generate all elemenets

     int index=0;
     for(auto it1=operators_.begin(); it1!=operators_.end(); ++it1)
       {
	 auto [fac, vec] =get_normal_form(*it1);
	 if(is_zero_)
	   {
	 check_if_operator_exists(vec, mat_terms);
	   }
	   else{
		
		TI_map_.insert({print_op(vec), { "0",1}});
	   }

    for(auto it2=it1; it2!=operators_.end(); ++it2)
      {


	auto op_cp=*it2;
	for(int n=0; n<L_;n++ )
	  {

	    op_vec new_op;
	    
	    if(n>0)
	      {
	
		new_op=translation(op_cp, n, L_);
	      }
	    else{
	      new_op=op_cp;}

     auto v_x=dagger_operator(*it1);
     v_x.insert(v_x.end(), new_op.begin(),new_op.end());


     auto [fac_tot, vec_tot] =get_normal_form(v_x);
     

     check_if_operator_exists(vec_tot, mat_terms);


     	  }

     
     
       }
   }
  
	return;
  }
};
