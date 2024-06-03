#pragma once
#include "fusion.h"
#include"spins.hpp"
#include<unordered_map>
#include <memory>
using namespace mosek::fusion;
using namespace monty;
const double pi = std::acos(-1.0);
using basis_structure =std::unordered_map<int, std::vector<op_vec>>;

struct mom_ref
{
    mom_ref(Variable::t var, int i1, op_vec vec):var_(var),i1_(i1),vec_(vec){}
  Variable::t var_;
   int i1_;
  op_vec vec_;
  mom_ref()=default;
};
op_vec translation(op_vec op, int j, int L){
  op_vec vec;
  for(int i=0; i<op.size(); i++)
    {
      vec.push_back(op[i].get_translated(j,L));
    }
  return vec;
}
std::vector<op_vec> generate_all_translations(op_vec op, int L)
{
  std::vector<op_vec> all_T;
  all_T.push_back(op);
  for(int i=1; i<L; i++)
    {
      auto new_op=translation(op, i,L);
      //      auto coeff=bubbleSort(new_op, new_op.size());
      all_T.push_back(new_op);
	
    }
  return all_T;
}

// desription:
//TI_map key: operator and value is pair, containing the key in total_refs, and the value to get it (is one for commuting variables)
// total_refs contains key: operator string, value mom_ref with variable reference, index, and op_vec
class momentum_basis{
public:
  int L_;
  basis_structure  operators_;
  momentum_basis(int L, basis_structure operators):L_(L),operators_(operators){}; 
};

class momentum_block_base{
public:
  int L_;
  std::vector<op_vec> operators_;
  Eigen::MatrixXcd FT;
  Model::t M_;
  std::string sector_;
   Variable::t XT_;
    std::vector<Variable::t> blocks_;
  std::map<std::string, Variable::t> d_vector_map;
  std::map<std::string, Variable::t> G_variables_real;
  std::map<std::string, Variable::t> G_variables_imag;
  std::map<std::string, std::pair<std::string, std::complex<double>>> TI_map;

  std::map<std::string, mom_ref> total_refs_;
  // convention, total_refs contains all elements appearing with prefact 1
 void generate_G_elements(op_vec op1, op_vec op2, std::string key, int i)
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

     // if(vec.size()<1)
     //   {
     // 	  M_->constraint( G_variables_real.at(key+"_real")->index(i), Domain::equalsTo(1.0*fac.real()));
     // 	 M_->constraint( G_variables_imag.at(key+"_imag")->index(i), Domain::equalsTo(1.*fac.imag()));
     // 	 //std::cout<< print_op(v_x)<<std::endl;
     // 	 //std::cout<< fac.real() << " "<<fac.imag()<<std::endl;
     //   }
     // else{

       auto ti_key=TI_map.at(print_op(vec)).first;
       auto ti_val=TI_map.at(print_op(vec)).second;
       auto el=total_refs_.at(ti_key);
       cpx total_fac=fac*ti_val;
       M_->constraint( Expr::add(G_variables_real.at(key+"_real")->index(i),Expr::mul(-1.*total_fac.real(),el.var_)), Domain::equalsTo(0.0));
     M_->constraint( Expr::add(G_variables_imag.at(key+"_imag")->index(i),Expr::mul(-1.*total_fac.imag(),el.var_)), Domain::equalsTo(0.0));
     // }
     
     return;
    
  }

  




  momentum_block_base(int L,std::vector<op_vec> operators, Model::t M, std::string sector=""): L_(L), operators_(operators), M_(M), sector_(sector){
    FT=Eigen::MatrixXcd(L,L);
    for(int i=0; i<L; i++)
      {
    for(int j=0; j<L; j++)
      {
	std::complex<double> phase(0.,-2.*i*j*pi/L);
	
       FT(i,j)=std::exp(phase);
      }
      }


  }
   
    void check_if_operator_exists(op_vec op,std::map<std::string, op_vec>& mat_terms){
       if(op.size()<1)
      	{
	  
       	  auto it=mat_terms.find(print_op(op));
       	  if(it != mat_terms.end())
       	       {
     		
       	       }
	  else{ mat_terms.insert({print_op(op), op });
       		 TI_map.insert({print_op(op), { print_op(op),1}});}
      // 	  return;
       	}
       	  else
	    {
      auto all_t=generate_all_translations(op, L_);
	 bool found=false;
	 //std::cout<< print_op(op)<<std::endl;
	 //std::cout<< "start "<<std::endl;
	 for(auto op_t: all_t)
	   {
	     //	     int coeff_new= sort(op_t, op_t.size());
	     	auto [fac, op_nf] =get_normal_form(op_t);
	       auto it=mat_terms.find(print_op(op_nf));
	       //std::cout<<print_op(op_t)<<std::endl;
	     if(it != mat_terms.end())
	       {
		 
		 TI_map.insert({print_op(op), { it->first,fac}});
		 found=true;
		 break;
	       }    
	   }
	 if(!found){mat_terms.insert({print_op(op), op });
	   TI_map.insert({print_op(op), { print_op(op),1}});
	 }}
	 //std::cout<< "end "<<std::endl;
	 return;
  }
  void generate_XT(){
    // generate all elemenets
     std::map<std::string, op_vec> mat_terms;
     std::map<std::string, op_vec> extra_terms;
     
     int index=0;
     for(int i=0; i<operators_.size(); i++)
       {
	 auto [fac, vec] =get_normal_form(operators_[i]);
	 check_if_operator_exists(vec, mat_terms);
	 
    for(int j=0; j<operators_.size(); j++)
      {
        auto op_dagger=dagger_operator(operators_[j]);
        
	auto [fac, vec_dagger] =get_normal_form(op_dagger);
	check_if_operator_exists(vec_dagger, mat_terms);
	auto op_cp=operators_[i];
	for(int n=0; n<L_;n++ )
	  {
	    op_vec new_op;
	    
	    if(n>0)
	      {
	
		new_op=translation(op_cp, n, L_);
	      }
	    else{
	      new_op=op_cp;}
	    //	    std::cout<< "start n="<<n <<" "<< print_op(op_dagger)<<print_op(op_cp)<<std::endl;
     auto v_x=op_dagger;
     v_x.insert(v_x.end(), new_op.begin(),new_op.end());
     auto [fac_tot, vec_tot] =get_normal_form(v_x);
     check_if_operator_exists(vec_tot, mat_terms);

     	  }

     
     
       }
   }
     //  std::cout<< "size "<< mat_terms.size()<< " and "<<TI_map.size()<<std::endl;
     // for(auto a : mat_terms)
     // {std::cout<< a.first<<std::endl;}
      
      for(auto v:TI_map )
        {
	  //std::cout<<v.first<< "  "<< v.second.first << " "<< v.second.second<<std::endl; 
        }
      	 int new_index=0;
	 XT_=M_->variable("T", mat_terms.size(),Domain::inRange(-1., 1));
	 for(auto a: mat_terms)
	   {
	     //	     std::cout<<"F "<< a.first<< "  "<<std::endl;
	     total_refs_.insert({a.first, mom_ref(XT_->index(new_index),new_index, a.second)});
	     // if(check_zero(a.second))
	     //   {
		 
	     // 	  M_->constraint( XT_->index(new_index), Domain::equalsTo(0.0));
	     //   }
	     
	     new_index+=1;
	   }
       
	return;
  }


};
class momentum_zero_block: public momentum_block_base{
public:
  //  std::unordered_map<int, >;
 
 

  momentum_zero_block(int L,std::vector<op_vec> operators, Model::t M, std::string sector=""):  momentum_block_base(L, operators, M, sector){
  //  void initialize_problem();
    generate_XT();
    // relate_elements();
            generate_block();
    
  }
   void  generate_block(){
    std::string block_name="X0";
    int dim_0=operators_.size()+1; //dimension of 0th block
    int dim_x=operators_.size(); //dimension of other blocks
    blocks_.push_back(M_->variable(block_name, Domain::inPSDCone(2*(dim_0))));
    M_->constraint( blocks_[0]->index(0,0), Domain::equalsTo(1.0));
    // imaginary parts
    M_->constraint( blocks_[0]->index(0,dim_0), Domain::equalsTo(0.0));

    // fix one
    auto el=total_refs_.at(print_op({}));
    M_->constraint( el.var_, Domain::equalsTo(1.0));
      for(int i=1; i<L_; i++)
	{
	  std::string block_name="X"+std::to_string(i);
	  blocks_.push_back(M_->variable(block_name, Domain::inPSDCone(2*dim_x)));
	}
      // The "c" terms first row and column in block 0
       for(int i=0; i<operators_.size(); i++)
       {
	 auto op=operators_[i];
	 auto ti_key=TI_map.at(print_op(op)).first;
	 auto el=total_refs_.at(ti_key);
        M_->constraint(Expr::add(blocks_[0]->index(0,i+1),Expr::mul(-1.*1*std::sqrt(L_),el.var_)), Domain::equalsTo(0.0));
	M_->constraint( Expr::add(blocks_[0]->index(i+1,0),Expr::mul(-1.,blocks_[0]->index(0,i+1))), Domain::equalsTo(0.0));
        // imaginary parts
	 M_->constraint( blocks_[0]->index(0,dim_0+i+1), Domain::equalsTo(0.0));
	 M_->constraint( blocks_[0]->index(i+1,dim_0), Domain::equalsTo(0.0));
        }
       
  // // // generate all elements
  //  //      std::cout<<FT<<std::endl;
       for(int i=0; i<operators_.size(); i++)
       {
      
   for(int j=i; j<operators_.size(); j++)
     {
      
       std::string g_key="G_"+sector_+"_"+std::to_string(i)+"/"+std::to_string(j);
       std::string g_key_2="G_"+sector_+"_"+std::to_string(j)+"/"+std::to_string(i);
       G_variables_real.insert({g_key+"_real", M_->variable(g_key+"_real",L_,Domain::inRange(-1., 1))});
       G_variables_imag.insert({g_key+"_imag", M_->variable(g_key+"_imag",L_,Domain::inRange(-1., 1))});
       std::vector<Expression::t> expressions_1_real(L_, Expr::constTerm(0.));
       std::vector<Expression::t> expressions_1_imag(L_, Expr::constTerm(0.));
       std::vector<Expression::t> expressions_2_real(L_,Expr::constTerm(0.));
       std::vector<Expression::t> expressions_2_imag(L_,Expr::constTerm(0.));
			  if(i!=j)
			    {
			      G_variables_real.insert({g_key_2+"_real", M_->variable(g_key_2+"_real",L_,Domain::inRange(-1., 1))});
			      G_variables_imag.insert({g_key_2+"_imag", M_->variable(g_key_2+"_imag",L_,Domain::inRange(-1., 1))});
			     
			    }
			  for(int mat_pos=0; mat_pos<L_; mat_pos++)
			    {
			      auto expr_real=Expr::constTerm(0.);
			      auto expr_imag=Expr::constTerm(0.);
			      int shift=1;
			      int dim=dim_0;
			      if(mat_pos>0){shift=0; dim=dim_x;}
			      generate_G_elements(operators_[i], operators_[j], g_key, mat_pos);
			      if(i!=j)
				{ generate_G_elements(operators_[j], operators_[i], g_key_2, mat_pos);}
			      
			  for(int pos=0; pos<L_; pos++)
			    {
			      	      int position_in_G=pos;//?
			      // int position_in_G=(L_-pos)%L_;//?
  // 			      // maybe this position is not correct
  // 			      
   			      double var_real=FT(mat_pos,pos).real();
			      double var_imag=FT(mat_pos,pos).imag();
			      expr_real=Expr::add(expr_real, Expr::mul(var_real,G_variables_real.at(g_key+"_real")->index(position_in_G) ));
			      expr_real=Expr::add(expr_real, Expr::mul(var_real, Expr::mul(-1.*var_imag,G_variables_imag.at(g_key+"_imag")->index(position_in_G))) );
			      expr_imag=Expr::add(expr_imag,Expr::mul(var_imag,G_variables_real.at(g_key+"_real")->index(position_in_G) ) );
			      expr_imag=Expr::add(expr_imag,Expr::mul(var_real,G_variables_imag.at(g_key+"_imag")->index(position_in_G) ));
			
			    }
   			
			  
			   if(i==j)
   			    {
			      
   			   M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,j+shift),Expr::mul(-1.,expr_real)), Domain::equalsTo(0.0));	  
			  // imaginary part
   			  M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,dim+j+shift),Expr::mul(-1.,expr_imag)), Domain::equalsTo(0.0));

			   // fix diagonal part of H_I is zero
			   M_->constraint(blocks_[mat_pos]->index(i+shift,dim+j+shift), Domain::equalsTo(0.0));
   			  
		         

   			     }
			  

    			    else
   			       {
		
				 M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,j+shift),Expr::mul(-1.,expr_real)), Domain::equalsTo(0.0));
				 M_->constraint( Expr::add(blocks_[mat_pos]->index(j+shift,i+shift),Expr::mul(-1.,expr_real)), Domain::equalsTo(0.0));	  
			  // imaginary part
   			  M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,dim+j+shift),Expr::mul(-1.,expr_imag)), Domain::equalsTo(0.0));
			  M_->constraint( Expr::add(blocks_[mat_pos]->index(j+shift,dim+i+shift),Expr::mul(1.,expr_imag)), Domain::equalsTo(0.0));



      			    }
					
				
			  	    }
			  

			    }
       }
       for(int i=0; i<blocks_.size();i++)
	 {
	   int index_shift=0;
	   if(i==0)
	     {index_shift=dim_0;}
	   else{index_shift=dim_x;}
	   for(int l=0; l<index_shift;l++)
	     {
	       for(int m=l;m<index_shift; m++)
		 {
		   
		   M_->constraint( Expr::add(blocks_[i]->index(l,m),Expr::mul(-1.,blocks_[i]->index(l+index_shift,m+index_shift))), Domain::equalsTo(0.0));
		    M_->constraint( Expr::add(blocks_[i]->index(l,m+index_shift),Expr::mul(1.,blocks_[i]->index(m,l+index_shift))), Domain::equalsTo(0.0));

		 }
	     }

	 }
       
      return;}

};
