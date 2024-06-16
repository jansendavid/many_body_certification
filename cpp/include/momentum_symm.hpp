#pragma once
#include "fusion.h"
#include"electrons.hpp"
#include<unordered_map>
#include <memory>
#include"functions.hpp"
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
  std::unordered_map<std::string, Variable::t> d_vector_map;
  std::unordered_map<std::string, Variable::t> G_variables;
  std::unordered_map<std::string, std::pair<std::string, int>> TI_map;

  std::unordered_map<std::string, mom_ref> total_refs_;
  // convention, total_refs contains all elements appearing with prefact 1
  
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
  void relate_elements()
  {
    // for(auto it=total_refs_.begin();it!=total_refs_.end();++it )
    //   {	
	
    // 	for(int j=1; j<L_; j++)
    // 	  {
    // 	    auto string=it->second.vec_;
    // 	    auto new_op=translation(string, j,L_);
    // 	    auto coeff=bubbleSort(new_op, new_op.size());
    // 	   auto it3=it;
    // 	   ++it3;
    // 	   for(auto it2=it3;it2!=total_refs_.end(); ++it2 )
    // 	     {
    // 	       if(print_op(it2->second.vec_)== print_op(new_op))
    // 		 {
    // 		   //std::cout<< "found trans "<<print_op(it->second.vec_)<< " to "<< print_op(new_op)<< " coeff "<< coeff<<std::endl;
    // 		   M_->constraint( Expr::add(it->second.var_,Expr::mul(-1*coeff,it2->second.var_ ) ), Domain::equalsTo(0.0));
    // 		   break;
    // 		 }
    // 	     }
    // 	  }
	
    //   }
      return; }
  void generate_G_elements(op_vec op1, op_vec op2, std::string key, int i)
  {
    //    std::cout<< "generate G elements "<<key<<std::endl;

    auto op_dagger=dagger_operator(op1);
    op_vec new_op;
    //     std::cout<< "left "<< print_op(op_dagger)<<std::endl;
    //std::cout<< "right "<< print_op(op2)<<std::endl;
    if(i>0)
      {
	new_op=translation(op2, i, L_);
      }
    else{
      new_op=op2;
    }
    auto v_x=op_dagger;

    v_x.insert(v_x.end(), new_op.begin(),new_op.end());
    int coeff= bubbleSort(v_x, v_x.size());
    bool dummy_variable=false;
    auto all_terms=generate_all_terms(v_x,dummy_variable);
    //    std::cout<< "term "<< print_op(v_x)<<std::endl;
     std::vector<Expression::t> expressions={};
	    for(int l=0; l<all_terms.size(); l++)
	      {
		auto [ key,coeff_new]=TI_map.at(print_op(all_terms[l].second));
	   	auto el=total_refs_.at(key);
	   	expressions.push_back(Expr::mul(-1.*all_terms[l].first*coeff_new,el.var_));
		//std::cout<< "expression "<< print_op(all_terms[l].second) << "  coeff "<<-1.*all_terms[l].first<<std::endl;
	   }	    
	    auto expressions_pt=monty::new_array_ptr<Expression::t>(expressions);
	    Expression::t exp=Expr::add(expressions_pt);
	    if(dummy_variable){
	      //std::cout<< "dummy "<<std::endl;
	      exp=Expr::add(-1., exp);}
	    //std::cout<< "coeff "<<coeff<<std::endl;
	    M_->constraint( Expr::add(Expr::mul(coeff,G_variables.at(key)->index(i)),exp), Domain::equalsTo(0.0));
    return;	
  }

  
};


class momentum_zero_block: public momentum_block_base{
public:
  //  std::unordered_map<int, >;
 
 

  momentum_zero_block(int L,std::vector<op_vec> operators, Model::t M, std::string sector=""):  momentum_block_base(L, operators, M, sector){
  //  void initialize_problem();
    generate_XT();
    relate_elements();
    generate_block();
    
  }
  void check_if_operator_exists(op_vec op,std::unordered_map<std::string, op_vec>& mat_terms){
    auto all_t=generate_all_translations(op, L_);
	 bool found=false;
	 //std::cout<< print_op(op)<<std::endl;
	 //std::cout<< "start "<<std::endl;
	 for(auto op_t: all_t)
	   {
	     int coeff_new= bubbleSort(op_t, op_t.size());
	       auto it=mat_terms.find(print_op(op_t));
	       //std::cout<<print_op(op_t)<<std::endl;
	     if(it != mat_terms.end())
	       {
		 
		 TI_map.insert({print_op(op), { it->first,coeff_new}});
		 found=true;
		 break;
	       }    
	   }
	 if(!found){mat_terms.insert({print_op(op), op });
	   TI_map.insert({print_op(op), { print_op(op),1}});
	 }
	 //std::cout<< "end "<<std::endl;
	 return;
  }
  void generate_XT(){
    // generate all elemenets
     std::unordered_map<std::string, op_vec> mat_terms;
     std::unordered_map<std::string, op_vec> extra_terms;
     
     int index=0;
     for(int i=0; i<operators_.size(); i++)
       {
	 check_if_operator_exists(operators_[i], mat_terms);
	 
   // 	  index+=1;
    for(int j=0; j<operators_.size(); j++)
      {
        auto op_dagger=dagger_operator(operators_[j]);
        int coeff= bubbleSort(op_dagger, op_dagger.size());
       
        check_if_operator_exists(op_dagger, mat_terms);
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
   // 	    //	    std::cout<< "start n="<<n <<" "<< print_op(op_dagger)<<print_op(op_cp)<<std::endl;
     auto v_x=op_dagger;
     v_x.insert(v_x.end(), new_op.begin(),new_op.end());
     int coeff= bubbleSort(v_x, v_x.size());
     bool dummy_variable=false;
     //std::cout<< print_op(v_x)<<std::endl;
     auto all_terms=generate_all_terms(v_x,dummy_variable);
     //     std::cout<< "size "<< all_terms.size()<<std::endl;
    for(auto op: all_terms)
      {
        check_if_operator_exists(op.second, mat_terms);
	  }


    	  }

     
     
      }
  }
     //  std::cout<< "size "<< mat_terms.size()<< " and "<<TI_map.size()<<std::endl;
     // for(auto a : mat_terms)
     //   {std::cout<< a.first<<std::endl;}
     // std::cout<<std::endl;
     // for(auto v:TI_map )
     //   {
     // 	 std::cout<<v.first<< "  "<< v.second.first << " "<< v.second.second<<std::endl; 
     //   }
	 int new_index=0;
	 XT_=M_->variable("T", mat_terms.size(),Domain::inRange(-1., 1));
	 for(auto a: mat_terms)
	   {
	     total_refs_.insert({a.first, mom_ref(XT_->index(new_index),new_index, a.second)});
	     if(check_zero(a.second))
	       {
		 
		  M_->constraint( XT_->index(new_index), Domain::equalsTo(0.0));
	       }
	     
	     new_index+=1;
	   }
       
	return;
  }
  void  generate_block(){
    std::string block_name="X0";
    int dim_0=operators_.size()+1; //dimension of 0th block
    int dim_x=operators_.size(); //dimension of other blocks
    blocks_.push_back(M_->variable(block_name, Domain::inPSDCone(2*(dim_0))));
    M_->constraint( blocks_[0]->index(0,0), Domain::equalsTo(1.0));
    M_->constraint( blocks_[0]->index(dim_0,dim_0), Domain::equalsTo(1.0));
    // imaginary parts
    M_->constraint( blocks_[0]->index(0,dim_0), Domain::equalsTo(0.0));
    M_->constraint( blocks_[0]->index(dim_0,0), Domain::equalsTo(0.0));
      for(int i=1; i<L_; i++)
	{
	  std::string block_name="X"+std::to_string(i);
	  blocks_.push_back(M_->variable(block_name, Domain::inPSDCone(2*dim_x)));
	}
      
       for(int i=0; i<operators_.size(); i++)
       {
	 auto op=operators_[i];
       auto op_dagger=dagger_operator(operators_[i]);
       int coeff= bubbleSort(op_dagger, op_dagger.size());
       auto [ key,coeff_new]=TI_map.at(print_op(op));
       auto el=total_refs_.at(print_op(op));
       
       // std::cout<< print_op(op)<< " i "<< i<<print_op(el.vec_)<< "inedc "<< el.i1_ <<std::endl;
       //std::cout<< print_op(op_dagger)<<" i "<<i<<std::endl;
	 // convention, in XT the terms like c_ic_j are stored
	 // the matrix element is prefac*XT element
	 //     std::cout<< "i= "<<0<< ","<<i+1<< "  op "<< print_op(op)<<std::endl;
	      
       M_->constraint(Expr::add(blocks_[0]->index(0,i+1),Expr::mul(-1.*1*std::sqrt(L_)*coeff_new,el.var_)), Domain::equalsTo(0.0));
	 M_->constraint( Expr::add(blocks_[0]->index(dim_0,dim_0+i+1),Expr::mul(-1.,blocks_[0]->index(0,i+1))), Domain::equalsTo(0.0));
       // imaginary parts
	 M_->constraint( blocks_[0]->index(0,dim_0+i+1), Domain::equalsTo(0.0));
	 M_->constraint( blocks_[0]->index(dim_0+i+1,0), Domain::equalsTo(0.0));
	 M_->constraint( blocks_[0]->index(dim_0,i+1), Domain::equalsTo(0.0));
	 M_->constraint( blocks_[0]->index(i+1,dim_0), Domain::equalsTo(0.0));
	 auto [ key_dag,coeff_new_dag]=TI_map.at(print_op(op_dagger));
	 el=total_refs_.at(key_dag);
	 //	      std::cout<< "i= "<<i+1<< ","<<0<< "  op "<< print_op(op_dagger)<<std::endl;
	 // convention, in XT the terms like c_ic_j are stored
	 // the matrix element is prefac*XT element
	 M_->constraint( Expr::add(blocks_[0]->index(i+1,0),Expr::mul(-1.*1*std::sqrt(L_)*coeff_new_dag,el.var_)), Domain::equalsTo(0.0));
	 M_->constraint( Expr::add(blocks_[0]->index(dim_0+i+1,dim_0),Expr::mul(-1.,blocks_[0]->index(i+1,0) )), Domain::equalsTo(0.0));
	 }

       
  // // generate all elements
   //      std::cout<<FT<<std::endl;
    for(int i=0; i<operators_.size(); i++)
       {
      
   for(int j=i; j<operators_.size(); j++)
     {
       std::string g_key="G_"+sector_+"_"+std::to_string(i)+"/"+std::to_string(j);
       std::string g_key_2="G_"+sector_+"_"+std::to_string(j)+"/"+std::to_string(i);
       G_variables.insert({g_key, M_->variable(g_key,L_,Domain::inRange(-1., 1))});
       std::vector<Expression::t> expressions_1_real(L_, Expr::constTerm(0.));
       std::vector<Expression::t> expressions_1_imag(L_, Expr::constTerm(0.));
       std::vector<Expression::t> expressions_2_real(L_,Expr::constTerm(0.));
       std::vector<Expression::t> expressions_2_imag(L_,Expr::constTerm(0.));
			  if(i!=j)
			    {
			             G_variables.insert({g_key_2, M_->variable(g_key_2,L_,Domain::inRange(-1., 1))});

			    }
			  for(int mat_pos=0; mat_pos<L_; mat_pos++)
			    {
			      int shift=1;
			      int dim=dim_0;
			      if(mat_pos>0){shift=0; dim=dim_x;}
			      generate_G_elements(operators_[i], operators_[j], g_key, mat_pos);
			      if(i!=j)
				{ generate_G_elements(operators_[j], operators_[i], g_key_2, mat_pos);}
			      
			  for(int pos=0; pos<L_; pos++)
			    {
			      int position_in_G=pos;//?
			      // maybe this position is not correct
			      //			      std::cout<< "po "<<pos<<"  "<<FT(pos,mat_pos)<< std::endl;
			      double var_real=FT(mat_pos,pos).real();
				double var_imag=FT(mat_pos,pos).imag();
			  
			      expressions_1_real[mat_pos]=Expr::add(expressions_1_real[mat_pos], Expr::mul(var_real,G_variables.at(g_key)->index(position_in_G) ));
			      expressions_1_imag[mat_pos]=Expr::add(expressions_1_imag[mat_pos], Expr::mul(var_imag,G_variables.at(g_key)->index(position_in_G) ));
			      if(i!=j)
				{
				 expressions_2_real[mat_pos]=Expr::add(expressions_2_real[mat_pos], Expr::mul(var_real,G_variables.at(g_key_2)->index(position_in_G) ));
			      expressions_2_imag[mat_pos]=Expr::add(expressions_2_imag[mat_pos], Expr::mul(var_imag,G_variables.at(g_key_2)->index(position_in_G) ));
				}
			    }
			  auto expressions_pt_1_real=monty::new_array_ptr<Expression::t>(expressions_1_real);
			  auto expressions_pt_1_imag=monty::new_array_ptr<Expression::t>(expressions_1_imag);
			  Expression::t exp_1_real=Expr::add(expressions_pt_1_real);
			  Expression::t exp_1_imag=Expr::add(expressions_pt_1_imag);
			  
			  if(i==j)
			    {
			      
			   M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,j+shift),Expr::mul(-1.,exp_1_real)), Domain::equalsTo(0.0));
			  M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,j+shift),Expr::mul(-1.,blocks_[mat_pos]->index(dim+i+shift,dim+j+shift))), Domain::equalsTo(0.0));

			  M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,dim+j+shift),Expr::mul(-1.,exp_1_imag)), Domain::equalsTo(0.0));
			  M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,dim+j+shift),Expr::mul(1.,blocks_[mat_pos]->index(dim+j+shift,i+shift))), Domain::equalsTo(0.0));
			  // diagonal part of H_I is zero
			  M_->constraint(blocks_[mat_pos]->index(i+shift,dim+j+shift), Domain::equalsTo(0.0));
			  M_->constraint(blocks_[mat_pos]->index(dim+j+shift,i+shift), Domain::equalsTo(0.0));




			  // if(mat_pos==1)
			  //   {
			  //     std::cout<< "real mat nr, "<< mat_pos<< " shape: "<< 2*dim <<"x"<< 2*dim<<std::endl;
			  //     std::cout<< i+shift<<","<<j+shift<< " = "<<dim+i+shift <<" , "<<dim+j+shift<<std::endl;
			  //     std::cout<< "complec mat nr, "<< mat_pos<< " shape: "<< 2*dim <<"x"<< 2*dim<<std::endl;
			  //     std::cout<<"fix zero "<< i+shift<<","<<dim+j+shift<< " and plus "<<dim+j+shift <<" , "<<i+shift<<std::endl;}
			     }
			  
			  
			    else{
			      // real part
			      auto expressions_pt_2_real=monty::new_array_ptr<Expression::t>(expressions_2_real);
			      auto expressions_pt_2_imag=monty::new_array_ptr<Expression::t>(expressions_2_imag);
			      Expression::t exp_2_real=Expr::add(expressions_pt_2_real);
			      Expression::t exp_2_imag=Expr::add(expressions_pt_2_imag);
			      M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,j+shift),Expr::mul(-1.,exp_1_real)), Domain::equalsTo(0.0));
			      M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,j+shift),Expr::mul(-1.,blocks_[mat_pos]->index(dim+i+shift,dim+j+shift))), Domain::equalsTo(0.0));
			      M_->constraint( Expr::add(blocks_[mat_pos]->index(j+shift,i+shift),Expr::mul(-1.,blocks_[mat_pos]->index(dim+j+shift,dim+i+shift))), Domain::equalsTo(0.0));
			     
			      M_->constraint( Expr::add(blocks_[mat_pos]->index(j+shift,i+shift),Expr::mul(-1.,blocks_[mat_pos]->index(i+shift,j+shift))), Domain::equalsTo(0.0));			  

			      //imaginary part 			  
    
			      M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,dim+j+shift),Expr::mul(-1.,exp_1_imag)), Domain::equalsTo(0.0));
			      M_->constraint( Expr::add(blocks_[mat_pos]->index(dim+j+shift,i+shift),Expr::mul(-1.,exp_2_imag)), Domain::equalsTo(0.0));
			      M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,dim+j+shift),Expr::mul(-1.,blocks_[mat_pos]->index(dim+j+shift,i+shift))), Domain::equalsTo(0.0));
			      M_->constraint( Expr::add(blocks_[mat_pos]->index(j+shift,dim+i+shift),Expr::mul(-1.,blocks_[mat_pos]->index(dim+i+shift,j+shift))), Domain::equalsTo(0.0));
			      M_->constraint( Expr::add(blocks_[mat_pos]->index(j+shift,dim+i+shift),Expr::mul(1.,blocks_[mat_pos]->index(i+shift,dim+j+shift))), Domain::equalsTo(0.0));
			  // if(mat_pos==1)
			  //   {
			  //     //  std::cout<< "real mat nr, "<< mat_pos<< " shape: "<< 2*dim <<"x"<< 2*dim<<std::endl;
			  //     // std::cout<< i+shift<<","<<j+shift<< " = "<<dim+i+shift <<" , "<<dim+j+shift<<std::endl;
			  //     // std::cout<<"real "<< j+shift<<","<<i+shift<< " and "<<dim+j+shift <<" , "<<dim+i+shift<<std::endl;
			  //     std::cout<< "complec mat nr, "<< mat_pos<< " shape: "<< 2*dim <<"x"<< 2*dim<<std::endl;
			  //    std::cout<< " "<<i+shift<<","<<dim+j+shift<< " =  "<< dim+j+shift<<" , "<<i+shift<<std::endl;		      
			      
			  //     std::cout<< " "<<j+shift<<","<<dim+i+shift<< " =  "<< dim+i+shift<<" , "<<j+shift<<std::endl;

			  //     //	     std::cout<< "set equal "<< j+shift<<","<<i+shift<< " and "<<i+shift <<" , "<<j+shift<<std::endl;
			  //    std::cout<< "set und equal "<< j+shift<<","<<dim+i+shift<< " and "<<i+shift<<","<<dim+j+shift<<std::endl;
					    
			  //   }


    			    }
			  if(mat_pos==1)
			    {
			      //    std::cout<< "end"<<std::endl;
			    }
									
				
			    }
			  

			    }
       }
  
      
      return;}

  
};
