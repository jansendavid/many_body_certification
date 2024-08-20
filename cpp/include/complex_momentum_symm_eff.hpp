#pragma once
#include "fusion.h"
#include"spins.hpp"
#include"symmetries.hpp"
#include<unordered_map>
#include <memory>
#include"complex_momentum_symm.hpp"
using namespace mosek::fusion;
using namespace monty;

class momentum_block_eff: public momentum_block_base{
public:

  bool is_zero_{0};
  std::vector<int> block_shifts;
  momentum_block_eff(int L,std::vector<op_vec> operators, Model::t M, bool is_zero, TI_map_type& TI_map,std::map<std::string, mom_ref>& total_refs, Eigen::MatrixXcd& FT, std::string sector_label="",  std::string permuts="xyz"):is_zero_(is_zero),  momentum_block_base(L, operators, M,TI_map, total_refs,FT, sector_label, permuts){
    
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

    void initialize_blocks()
  {
      
    int dim_x=operators_.size(); //dimension of other blocks

    // fix one
    if(total_refs_.find(print_op({}))!=total_refs_.end())
      {
    auto el=total_refs_.at(print_op({}));
    M_->constraint( el.var_, Domain::equalsTo(1.0));
      }
        // fix zero
    if(total_refs_.find("0")!=total_refs_.end())
      {
    auto el=total_refs_.at("0");
    M_->constraint( el.var_, Domain::equalsTo(0.0));
      }
      for(int i=0; i<L_; i++)
	{
	  std::string block_name="X"+std::to_string(i)+"_"+sector_label_;
	  blocks_.push_back(M_->variable(block_name, Domain::inPSDCone(2*dim_x)));
	  block_shifts.push_back(dim_x);
	}
       return;}
  void initialize_blocks_zero()
  {
      std::string block_name="X0_"+sector_label_;
    int dim_0=operators_.size()+1; //dimension of 0th block
    int dim_x=operators_.size(); //dimension of other blocks
    blocks_.push_back(M_->variable(block_name, Domain::inPSDCone(2*(dim_0))));
    M_->constraint( Expr::add(blocks_[0]->index(0,0),blocks_[0]->index(dim_0,dim_0)), Domain::equalsTo(1.0));

    block_shifts.push_back(dim_0);
    // imaginary parts, X^T_3[0,0]=X_3[0,0]
    M_->constraint( blocks_[0]->index(0,dim_0), Domain::equalsTo(0.0));

    // fix one
    if(total_refs_.find(print_op({}))!=total_refs_.end())
      {
    auto el=total_refs_.at(print_op({}));
    M_->constraint( el.var_, Domain::equalsTo(1.0));
      }
      for(int i=1; i<L_; i++)
	{
	  std::string block_name="X"+std::to_string(i)+"_"+sector_label_;
	  blocks_.push_back(M_->variable(block_name, Domain::inPSDCone(2*dim_x)));
	  block_shifts.push_back(dim_x);
	}
      // The "c" terms first row and column in block 0
      int i=0;
      for(auto it=operators_.begin(); it!=operators_.end(); ++it)
       {
	 auto op=*it;
	 auto [fac, vec] =get_normal_form(op);
	 auto ti_key=TI_map_.at(print_op(vec)).first;
	 auto el=total_refs_.at(ti_key);
	 if(std::abs(fac.imag())>1e-9){std::cout<< "errpr"<<std::endl;}
	 M_->constraint(Expr::add(Expr::add(blocks_[0]->index(0,i+1),blocks_[0]->index(dim_0,dim_0+i+1)),Expr::mul(-1.*fac.real()*std::sqrt(L_),el.var_)), Domain::equalsTo(0.0));

	 M_->constraint( Expr::add(Expr::add(blocks_[0]->index(i+1,0),blocks_[0]->index(dim_0+i+1,dim_0)),Expr::mul(-1.,blocks_[0]->index(0,i+1))), Domain::equalsTo(0.0));
	 //M_->constraint(Expr::add(blocks_[0]->index(0,i+1),Expr::mul(-1.*fac.real()*std::sqrt(L_),el.var_)), Domain::equalsTo(0.0));
	 // M_->constraint( Expr::add(blocks_[0]->index(i+1,0),Expr::mul(-1.,blocks_[0]->index(0,i+1))), Domain::equalsTo(0.0));
        // imaginary parts
	 M_->constraint( Expr::add(blocks_[0]->index(dim_0,i+1),Expr::mul(-1.,blocks_[0]->index(dim_0+i+1,0))), Domain::equalsTo(0.0));
	 
	 //M_->constraint( blocks_[0]->index(dim_0,i+1), Domain::equalsTo(0.0));
	 //M_->constraint( blocks_[0]->index(dim_0+i+1,0), Domain::equalsTo(0.0));
	 i++;
        } 
       return;}
   void  generate_block(){
  
     

    int i=0;
    for(auto it1=operators_.begin(); it1!=operators_.end(); ++it1)
       {
	 int j=i;
	 for(auto it2=it1; it2!=operators_.end(); ++it2)
     {
      
       std::string g_key="G_"+sector_label_+"_"+std::to_string(i)+"/"+std::to_string(j);
       std::string g_key_2="G_"+sector_label_+"_"+std::to_string(j)+"/"+std::to_string(i);



			  for(int mat_pos=0; mat_pos<L_; mat_pos++)
			    {
		
			      auto expr_real=Expr::constTerm(0.);
			      auto expr_imag=Expr::constTerm(0.);

			      
			      // determines if first block of zeroth moment blocks
			      int shift=block_shifts[mat_pos]%operators_.size();
			      
			      // gives the shift between real and complex components
			      int dim=block_shifts[mat_pos];
			      

			      for(int pos=0; pos<L_; pos++)
			    {
			      int position_in_G=pos;//?
			      auto construct= generate_single_G_element(*it1, *it2, g_key, pos);
			      //int position_in_G=(L_-pos)%L_;//?
  // 			      // maybe this position is not correct
  // 			      
			       double var_real=FT_(pos,mat_pos).real();
			       double var_imag=FT_(pos,mat_pos).imag();


			       
			        expr_real=Expr::add(expr_real, Expr::mul(var_real*construct.first.prefac_,construct.first.var_ ));
			       expr_real=Expr::add(expr_real,  Expr::mul(-1.*var_imag*construct.second.prefac_,construct.second.var_)) ;

			       expr_imag=Expr::add(expr_imag,Expr::mul(var_imag*construct.first.prefac_,construct.first.var_ ) );
			        expr_imag=Expr::add(expr_imag,Expr::mul(var_real*construct.second.prefac_,construct.second.var_ ));


		
			    }
			  // real part
			      M_->constraint( Expr::add( Expr::add(blocks_[mat_pos]->index(i+shift,j+shift), blocks_[mat_pos]->index(i+shift+dim,j+shift+dim)),Expr::mul(-1.,expr_real)), Domain::equalsTo(0.0));	  
			      //	  M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,j+shift),Expr::mul(-1.,expr_real)), Domain::equalsTo(0.0));	  
	        
			  // // imaginary part
			      M_->constraint( Expr::add(Expr::add(blocks_[mat_pos]->index(i+shift,dim+j+shift),Expr::mul(-1.,blocks_[mat_pos]->index(j+shift,dim+i+shift))),Expr::mul(-1.,expr_imag)), Domain::equalsTo(0.0));
			      //M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,dim+j+shift),Expr::mul(-1.,expr_imag)), Domain::equalsTo(0.0));

			  	    }
			  
			  
			  j+=1;
			    }
	 i+=1;
       }

       // for(int i=0; i<blocks_.size();i++)
       // 	 {
	
       // 	   int index_shift=block_shifts[i];

       // 	   for(int l=0; l<index_shift;l++)
       // 	     {
       // 	       for(int m=l;m<index_shift; m++)
       // 		 {
       // 		   //
       // 		   		   		   M_->constraint( Expr::add(blocks_[i]->index(l,m),Expr::mul(-1.,blocks_[i]->index(l+index_shift,m+index_shift))), Domain::equalsTo(0.0));
       // 						   M_->constraint( Expr::add(blocks_[i]->index(l,m+index_shift),Expr::mul(1.,blocks_[i]->index(m,l+index_shift))), Domain::equalsTo(0.0));
		    

       // 		 }
       // 	     }

       // 	 }

      return;}

};



class momentum_basis_eff{
  // note, the first sector must contain the unit element
public:
  int L_;
  basis_structure  operators_;
  Model::t M_;
  std::map<int, momentum_block_eff> sectors_;
  std::string sector_;
   Variable::t XT_;
  std::vector<Variable::t> blocks_;
  std::map<std::string, std::pair<std::string, std::complex<double>>> TI_map_;
  std::map<std::string, mom_ref> total_refs_;
  Eigen::MatrixXcd FT_;
  
  momentum_basis_eff(int L, basis_structure operators,Model::t M, std::string permuts="xyz"):L_(L),operators_(operators), M_(M)
  {
        FT_=Eigen::MatrixXcd(L,L);
    for(int i=0; i<L; i++)
      {
    for(int j=0; j<L; j++)
      {
	std::complex<double> phase(0.,-2.*i*j*pi/L_);
	
       FT_(i,j)=std::exp(phase);
      }
      }

    for(auto it=operators.begin(); it!=operators.end(); ++it)
      {
	bool block_zero=false;
	if(it->first==0)
	  {
	    
	    block_zero=true;
	  }
	else{
	    block_zero=false;
	}
	auto Block=momentum_block_eff(L_,it->second,M_, block_zero,TI_map_,total_refs_, FT_, std::to_string(it->first), permuts);
	sectors_.insert({it->first, Block});

      }
     initialize_XT();
     std::cout<< "size TI map "<< TI_map_.size()<<std::endl;
     std::cout<< "size total refs "<< total_refs_.size()<<std::endl;

     auto it=sectors_.begin();
     it->second.initialize_blocks_zero();
     ++it;
     while(it!=sectors_.end())
	 {
	 it->second.initialize_blocks();
	 ++it;
	 }
     
     auto it_2=sectors_.begin();
     it_2->second.generate_block();

     ++it_2;
       while(it_2!=sectors_.end())
	 {
	 it_2->second.generate_block();

	 ++it_2;
	 }

        return;
 
  }; 
  void initialize_XT()
  {
    std::map<std::string, op_vec> mat_terms;
   for(auto& b : sectors_)
       {

	 b.second.generate_TI_map(mat_terms);

       }
   int new_index=0;
   XT_=M_->variable("T", mat_terms.size(),Domain::inRange(-1., 1));
   for(auto a: mat_terms)
     {
       
       total_refs_.insert({a.first, mom_ref(XT_->index(new_index),new_index, a.second)});
	    
       
       new_index+=1;
     }
  }

};
