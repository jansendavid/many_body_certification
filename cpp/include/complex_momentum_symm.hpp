#pragma once
#include "fusion.h"
#include"spins.hpp"
#include<unordered_map>
#include <memory>
using namespace mosek::fusion;
using namespace monty;
const double pi = std::acos(-1.0);
using basis_structure =std::map<int, std::vector<op_vec>>;
using TI_map_type =std::map<std::string, std::pair<std::string, std::complex<double>>>;

struct mom_ref
{
    mom_ref(Variable::t var, int i1, op_vec vec):var_(var),i1_(i1),vec_(vec){}
  Variable::t var_;
   int i1_;
  op_vec vec_;
  mom_ref()=default;
};
op_vec translation_y(op_vec op, int j, int L){
  op_vec vec;
  for(int i=0; i<op.size(); i++)
    {
      vec.push_back(op[i].get_translated_y(j,L));
    }
  return vec;
}
std::vector<op_vec> generate_all_translations_y(op_vec op, int L, int inc=1)
{
  std::vector<op_vec> all_T;
  all_T.push_back(op);
  for(int i=1; i<L; i+=inc)
    {
      auto new_op=translation_y(op, i,L);
      //      auto coeff=bubbleSort(new_op, new_op.size());
      auto [fac, vec] =get_normal_form(new_op);
      if(std::abs(fac.imag())>1e-9 or std::abs(fac.real()-1)>1e-9 ){std::cout<< "error in all trans"<<std::endl;}
      all_T.push_back(vec);
	
    }
  return all_T;
}


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
  if(print_op(op)=="1"){return all_T;}
  

  for(int i=1; i<L; i++)
    {
      auto new_op=translation(op, i,L);
      //    do we need auto
      auto [fac, vec] =get_normal_form(new_op);
      if(std::abs(fac.imag())>1e-9 or std::abs(fac.real()-1)>1e-9 ){std::cout<< "error in all trans"<<std::endl;}
      all_T.push_back(vec);
	
    }
  return all_T;
}
  void display(char a[], int n) 
{ 
  for (int i = 0; i < n; i++) { 
    std::cout << a[i] << " "; 
  } 
  std::cout << std::endl; 
}
 void findPermutations(char a[], int n) 
{ 
 
  // Sort the given array 
  std::sort(a, a + n); 
 
  // Find all possible permutations 
  std::cout << "Possible permutations are:\n"; 
  do { 
    display(a, n); 
  } while (std::next_permutation(a, a + n)); 
} 
std::vector<op_vec> generate_all_permutations_xy(op_vec op)
{
 
 
 std::vector<op_vec> all_P;
 all_P.push_back(op);


  std::vector<std::map<std::string, std::string>> permutations(1);
 permutations[0].insert({"x","y"});
 permutations[0].insert({"y","x"});
 permutations[0].insert({"z","z"});

 for(auto& a: permutations)
   {
     auto new_op=op;
     std::for_each(new_op.begin(), new_op.end(), [a](spin_op &n){
       auto old_d=n.dir_;
       n.dir_=a.at(old_d); });
     all_P.push_back(new_op);
   }
 

  return all_P;
}
std::vector<op_vec> generate_all_permutations_xyz(op_vec op)
{
 
 
 std::vector<op_vec> all_P;
 all_P.push_back(op);
 std::vector<std::map<std::string, std::string>> permutations(5);
 permutations[0].insert({"x","x"});
 permutations[0].insert({"y","z"});
 permutations[0].insert({"z","y"});

 permutations[1].insert({"x","y"});
 permutations[1].insert({"y","x"});
 permutations[1].insert({"z","z"});
 
 permutations[2].insert({"x","y"});
 permutations[2].insert({"y","z"});
 permutations[2].insert({"z","x"});

 permutations[3].insert({"x","z"});
 permutations[3].insert({"y","x"});
 permutations[3].insert({"z","y"});

 permutations[4].insert({"x","z"});
 permutations[4].insert({"y","y"});
 permutations[4].insert({"z","x"});


 for(auto& a: permutations)
   {
     auto new_op=op;
     std::for_each(new_op.begin(), new_op.end(), [a](spin_op &n){
       auto old_d=n.dir_;
       n.dir_=a.at(old_d); });
     all_P.push_back(new_op);
   }
 

  return all_P;
}
// desription:
//TI_map key: operator and value is pair, containing the key in total_refs, and the value to get it (is one for commuting variables)
// total_refs contains key: operator string, value mom_ref with variable reference, index, and op_vec


class momentum_block_base{
public:
  int L_;
  std::vector<op_vec> operators_;
  Eigen::MatrixXcd& FT_;
  Model::t M_;
  std::string sector_label_;
    std::string permuts_;
  std::vector<Variable::t> blocks_;
  std::map<std::string, Variable::t> G_variables_real;
  std::map<std::string, Variable::t> G_variables_imag;
  TI_map_type& TI_map_;
  std::map<std::string, mom_ref>& total_refs_;

  
  momentum_block_base(int L,std::vector<op_vec> operators, Model::t M, TI_map_type& TI_map,std::map<std::string, mom_ref>& total_refs,Eigen::MatrixXcd& FT,  std::string sector_label="",   std::string permuts="xyz"): L_(L), operators_(operators), M_(M),TI_map_(TI_map), total_refs_(total_refs), FT_(FT), sector_label_(sector_label), permuts_(permuts){
    if(permuts!="xyz" and permuts!="xy")
      {std::cout<< "permutation error"<<std::endl;}

      

  }

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
     auto [ti_key,ti_val]=TI_map_.at(print_op(vec));
     
       auto el=total_refs_.at(ti_key);
       cpx total_fac=fac;//fac*ti_val;
       M_->constraint( Expr::add(G_variables_real.at(key+"_real")->index(i),Expr::mul(-1.*total_fac.real(),el.var_)), Domain::equalsTo(0.0));
       
        M_->constraint( Expr::add(G_variables_imag.at(key+"_imag")->index(i),Expr::mul(-1.*total_fac.imag(),el.var_)), Domain::equalsTo(0.0));
       
     return;
    
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
		 
		 TI_map_.insert({print_op(op), { it->first,1.}});
		 found=true;
		 break;
	       }    
	   }
	   }
	   }
	 if(!found){mat_terms.insert({print_op(op), op });
	   TI_map_.insert({print_op(op), { print_op(op),1}});
	 }}
	    

	 return;
  }



};
class momentum_block: public momentum_block_base{
public:

  bool is_zero_{0};
  std::vector<int> block_shifts;
  momentum_block(int L,std::vector<op_vec> operators, Model::t M, bool is_zero, TI_map_type& TI_map,std::map<std::string, mom_ref>& total_refs, Eigen::MatrixXcd& FT, std::string sector_label="",  std::string permuts="xyz"):is_zero_(is_zero),  momentum_block_base(L, operators, M,TI_map, total_refs,FT, sector_label, permuts){
    
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
    M_->constraint( blocks_[0]->index(0,0), Domain::equalsTo(1.0));

    block_shifts.push_back(dim_0);
    // imaginary parts
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
	 M_->constraint(Expr::add(blocks_[0]->index(0,i+1),Expr::mul(-1.*fac.real()*std::sqrt(L_),el.var_)), Domain::equalsTo(0.0));
	 	 M_->constraint( Expr::add(blocks_[0]->index(i+1,0),Expr::mul(-1.,blocks_[0]->index(0,i+1))), Domain::equalsTo(0.0));
        // imaginary parts
	 M_->constraint( blocks_[0]->index(dim_0,i+1), Domain::equalsTo(0.0));
	 M_->constraint( blocks_[0]->index(dim_0+i+1,0), Domain::equalsTo(0.0));
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
       G_variables_real.insert({g_key+"_real", M_->variable(g_key+"_real",L_,Domain::inRange(-1., 1))});
       G_variables_imag.insert({g_key+"_imag", M_->variable(g_key+"_imag",L_,Domain::inRange(-1., 1))});
       
       //       std::cout<< g_key << "  " << g_key_2 << std::endl;
       //       std::cout<< print_op(*it1) << "  " << print_op(*it2) << std::endl;


			  for(int mat_pos=0; mat_pos<L_; mat_pos++)
			    {
		
			      auto expr_real=Expr::constTerm(0.);
			      auto expr_imag=Expr::constTerm(0.);

			      
			      // determines if first block of zeroth moment blocks
			      int shift=block_shifts[mat_pos]%operators_.size();
			      
			      // gives the shift between real and complex components
			      int dim=block_shifts[mat_pos];
			      
			      generate_G_elements(*it1, *it2, g_key, mat_pos);
			      
			      
			  for(int pos=0; pos<L_; pos++)
			    {
			      int position_in_G=pos;//?
			      //int position_in_G=(L_-pos)%L_;//?
  // 			      // maybe this position is not correct
  // 			      
			       double var_real=FT_(pos,mat_pos).real();
			       double var_imag=FT_(pos,mat_pos).imag();
			       expr_real=Expr::add(expr_real, Expr::mul(var_real,G_variables_real.at(g_key+"_real")->index(position_in_G) ));
			       expr_real=Expr::add(expr_real,  Expr::mul(-1.*var_imag,G_variables_imag.at(g_key+"_imag")->index(position_in_G))) ;

			       expr_imag=Expr::add(expr_imag,Expr::mul(var_imag,G_variables_real.at(g_key+"_real")->index(position_in_G) ) );
			       expr_imag=Expr::add(expr_imag,Expr::mul(var_real,G_variables_imag.at(g_key+"_imag")->index(position_in_G) ));

		
			    }
			  // real part
			  M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,j+shift),Expr::mul(-1.,expr_real)), Domain::equalsTo(0.0));	  
	        
			  // // imaginary part
				 
				M_->constraint( Expr::add(blocks_[mat_pos]->index(i+shift,dim+j+shift),Expr::mul(-1.,expr_imag)), Domain::equalsTo(0.0));
				if(i!=j)
				  {
			  // real part
			  M_->constraint( Expr::add(blocks_[mat_pos]->index(j+shift,i+shift),Expr::mul(-1.,expr_real)), Domain::equalsTo(0.0));	  
	        
			  // // imaginary part
				 
				M_->constraint( Expr::add(blocks_[mat_pos]->index(j+shift,dim+i+shift),Expr::mul(1.,expr_imag)), Domain::equalsTo(0.0));
				  }
			  	    }
			  
			  
			  j+=1;
			    }
	 i+=1;
       }

       for(int i=0; i<blocks_.size();i++)
	 {
	
	   int index_shift=block_shifts[i];

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



class momentum_basis{
  // note, the first sector must contain the unit element
public:
  int L_;
  basis_structure  operators_;
  Model::t M_;
  std::map<int, momentum_block> sectors_;
  std::string sector_;
   Variable::t XT_;
  std::vector<Variable::t> blocks_;
  std::map<std::string, std::pair<std::string, std::complex<double>>> TI_map_;
  std::map<std::string, mom_ref> total_refs_;
  Eigen::MatrixXcd FT_;
  
  momentum_basis(int L, basis_structure operators,Model::t M, std::string permuts="xyz"):L_(L),operators_(operators), M_(M)
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
	auto Block=momentum_block(L_,it->second,M_, block_zero,TI_map_,total_refs_, FT_, std::to_string(it->first), permuts);
	sectors_.insert({it->first, Block});

      }
     initialize_XT();
     std::cout<< "size TI map "<< TI_map_.size()<<std::endl;
     std::cout<< "size total refs "<< total_refs_.size()<<std::endl;
     // for(auto it=TI_map_.begin(); it!=TI_map_.end(); it++)
     //   {
     // 	 //	 std::cout<< it->first <<" --> "<<it->second.first<<std::endl;
     //   }
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
      //it->second.connect_Gelements();
     ++it_2;
       while(it_2!=sectors_.end())
	 {
	 it_2->second.generate_block();
	 //it->second.connect_Gelements();
	 ++it_2;
	 }
       // for(auto it=total_refs_.begin(); it!=total_refs_.end();++it)
       // 	 {
       // 	   auto n=TI_map_.find(it->first);
       // 	   if(n->first!=it->first){std::cout<< "dede"<<std::endl;}


       // 	 }
       // //       generate_TI_in_y();
      
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
  // void generate_TI_in_y()
  // {
  //   for(auto it=total_refs_.begin(); it!=total_refs_.end(); ++it)
  //     {
  // 	auto vec=it->second.vec_;
  // 	auto all=generate_all_translations_y(vec, L_);
  // 	for(auto a : all )
  // 	  {
  // 	    if(total_refs_.find(print_op(a))!=total_refs_.end() and (a!=vec))
  // 	      {
  // 		auto el1=total_refs_.at(print_op(a));
  // 		M_->constraint( Expr::add(el1.var_,Expr::mul(-1.,it->second.var_)), Domain::equalsTo(0.0));
  // 	      }
  // 	  }
  //     }
    

  // }
   
};
