#pragma once
#include "fusion.h"
#include"spins.hpp"
#include<unordered_map>
#include <memory>
#include"symmetries.hpp"
using namespace mosek::fusion;
using namespace monty;
const double pi = std::acos(-1.0);
using basis_structure =std::map<int, std::vector<op_vec>>;
using TI_map_type =std::map<std::string, std::pair<std::string, std::complex<double>>>;
std::shared_ptr<ndarray<int,1>>    nint(const std::vector<int> &X)    { return new_array_ptr<int>(X); }
std::shared_ptr<ndarray<double,1>> ndou(const std::vector<double> &X) { return new_array_ptr<double>(X); }
int getIndex(std::vector<int> v, int K) 
{ 
    auto it = find(v.begin(), v.end(), K); 
  
    // If element was found 
    if (it != v.end())  
    { 
      
        // calculating the index 
        // of K 
        int index = it - v.begin(); 
        return index; 
    } 
    else { 
        // If the element is not 
        // present in the vector 
        return -1; 
    } 
} 
// desription:
//TI_map key: operator and value is pair, containing the key in total_refs, and the value to get it (is one for commuting variables)
// total_refs contains key: operator string, value mom_ref with variable reference, index, and op_vec


class momentum_block_base_dual{
public:
  int L_;
  std::vector<op_vec> operators_;
  Eigen::MatrixXcd& FT_;
  Model::t M_;
  std::string sector_label_;
    std::string permuts_;
  std::vector<std::vector<Variable::t>> blocks_variables_;
  // constants are stored
  std::vector<Expression::t> first_blocks_;
  TI_map_type& TI_map_;
  std::map<std::string, int>& total_refs_;

  
  momentum_block_base_dual(int L,std::vector<op_vec> operators, Model::t M, TI_map_type& TI_map,std::map<std::string, int>& total_refs,Eigen::MatrixXcd& FT,  std::string sector_label="",   std::string permuts="xyz"): L_(L), operators_(operators), M_(M),TI_map_(TI_map), total_refs_(total_refs), FT_(FT), sector_label_(sector_label), permuts_(permuts){
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
class momentum_block_dual: public momentum_block_base_dual{
public:

  bool is_zero_{0};
  std::vector<int> block_shifts;
  
  momentum_block_dual(int L,std::vector<op_vec> operators, Model::t M, bool is_zero, TI_map_type& TI_map,std::map<std::string, int>& total_refs, Eigen::MatrixXcd& FT, std::string sector_label="",  std::string permuts="xyz"):is_zero_(is_zero),  momentum_block_base_dual(L, operators, M,TI_map, total_refs,FT, sector_label, permuts){
    
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

  void initialize_blocks(Expression::t& second_block,std::vector<Variable::t>& variables_,std::vector<double>& b)
  {
      
     int dim_x=operators_.size(); //dimension of other blocks

    // // fix one
       for(int i=0; i<L_; i++)
     	{
    // 	  std::string block_name="X"+std::to_string(i)+"_"+sector_label_;
    // 	  blocks_.push_back(M_->variable(block_name, Domain::inPSDCone(2*dim_x)));
	  std::vector<double> val={};
     	std::vector<int> row={};
     	std::vector<int> col={};
	Matrix::t Alpha  = Matrix::dense(Matrix::sparse(2*(dim_x),2*(dim_x), nint(row), nint(col), ndou(val)));
     	first_blocks_.push_back(Expr::constTerm(Alpha));
     	  block_shifts.push_back(dim_x);
     	}
       return;
  }
  void initialize_blocks_zero( Expression::t& second_block,std::vector<Variable::t>& variables, std::vector<double>& b)
  {
        std::string block_name="X0_"+sector_label_;
      int dim_0=operators_.size()+1; //dimension of 0th block
      int dim_x=operators_.size(); //dimension of other blocks
      	std::vector<double> val={1.,1.};
     	std::vector<int> row={0, dim_0};
     	std::vector<int> col={0, dim_0};
     	Matrix::t Alpha  = Matrix::dense(Matrix::sparse(2*(dim_0),2*(dim_0), nint(row), nint(col), ndou(val)));
	// initialize variables
	
	variables.push_back(M_->variable("ones", 1));
	
     	first_blocks_.push_back(Expr::mul(variables.back()->index(0),Alpha));

	b.push_back(1);
	//	variables.push_back(variables_in_r_);
	block_shifts.push_back(dim_0);

 //    // // // fix one
     
      	{
	  auto el=total_refs_.at(print_op({}));
	  std::vector<double> val={1.};
     	std::vector<int> row={el};
     	std::vector<int> col={0};
     	Matrix::t Alpha  = Matrix::sparse(total_refs_.size(),1, nint(row), nint(col), ndou(val));
	
      	variables.push_back(M_->variable("ones_2", 1));
      b.push_back(2);
      second_block=Expr::add(second_block,Expr::mul(variables.back()->index(0),Alpha));

	}
 // //      //    M_->constraint( el.var_, Domain::equalsTo(1.0));
     
     
 // //     // fix zero
 // //       //     // fix zero
     if(total_refs_.find("0")!=total_refs_.end())
       {
	 
     auto el=total_refs_.at("0");

     std::vector<double> val={1.};
     	std::vector<int> row={el};
     	std::vector<int> col={0};
     	Matrix::t Alpha  = Matrix::sparse(total_refs_.size(),1, nint(row), nint(col), ndou(val));
      b.push_back(1);
      variables.push_back(M_->variable("zero", 1));
      second_block=Expr::add(second_block,Expr::mul(variables.back()->index(0), Alpha));

      // M_->constraint( el.var_, Domain::equalsTo(0.0));
       }
       for(int i=1; i<L_; i++)
     	{
    // 	  std::string block_name="X"+std::to_string(i)+"_"+sector_label_;
    // 	  blocks_.push_back(M_->variable(block_name, Domain::inPSDCone(2*dim_x)));
	  std::vector<double> val={};
     	std::vector<int> row={};
     	std::vector<int> col={};
	Matrix::t Alpha  = Matrix::dense(Matrix::sparse(2*(dim_x),2*(dim_x), nint(row), nint(col), ndou(val)));
     	first_blocks_.push_back(Expr::constTerm(Alpha));
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
	 	 {
		   std::vector<double> val={1./2,1./2,1./2,1./2};
		   std::vector<int> row={0, i+1,dim_0,i+1+dim_0 };
	     std::vector<int> col={i+1,0,i+1+dim_0,dim_0};
        Matrix::t Beta_i  = Matrix::sparse(2*dim_0,2*dim_0, nint(row), nint(col), ndou(val));
	 

	 std::vector<double> val_e={-1.*fac.real()*std::sqrt(L_)};
     	std::vector<int> row_e={el};
     	std::vector<int> col_e={0};
     	Matrix::t Alpha_e  = Matrix::sparse(total_refs_.size(),1, nint(row_e), nint(col_e), ndou(val_e));
	
	variables.push_back(M_->variable("C_r"+std::to_string(i), 1));
	 b.push_back(-1.*fac.real()*std::sqrt(L_));
	 second_block=Expr::add(second_block,Expr::mul(variables.back()->index(0),Alpha_e));

      
	 first_blocks_[0]=Expr::add(first_blocks_[0], Expr::mul(Beta_i,variables.back()->index(0)));

	
        }
		 

	 
 {
   if(std::abs(fac.imag())>1e-9){std::cout<< "error c elements"<<std::endl;}
         std::vector<double> val={1./2,1./2,-1./2, -1./2};
	 std::vector<int> row={dim_0, i+1,dim_0,i+1+dim_0 };
	 std::vector<int> col={i+1,dim_0,i+1+dim_0,dim_0};
	 Matrix::t Beta_i  = Matrix::sparse(2*dim_0,2*dim_0, nint(row), nint(col), ndou(val));
	  std::vector<double> val_e={-1.*fac.imag()};
     	std::vector<int> row_e={el};
     	std::vector<int> col_e={0};
     	Matrix::t Alpha_e  = Matrix::sparse(total_refs_.size(),1, nint(row_e), nint(col_e), ndou(val_e));
	

	 b.push_back(-1.*fac.imag());	
	 variables.push_back(M_->variable("C_i"+std::to_string(i), 1));

	 second_block=Expr::add(second_block,Expr::mul(variables.back()->index(0),Alpha_e));
	 first_blocks_[0]=Expr::add(first_blocks_[0], Expr::mul(Beta_i,variables.back()->index(0)));

	 
	 }
     	 i++;
       } 
       return;}
  class loop_object{
  public:
    // this object will after the loop create will be element of a vector
    //M_->variable("C_r"+std::to_string(mat_pos)+g_key, size); that it will append to variables_
    // Matrix is of size size and afterwards we add blokc[i]+=matrix[i]*variable[i] for i : mat_pos
    // afterwards we create <prefac, variable> that we add to the "diagonal index"
    std::vector<Matrix::t> matrix_vec_;
    int size_{0};
    std::vector<int> mat_pos_;
    std::vector<double> prefac_vec_;
    

  };
  void  generate_block(Expression::t& second_block,std::vector<Variable::t>& variables, std::vector<double>& b , std::vector<int>& count)
  {
     std::vector<loop_object>  OBJ(total_refs_.size());
    int i=0;
    std::cout<< "start "<<std::endl;
    auto temp=std::vector<Expression::t>(total_refs_.size(),Expr::constTerm(0.0));
    int XX=0;
    std::vector<Variable::t> variables_temp_real;
    std::vector<Variable::t> variables_temp_imag;
        std::vector<double> b_temp_real;
    std::vector<double> b_temp_imag;
    
    for(auto it1=operators_.begin(); it1!=operators_.end(); ++it1)
       {
	 int j=i;
	 for(auto it2=it1; it2!=operators_.end(); ++it2)
     {
      
       std::string g_key="G_"+sector_label_+"_"+std::to_string(i)+"/"+std::to_string(j);
       std::string g_key_2="G_"+sector_label_+"_"+std::to_string(j)+"/"+std::to_string(i);



			  for(int mat_pos=0; mat_pos<L_; mat_pos++)
			    {
		

			      
			      // determines if first block of zeroth moment blocks
			      int shift=block_shifts[mat_pos]%operators_.size();
			      
			      // gives the shift between real and complex components
			      int dim=block_shifts[mat_pos];
			      
			 
			      //here we store the extra term since <A_iX>=b_i, since the diagonal elements are, e.g.,  (1+sigma)
			      double sum_real{0};
			      double sum_imag{0};
			      variables_temp_real.push_back(M_->variable("C_r"+std::to_string(mat_pos)+g_key, 1));
			      variables_temp_imag.push_back(M_->variable("C_i"+std::to_string(mat_pos)+g_key, 1));
	std::vector<int> row_e_real;
	std::vector<double> val_e_real;
	std::vector<int> row_e_im;
	std::vector<double> val_e_im;	
			      for(int pos=0; pos<L_; pos++)
			    {
			      int position_in_G=pos;//?
			      auto construct= generate_single_G_element(*it1, *it2, g_key, pos);
  // 			      
			       double var_real=FT_(pos,mat_pos).real();
			       double var_imag=FT_(pos,mat_pos).imag();
			       auto index_r_1=getIndex(row_e_real,  construct.first.pos_);
			       auto index_r_2=getIndex(row_e_real,  construct.second.pos_);
			       if(index_r_1 <0)
				 {
				   row_e_real.push_back(construct.first.pos_);
				   val_e_real.push_back((-1)*var_real*construct.first.prefac_);
				 }
			       else{
				 val_e_real[index_r_1]+=(-1)*var_real*construct.first.prefac_;
			       }
			     if(index_r_2 <0)
				 {
				   row_e_real.push_back(construct.second.pos_);
				   val_e_real.push_back((-1.)*(-1)*var_imag*construct.second.prefac_);
				 }
			       else{
				 val_e_real[index_r_2]+=(-1.)*(-1)*var_imag*construct.second.prefac_;
			       }
			       //			       temp[construct.first.pos_]=Expr::add(temp[construct.first.pos_],Expr::mul(variables_temp_real.back()->index(0),(-1)*var_real*construct.first.prefac_));
	 sum_real+=(-1)*var_real*construct.first.prefac_;
	 //temp[construct.second.pos_]=Expr::add(temp[construct.second.pos_],Expr::mul(variables_temp_real.back()->index(0),(-1.)*(-1)*var_imag*construct.second.prefac_));

	
	     sum_real+=(-1)*(-1.)*var_imag*construct.second.prefac_;
	     
	      
			       auto index_i_1=getIndex(row_e_im,  construct.first.pos_);
			       auto index_i_2=getIndex(row_e_im,  construct.second.pos_);
			       if(index_i_1 <0)
				 {
				   row_e_im.push_back(construct.first.pos_);
				   val_e_im.push_back((-1)*var_imag*construct.first.prefac_);
				 }
			       else{
				 val_e_im[index_i_1]+=(-1)*var_imag*construct.first.prefac_;
			       }
			     if(index_i_2 <0)
				 {
				   row_e_im.push_back(construct.second.pos_);
				   val_e_im.push_back((-1)*var_real*construct.second.prefac_);
				 }
			       else{
				 val_e_im[index_i_2]+=(-1)*var_real*construct.second.prefac_;
			       }
			     
	     //temp[construct.first.pos_]=Expr::add(temp[construct.first.pos_],Expr::mul(variables_temp_imag.back()->index(0),(-1)*var_imag*construct.first.prefac_));
	     sum_imag+=(-1)*var_imag*construct.first.prefac_;

	     //     temp[construct.second.pos_]=Expr::add(temp[construct.second.pos_],Expr::mul(variables_temp_imag.back()->index(0),(-1)*var_real*construct.second.prefac_));
	     sum_imag+=(-1)*var_real*construct.second.prefac_;

 			    }
	
   
// // // // // 			      // note no multiplictaion with -1 in bs since was done in the above code
			    
// std::vector<double> val_e={-1.*fac.imag()};
  
		  {
		    std::vector<int> col_e(row_e_real.size(),0);
		    Matrix::t Alpha_e_1  = Matrix::sparse(total_refs_.size(),1, nint(row_e_real), nint(col_e), ndou(val_e_real));
		    second_block=Expr::add(second_block, Expr::mul(variables_temp_real.back()->index(0), Alpha_e_1));
		  }

		  {
		    std::vector<int> col_e(row_e_im.size(),0);
		  
		    Matrix::t Alpha_e_1  = Matrix::sparse(total_refs_.size(),1, nint(row_e_im), nint(col_e), ndou(val_e_im));
		    second_block=Expr::add(second_block, Expr::mul(variables_temp_imag.back()->index(0), Alpha_e_1));
		  }
		  
	{
	  	 
	   
	  std::vector<double> val={1./2,1./2};
	  std::vector<int> row={i+shift, i+shift+dim };
	  std::vector<int> col={j+shift,j+shift+dim};
        Matrix::t Beta_i  = Matrix::sparse(2*dim,2*dim, nint(row), nint(col), ndou(val));
	Matrix::t Beta_i_t=Matrix::sparse(2*dim,2*dim, nint(col), nint(row), ndou(val));
	first_blocks_[mat_pos]=Expr::add(first_blocks_[mat_pos], Expr::mul(Beta_i,variables_temp_real.back()->index(0)));
	first_blocks_[mat_pos]=Expr::add(first_blocks_[mat_pos], Expr::mul(Beta_i_t,variables_temp_real.back()->index(0)));
        
	}
	 
	{
	  std::vector<double> val={1./2,-1./2};
	 std::vector<int> row={i+shift,j+shift };
	 std::vector<int> col={dim+j+shift,dim+i+shift};
	 Matrix::t Beta_i  = Matrix::sparse(2*dim,2*dim, nint(row), nint(col), ndou(val));
	 Matrix::t Beta_i_t  = Matrix::sparse(2*dim,2*dim, nint(col), nint(row), ndou(val));
	 first_blocks_[mat_pos]=Expr::add(first_blocks_[mat_pos], Expr::mul(Beta_i,variables_temp_imag.back()->index(0)));
	 first_blocks_[mat_pos]=Expr::add(first_blocks_[mat_pos], Expr::mul(Beta_i_t,variables_temp_imag.back()->index(0)));
 	}

    	 b_temp_real.push_back(sum_real);

    	 b_temp_imag.push_back(sum_imag);	

 			  	    }
			  
			  
			  j+=1;
			    }
	 i+=1;
       }
    variables.insert (variables.end(), variables_temp_real.begin(), variables_temp_real.end());
    variables.insert (variables.end(), variables_temp_imag.begin(), variables_temp_imag.end());

    b.insert (b.end(), b_temp_real.begin(), b_temp_real.end());
    b.insert (b.end(), b_temp_imag.begin(),b_temp_imag.end());
    
  

      return;}

  void set_matrices_in_PSDcone()
  {
    for(auto& mat :first_blocks_)
      {
	M_->constraint(Expr::mul(-1.,mat), Domain::inPSDCone());
      }
  }

 };



class momentum_basis_dual{
  // note, the first sector must contain the unit element
public:
  int L_;
  basis_structure  operators_;
  Model::t M_;
  std::map<int, momentum_block_dual> sectors_;
  std::string sector_;
   Variable::t XT_;
  std::vector<Variable::t> blocks_;
  std::map<std::string, std::pair<std::string, std::complex<double>>> TI_map_;
  std::map<std::string, int> total_refs_; // reference to the the element
  Eigen::MatrixXcd FT_;
  Expression::t second_block_;
  //std::vector<Variable::t> variables_;
  std::vector<Variable::t> variables_;
  std::vector<double> b_;
  std::vector<int> count;
  
 
  
  momentum_basis_dual(int L, basis_structure operators,Model::t M, std::string permuts="xyz"):L_(L),operators_(operators), M_(M)
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
	auto Block=momentum_block_dual(L_,it->second,M_, block_zero,TI_map_,total_refs_, FT_, std::to_string(it->first), permuts);
	sectors_.insert({it->first, Block});

      }
      initialize_XT();
     std::cout<< "size TI map "<< TI_map_.size()<<std::endl;

     second_block_=Expr::constTerm(total_refs_.size(), 0.0);


      auto it=sectors_.begin();
      it->second.initialize_blocks_zero(second_block_, variables_, b_);
      ++it;
     while(it!=sectors_.end())
	 {
	   it->second.initialize_blocks(second_block_,variables_, b_);
	 ++it;
	 }
     
       auto it_2=sectors_.begin();
       it_2->second.generate_block(second_block_,variables_, b_, count);
       it_2->second.set_matrices_in_PSDcone();

      ++it_2;
         while(it_2!=sectors_.end())
	 {
	   it_2->second.generate_block( second_block_,variables_, b_, count);
		 it_2->second.set_matrices_in_PSDcone();

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
   for(auto a: mat_terms)
     {
       
       total_refs_.insert({a.first, new_index});
	    
       
       new_index+=1;
     }
  }
  
   Expression::t get_costfunction()
  {
	  auto h=Expr::constTerm(0.);
	  std::cout<< "b size "<<b_.size()<<std::endl;
	  for(int i=0; i<b_.size(); i++)
	 {
	   h=Expr::add(h, Expr::mul(b_[i], variables_[i]->index(0)));
	 

	 }
	 return h;
  }

};
