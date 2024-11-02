#pragma once
#include "fusion.h"
#include"spins.hpp"
#include<unordered_map>
#include <memory>
#include"symmetries.hpp"
#include"complex_momentum_parent.hpp"
#include <chrono>
using namespace mosek::fusion;
using namespace monty;
const double pi = std::acos(-1.0);


// desription:
//TI_map key: operator and value is pair, containing the key in total_refs, and the value to get it (is one for commuting variables)
// total_refs contains key: operator string, value mom_ref with variable reference, index, and op_vec



struct matrix_organizer{
  std::vector<int_pair> matrix_positions;
  std::vector<double> matrix_values;
  
  std::vector<double> b;
  Variable::t variables;
  int variable_index{0};
  
  void add_values(int_pair position,double value)
  {

    auto index=getIndex(matrix_positions, position);
    
    
    // Check if the target value was found
    if(index<0)
      {
	matrix_positions.push_back(position);
	matrix_values.push_back(value);
	
      }
    else
      {
	matrix_values[index]+=value;
      }
    
  }
  Matrix::t make_matrix(int dim1, int dim2)
  {
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> T;
    int i=0;
    for(auto& p:matrix_positions )
      {
	rows.push_back(p.first);
	cols.push_back(p.second);
	T.push_back(matrix_values[i]);

	i++;
      }
    
  
  return Matrix::sparse(dim1, dim2, nint(rows), nint(cols), ndou(T));
  }
};

class momentum_block_dual: public momentum_block_child{
public:
  momentum_block_dual(int L,std::vector<op_vec> operators, Model::t M, bool is_zero, TI_map_type& TI_map,std::map<std::string, int>& total_refs, Eigen::MatrixXcd& FT, std::string sector_label="",  std::string permuts="xyz"):  momentum_block_child(L, operators, M,is_zero, TI_map, total_refs,FT, sector_label, permuts){
    
  }
  void initialize_blocks()
  {
      
     int dim_x=operators_.size(); //dimension of other blocks

    // // fix one
       for(int i=0; i<L_; i++)
     	{
    
	  first_blocks_.push_back({});
     	  block_shifts.push_back(dim_x);
     	}
       return;
  }
  void initialize_blocks_zero( matrix_organizer& mat_org, std::vector<Expression::t>& second_block)
  {
        std::string block_name="X0_"+sector_label_;
	int dim_0=operators_.size()+1; //dimension of 0th block
	int dim_x=operators_.size(); //dimension of other blocks
	
      	std::vector<double> val={1.,1.};
     	std::vector<int> row={0, dim_0};
     	std::vector<int> col={0, dim_0};
     	Matrix::t Alpha  = (Matrix::sparse(2*(dim_0),2*(dim_0), nint(row), nint(col), ndou(val)));
	first_blocks_.push_back({});
	first_blocks_[0].push_back(Expr::mul(mat_org.variables->index(mat_org.variable_index),Alpha));
	mat_org.b[ mat_org.variable_index]+=1;
	mat_org.variable_index+=1;

	block_shifts.push_back(dim_0);


 // //    // // // fix one
     
      	{

	  auto el=total_refs_.at(print_op({}));
	  
	  mat_org.add_values( int_pair(el,mat_org.variable_index), 1.);
	  mat_org.b[ mat_org.variable_index]+=2;
	 
	    second_block[el]=Expr::add(second_block[el], mat_org.variables->index(mat_org.variable_index));
	  mat_org.variable_index+=1;

	}

     
     
 // // //     // fix zero
 // //       //     // fix zero
     if(total_refs_.find("0")!=total_refs_.end())
       {
	 
	 auto el=total_refs_.at("0");
  
	  mat_org.add_values( int_pair(el,mat_org.variable_index), 1.);
	  mat_org.b[ mat_org.variable_index]+=1;
	  //second_block[el]=Expr::add(second_block[el],mat_org.variables->index(mat_org.variable_index));
	  mat_org.variable_index+=1;
	 
       }
       for(int i=1; i<L_; i++)
	{

	  first_blocks_.push_back({});
      	   block_shifts.push_back(dim_x);
	   	}
 // //      // The "c" terms first row and column in block 0
       int i=0;
      
       std::cout<< "op size "<<operators_.size()<<std::endl;

       
       for(auto it=operators_.begin(); it!=operators_.end(); ++it)
       {
	 auto op=*it;
	 auto [fac, vec] =get_normal_form(op);
	 auto ti_key=TI_map_.at(print_op(vec)).first;
	 auto el=total_refs_.at(ti_key);

	 //	 if(std::abs(fac.imag())>1e-9){std::cout<< "errpr"<<std::endl;}
	 	 {
		    
		   std::vector<double> val={1./2,1./2,1./2,1./2};
		   std::vector<int> row={0, i+1,dim_0,i+1+dim_0 };
	     std::vector<int> col={i+1,0,i+1+dim_0,dim_0};
        Matrix::t Beta_i  = Matrix::sparse(2*dim_0,2*dim_0, nint(row), nint(col), ndou(val));
	
 first_blocks_[0].push_back(Expr::mul(mat_org.variables->index(mat_org.variable_index),Beta_i));
  mat_org.add_values( int_pair(el,mat_org.variable_index),-1.*fac.real()*std::sqrt(L_));
  mat_org.b[ mat_org.variable_index]+= -1.*fac.real()*std::sqrt(L_);

  
  mat_org.variable_index+=1;
	        }
		 	 
		 {
      // Why is not Tr[gamma M]=0?
   if(std::abs(fac.imag())>1e-9){std::cout<< "error c elements"<<std::endl;}
         std::vector<double> val={1./2,1./2,-1./2, -1./2};
	 std::vector<int> row={dim_0, i+1,dim_0,i+1+dim_0 };
	 std::vector<int> col={i+1,dim_0,i+1+dim_0,dim_0};
	 Matrix::t Beta_i  = Matrix::sparse(2*dim_0,2*dim_0, nint(row), nint(col), ndou(val));
first_blocks_[0].push_back(Expr::mul(mat_org.variables->index(mat_org.variable_index),Beta_i));
 mat_org.add_values( int_pair(el,mat_org.variable_index),  -1.*fac.imag());
 

 mat_org.b[ mat_org.variable_index]+= -1.*fac.imag();

 mat_org.variable_index+=1;
	

		 
	 }
     	 i++;
       }
       
         return;}
 
  void  generate_block(matrix_organizer& mat_org,std::vector<Expression::t>& second_block)
  {
    int i=0;
    

    const auto start{std::chrono::steady_clock::now()};
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
			      
			 
// 			      //here we store the extra term since <A_iX>=b_i, since the diagonal elements are, e.g.,  (1+sigma)
 			      double sum_real{0};
			      double sum_imag{0};
// 			 
			      // introduing two new variables

			      for(int pos=0; pos<L_; pos++)
			    {
			      int position_in_G=pos;//?
			      auto construct= generate_single_G_element(*it1, *it2, g_key, pos);
   			      
			       double var_real=FT_(pos,mat_pos).real();
			       double var_imag=FT_(pos,mat_pos).imag();
 			    
			       mat_org.add_values( int_pair(construct.first.pos_,mat_org.variable_index),(-1)*var_real*construct.first.prefac_);
			
			       mat_org.add_values( int_pair(construct.second.pos_,mat_org.variable_index),(-1.)*(-1)*var_imag*construct.second.prefac_);
			
 			       sum_real+=(-1)*var_real*construct.first.prefac_;


	
 			       sum_real+=(-1)*(-1.)*var_imag*construct.second.prefac_;


			       mat_org.add_values( int_pair( construct.first.pos_,mat_org.variable_index+1),(-1)*var_imag*construct.first.prefac_);
	  mat_org.add_values( int_pair( construct.second.pos_,mat_org.variable_index+1),(-1)*var_real*construct.second.prefac_);		       
// 	     		    
   
  	     sum_imag+=(-1)*var_imag*construct.first.prefac_;


  	     sum_imag+=(-1)*var_real*construct.second.prefac_;

   			    }
	
			      mat_org.b[mat_org.variable_index]+=(sum_real);
		 mat_org.b[mat_org.variable_index+1]+=(sum_imag);


			    

  

		  
  
	   {
	  std::vector<double> val={1./2,1./2, 1./2, 1./2};
	  std::vector<int> row={i+shift, j+shift,i+shift+dim,j+shift+dim };
	  std::vector<int> col={j+shift,i+shift,j+shift+dim,i+shift+dim};
        Matrix::t Beta_i  = Matrix::sparse(2*dim,2*dim, nint(row), nint(col), ndou(val));

	first_blocks_[mat_pos].push_back(Expr::mul(Beta_i,mat_org.variables->index(mat_org.variable_index)));
	
	}
	 
 	{
	  std::vector<double> val={1./2,1./2, -1./2, -1./2};
	 std::vector<int> row={i+shift,dim+j+shift,j+shift,dim+i+shift, };
	 std::vector<int> col={dim+j+shift,i+shift,dim+i+shift,j+shift};
	 Matrix::t Beta_i  = Matrix::sparse(2*dim,2*dim, nint(row), nint(col), ndou(val));
 
	 Matrix::t Beta_i_t  = Matrix::sparse(2*dim,2*dim, nint(col), nint(row), ndou(val));
	 first_blocks_[mat_pos].push_back(Expr::mul(Beta_i,mat_org.variables->index(mat_org.variable_index+1)));
	 
	}
 
 
		 mat_org.variable_index+=2;  
	
	

	

   			  	    }
			  
			  
			  j+=1;
			    }
     
	 i+=1;
     
       }

     const auto end{std::chrono::steady_clock::now()};
     const std::chrono::duration<double> elapsed_seconds{end - start};
    std::cout<< "time "<< elapsed_seconds.count()<<std::endl;
       return;}

  void set_matrices_in_PSDcone()
  {
    
    for(auto& block :first_blocks_ )
      {
        \
        if(block.size()>0)
        {
        M_->constraint(Expr::neg(Expr::add(new_array_ptr(block))), Domain::inPSDCone());
        }
    
       }
           \
    return;
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
  std::vector<Expression::t> second_block_;
  matrix_organizer mat_org_;
  std::vector<double> b_;
  Matrix::t C_;
  
 
  
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

      second_block_=std::vector<Expression::t>(total_refs_.size(),Expr::constTerm(0.0));
      
     int number_of_variables=3+2*operators_[0].size();
     // each symmetrie sector has L sub sectors with complex and real variable and 0.5*op.size()*(op.size()+1) constrains
     for(auto op: operators_)
       {
	 number_of_variables+=2*(0.5*op.second.size()*(op.second.size()+1))*L_;
	

       }
     std::cout<< "variable size "<<number_of_variables<<std::endl;
     mat_org_.variables=M->variable("variables",number_of_variables);
      mat_org_.b=std::vector<double>(number_of_variables, 0);

      auto it=sectors_.begin();
      it->second.initialize_blocks_zero(mat_org_, second_block_);
      ++it;
     while(it!=sectors_.end())
	 {
	   it->second.initialize_blocks();
	 ++it;
	 }
   
       	for(auto it=sectors_.begin(); it!=sectors_.end(); ++it)
	  {
	    it->second.generate_block(mat_org_, second_block_);
   
	      
	   
	 }

     

        return;
 
  };
   void prep_SDP()
  {
    auto total_matrix =mat_org_.make_matrix(total_refs_.size(),mat_org_.b.size());

    auto tot=Expr::mul(total_matrix,  mat_org_.variables);
    
  
       for(int i=0; i<total_refs_.size(); i++)
	 {
	   M_->constraint(tot->index(i),  Domain::lessThan(C_->get(i,0)));

    }
      	for(auto it=sectors_.begin(); it!=sectors_.end(); ++it)
	  {
	     it->second.set_matrices_in_PSDcone();
	  }
  }
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
   
    auto a = monty::new_array_ptr<double>(mat_org_.b);
    return Expr::dot(a, mat_org_.variables);
  }
  void set_C(Matrix::t& C)
  {
    C_=C;


  }
};
