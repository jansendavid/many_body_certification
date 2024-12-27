#pragma once
#include "fusion.h"
#include"spins.hpp"
#include"symmetries.hpp"
#include<unordered_map>
#include <memory>
#include"util.hpp"
#include"complex_momentum_parent.hpp"
#include<cassert>
using namespace mosek::fusion;
using namespace monty;

using namespace mosek::fusion;
using namespace monty;
using symmetry_sector =std::map<int,std::vector<std::vector<matrix_organizer>>>;
// implementing momentum symmetrie in x and y direction
class momentum_block_eff: public momentum_block_child{
public:
  std::vector<std::vector<Variable::t>> blocks_;
  int  sign_sector_{0};
  std::vector<std::vector<int>> block_shifts;
  momentum_block_eff(int L,std::vector<op_vec> operators, Model::t M, int sign_sector, TI_map_type& TI_map,std::map<std::string, int>& total_refs, Eigen::MatrixXcd& FT, std::string sector_label="",  std::string permuts="xyz"):sign_sector_(sign_sector),  momentum_block_child(L, operators, M,false,TI_map, total_refs,FT, sector_label, permuts){
    
  }


    void initialize_blocks_general()
  {
      
    int dim_x=operators_.size(); //dimension of other blocks
 
  for(int j=0; j<L_; j++)
  {
  block_shifts.push_back({});

      for(int i=0; i<L_; i++)
	{
	  block_shifts[j].push_back(dim_x);
	}
  }

       return;}
void initialize_blocks(std::map<std::string,symmetry_sector>& As)
{
  if(sign_sector_==0)
  {initialize_blocks_zero(As);}
  else{
    initialize_blocks_general();
  }
  return;
}
  void initialize_blocks_zero(std::map<std::string,symmetry_sector>& As)
  {
//        std::string block_name="Xb0_"+sector_label_;
    int dim_0=operators_.size()+1; //dimension of 0th block
    int dim_x=operators_.size(); //dimension of other blocks
//     blocks_.push_back({});
block_shifts.push_back({});
//     blocks_[0].push_back(M_->variable(block_name, Domain::inPSDCone(2*(dim_0))));

// initializing block shifts
block_shifts[0].push_back(dim_0);
    for(int i=1; i<L_; i++)
    {

  block_shifts[0].push_back(dim_x);

    }

      for(int i=1; i<L_; i++)
	{
    block_shifts.push_back({});
    for(int j=0; j<L_; j++)
    {
	  block_shifts[i].push_back(dim_x);
	}
  }



	      As["1"][sign_sector_][0][0].add_values({0,0},1./2);
	      As["1"][sign_sector_][0][0].add_values({dim_0,dim_0},1./2);
  

//   //     // The "c" terms first row and column in block 0
       int i=0;
      
      for(auto it=operators_.begin(); it!=operators_.end(); ++it)
       {
	 auto op=*it;
   // get normal form
	 auto [coeff, nf] =get_normal_form(op);
  // get translation invariant representation 
	 auto ti_key=TI_map_.at(print_op(nf)).first;
	 auto el=total_refs_.at(ti_key);
 	 
   if(std::abs(coeff.real())>1e-9)
   {
 
    		As[ti_key][sign_sector_][0][0].add_values({0,i+1},1./2*coeff.real()*std::sqrt(L_));
		As[ti_key][sign_sector_][0][0].add_values({i+1,0},1./2*coeff.real()*std::sqrt(L_));
		 	As[ti_key][sign_sector_][0][0].add_values({dim_0,i+1+dim_0},1./2*coeff.real()*std::sqrt(L_));
		 As[ti_key][sign_sector_][0][0].add_values({i+1+dim_0,dim_0},1./2*coeff.real()*std::sqrt(L_));
    
   }
   assert(std::abs(coeff.imag())<1e-9);

	 i++;
        } 
     
       return;}
   void  generate_block(std::map<std::string, symmetry_sector>& As){
  
// //           std::cout<< "start "<<std::endl;

//     const auto start{std::chrono::steady_clock::now()};

    int i=0;
    for(auto it1=operators_.begin(); it1!=operators_.end(); ++it1)
       {
	 int j=i;
	 for(auto it2=it1; it2!=operators_.end(); ++it2)
     {

			  for(int mat_pos_x=0; mat_pos_x<L_; mat_pos_x++)
			    {
             for(int mat_pos_y=0; mat_pos_y<L_; mat_pos_y++)
			    {
		
    

			      
// 			      // determines if first block of zeroth moment blocks
			      int shift=block_shifts[mat_pos_x][mat_pos_y]%operators_.size();
			      
// 			      // gives the shift between real and complex components
			      int dim=block_shifts[mat_pos_x][mat_pos_y];
			      

            for(int pos_y=0; pos_y<L_; pos_y++)
            {	      
            	 std::complex<double> FT_factor_y=FT_(pos_y,mat_pos_y);

			      for(int pos_x=0; pos_x<L_; pos_x++)
			    {
            	    		       std::complex<double> FT_factor_x=FT_(pos_x,mat_pos_x);
			 
//              // to do, correct so that all terms appearing here appear in map
              auto construct= generate_single_G_element_sos(*it1, *it2, pos_y, pos_x);
              std::complex<double> total_prefactor=construct.prefac_*FT_factor_x*FT_factor_y;
              if(std::abs(total_prefactor.real())>1e-9)
              {

                	As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values({i+shift,j+shift},1./2*total_prefactor.real());
		            As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values({i+shift+dim, j+shift+dim},1./2*total_prefactor.real());
if(i!=j)
{
   	As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values({j+shift,i+shift},1./2*total_prefactor.real());
		As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values({j+shift+dim, i+shift+dim},1./2*total_prefactor.real());
}
              }
                  if(std::abs(total_prefactor.imag())>1e-9)
              {
               	As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values({i+shift,dim+j+shift},1./2*total_prefactor.imag());
		            As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values({ j+shift,dim+i+shift},-1./2*total_prefactor.imag());
if(i!=j)
{
 
   	As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values({j+shift,dim+i+shift},1./2*total_prefactor.imag());
	  As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values({ i+shift,dim+j+shift},-1./2*total_prefactor.imag());
}
              }


           }
  



           }


          }
          }
 			  	    
			  
			  
			  j+=1;
			    }
	 i+=1;
       }


      return;
     }

  };



class momentum_basis_xy{
  // note, the first sector must contain the unit element
  // solves min(by), with sum_i y_i A_i <<C
public:
  int L_;
  basis_structure  operators_;
  Model::t M_;
  std::map<int, momentum_block_eff> sectors_;
  std::string sector_;

  
  std::map<std::string, std::pair<std::string, std::complex<double>>> TI_map_;
  std::map<std::string, int> total_refs_;
  Eigen::MatrixXcd FT_;
   std::vector<double> b_;
   // contains the matrices As, for each sign symmetrye we have LxL blocks
  std::map<std::string, symmetry_sector> As_;
  
  momentum_basis_xy(int L, basis_structure operators,Model::t M, std::string permuts="xyz"):L_(L),operators_(operators), M_(M)
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

	auto Block=momentum_block_eff(L_,it->second,M_, it->first,TI_map_,total_refs_, FT_, std::to_string(it->first), permuts);
	sectors_.insert({it->first, Block});

      }

      initialize_XT();
     std::cout<< "size TI map "<< TI_map_.size()<<std::endl;
     std::cout<< "size total refs "<< total_refs_.size()<<std::endl;


      for(auto it =total_refs_.begin(); it!=total_refs_.end(); it++)
      {
        As_.insert({it->first, symmetry_sector()});
        for(auto it_sign_sector=sectors_.begin(); it_sign_sector!=sectors_.end(); ++it_sign_sector)
{
As_[it->first][it_sign_sector->first]={};
  
        for(int i=0; i<L; i++)
        {
          As_[it->first][it_sign_sector->first].push_back({});
          for(int j=0; j<L; j++)
          {
          As_[it->first][it_sign_sector->first][i].push_back(matrix_organizer());

          }
        }

        }
 }
 
     for(auto& sector : sectors_)
	 {
	 sector.second.initialize_blocks(As_);

	 }

   for(auto it_2=sectors_.begin(); it_2!=sectors_.end(); ++it_2)
   {
    it_2->second.generate_block(As_);

   }
   std::cout<< "finished making the As matrices"<<std::endl;
 
        return;
 
  }; 
  void initialize_XT()
  {
    std::map<std::string, op_vec> mat_terms;
   for(auto& b : sectors_)
       {

	 b.second.generate_TI_map_xy(mat_terms);
   std::cout<<"sizes "<< mat_terms.size()<< " adn "<<TI_map_.size()<<std::endl;
 


       }
   
  

   int new_index=0;
   for(auto a: mat_terms)
     {
       
      
       total_refs_.insert({a.first, new_index});
	    
       
       new_index+=1;
     }
  }
void set_b( std::vector<double> b )
{

  b_=b;
}

};
class momentum_symmetry_solver_dual: public momentum_basis_xy{
public:
   Variable::t y_;
momentum_symmetry_solver_dual(int L, basis_structure operators,Model::t M, std::string permuts="xyz"):momentum_basis_xy(L, operators,M, permuts)
{
y_=M_->variable("T", total_refs_.size());
     // fix 1
     	auto el=total_refs_.at("1");
	  M_->constraint( y_->index(el), Domain::equalsTo(1.0));
    // fix zero
    el=total_refs_.at("0");
   M_->constraint( y_->index(el), Domain::equalsTo(0.0));
}
void fix_constrains(){
 

    // iterate over sign sector
    for(auto sign_symm_sector :sectors_)
    {
      
      for(int i=0; i<L_; i++)
      {
       for(int j=0; j<L_; j++)
      {
        std::vector<Expression::t> matrices;

        for(auto op :total_refs_ )
        {
          if(op.first=="0"){continue;}
          if(As_[op.first][sign_symm_sector.first][i][j].has_elements_)
          {
            int matrix_dimension=2*sign_symm_sector.second.block_shifts[i][j];
 
            matrices.push_back(Expr::mul(y_->index(op.second), As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension,matrix_dimension)));

          }
          
        }

        if(matrices.size()>0)
        {

        Expression::t ee=matrices[0];
        for(int n=1; n<matrices.size(); n++)
        {
          ee=Expr::add(ee, matrices[n]);
        }

        M_->constraint(ee, Domain::inPSDCone());
     
        }
      } 
      }

    }
    std::cout<< "Finished generating the PSD constraints"<<std::endl;
  return;
}
   Expression::t get_costfunction()
  {
   
    auto b_arr = monty::new_array_ptr<double>(b_);
    return Expr::dot(b_arr,y_);
  }

};


class momentum_symmetry_solver_sos: public momentum_basis_xy{
public:
   std::map<int,std::vector<std::vector<Expression::t>>> Xs_;
  // here we store the C matrices (the constants)
   std::map<int,std::vector<std::vector<Matrix::t>>> Cs_;
   std::map<int,std::vector<std::vector<Matrix::t>>> zeros_;
momentum_symmetry_solver_sos(int L, basis_structure operators,Model::t M, std::string permuts="xyz"):momentum_basis_xy(L, operators,M, permuts)
{
     
     for(auto sign_symm_sector :sectors_)
    {
        Xs_[sign_symm_sector.first]={};
        Cs_[sign_symm_sector.first]={};
        zeros_[sign_symm_sector.first]={};
        
  
        for(int i=0; i<L; i++)
        {
          Xs_[sign_symm_sector.first].push_back({});
          Cs_[sign_symm_sector.first].push_back({});
          zeros_[sign_symm_sector.first].push_back({});
          for(int j=0; j<L; j++)
          {
            int matrix_dimension=2*sign_symm_sector.second.block_shifts[i][j];
            auto X=M_->variable("X_"+std::to_string(sign_symm_sector.first)+"_"+std::to_string(i)+std::to_string(j), Domain::inPSDCone(matrix_dimension));
  
            // This is "minus" x, thus, we must replace all x with neg(x)
            Xs_[sign_symm_sector.first][i].push_back(Expr::neg(X));

          }
        }

        }
 
}
bool check_if_op(op_vec v0, std::string nam)
{

           auto all_t=generate_all_translations(v0, L_);
	  bool found=false;
	  
	     	 for(auto op_t: all_t)
	   {
          auto all_ty=generate_all_translations_y(op_t, L_,1);
	     for(auto op_ty: all_ty)
	      {
  if(nam==print_op(op_ty))
  {std::cout<< nam<<std::endl;
  found==true;
 
  }
     }
        }
return found;
}

void fix_constrains(){
 


    //std::map<std::string, Expression::t> expressions_;
    std::vector<Expression::t> expressions_(total_refs_.size(),Expr::constTerm( 0) );
    //Expression::t expressions_ =Expr::constTerm( total_refs_.size(),0.);

    for(auto sign_symm_sector :sectors_)
    {
      
      for(int i=0; i<L_; i++)
      {
       for(int j=0; j<L_; j++)
      {
      

        for(auto op :total_refs_ )
        {
          int matrix_dimension=2*sign_symm_sector.second.block_shifts[i][j];
            if(op.first=="0"){
              continue;
          //     if(As_[op.first][sign_symm_sector.first][i][j].has_elements_)
          //   {
              
          // //     auto zero=As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension,matrix_dimension);
          // //   zeros_[sign_symm_sector.first][i].push_back(zero);
          //   }
             }
           if(op.first=="1")
           {
             if(As_[op.first][sign_symm_sector.first][i][j].has_elements_)
           {
            auto C=As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension,matrix_dimension);
            Cs_[sign_symm_sector.first][i].push_back(C);
           // std::cout<< C->toString()<<std::endl;
            //ones+=1;
           }
            }
           else{
            if(As_[op.first][sign_symm_sector.first][i][j].has_elements_)
           {

            int el=total_refs_.at(op.first);
            //std::cout<< el <<" x "<<op.first <<" "<<As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension,matrix_dimension)->toString()<<std::endl;
            expressions_[el]=Expr::add(expressions_[el], Expr::dot(As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension,matrix_dimension),(Xs_[sign_symm_sector.first][i][j]) ));
            
           }
           }
          
         }


        }
      } 
      }

        for(auto sign_symm_sector :sectors_)
    {
     

      for(auto l : Cs_[sign_symm_sector.first])
      {
     
      }
   
    }
       
      for(auto a: total_refs_)
      {
        if(a.first!="1" && a.first!="0")
        {
        
          int el=total_refs_.at(a.first);
         
        M_->constraint( expressions_[el], Domain::equalsTo(-1*b_[el]));
         
        }
      }

    
    std::cout<< "Finished generating the PSD constraints ones "<<std::endl;
  return;
}
   Expression::t get_costfunction()
  {
    Expression::t ee=Expr::constTerm(0.);
       for(auto sign_symm_sector :sectors_)
    {
      
      for(int i=0; i<L_; i++)
      {
       for(int j=0; j<L_; j++)
      {
       
      ee=Expr::add(ee, Expr::dot(Cs_[sign_symm_sector.first][i][j],(Xs_[sign_symm_sector.first][i][j])));
     
    
      } } }
      return ee;
  }

};