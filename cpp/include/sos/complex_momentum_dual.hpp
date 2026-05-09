#pragma once
#include "fusion.h"
#include "spins.hpp"
#include "symmetries.hpp"
#include <unordered_map>
#include <memory>
#include "util.hpp"
#include "lattices.hpp"
#include <cassert>
#include "reduced_dms.hpp"
#include "operator_operations.hpp"
using namespace mosek::fusion;
using namespace monty;

using symmetry_sector = std::map<int, std::vector<std::vector<matrix_organizer>>>;

// // implementing momentum symmetrie in x and y direction
template <typename Lattice>
class momentum_block
{
public:
  std::vector<std::vector<Variable::t>> blocks_;
  int sign_sector_{0};
  std::vector<std::vector<int>> block_shifts;
  Lattice &lattice_;

  Eigen::MatrixXcd &FTx_;
  Eigen::MatrixXcd &FTy_;

  momentum_block(Lattice &lattice, Model::t M, int sign_sector, Eigen::MatrixXcd &FTy, Eigen::MatrixXcd &FTx, std::string sector_label = "") : lattice_(lattice), sign_sector_(sign_sector), FTy_(FTy), FTx_(FTx)
  {
    // std::cout << FTx_ << std::endl;
    // std::cout << FTy_ << std::endl;
  }
  void initialize_blocks_zero(std::map<std::string, symmetry_sector> &As)
  {

    int dim_0 = lattice_.states_[sign_sector_].size() + 1; // dimension of 0th block
    int dim_x = lattice_.states_[sign_sector_].size();     // dimension of other blocks
    std::cout<< "d im 1 "<<dim_0<<std::endl;
    std::cout<< "d im x "<<dim_x<<std::endl;
    block_shifts.push_back({});

    // initializing block shifts
    block_shifts[0].push_back(dim_0);
    for (int i = 1; i < lattice_.Lx_; i++)
    {

      block_shifts[0].push_back(dim_x);
    }

    for (int i = 1; i < lattice_.Ly_; i++)
    {
      block_shifts.push_back({});
      for (int j = 0; j < lattice_.Lx_; j++)
      {
        block_shifts[i].push_back(dim_x);
      }
    }

    As["1"][sign_sector_][0][0].add_values({0, 0}, 1. / 2);
    As["1"][sign_sector_][0][0].add_values({dim_0, dim_0}, 1. / 2);

    //   //     // The "c" terms first row and column in block 0
    int i = 0;

    for (auto it = lattice_.states_[sign_sector_].begin(); it != lattice_.states_[sign_sector_].end(); ++it)
    {
      auto op = *it;

      // get normal form
      auto [coeff, nf] = get_normal_form(op);
      // get translation invariant representation

      auto ti_key = lattice_.TI_map_.at(print_op(nf)).first;

      //auto el = this->lattice_.variable_map_.at(ti_key);

      if (std::abs(coeff.real()) > 1e-9)
      {
        //std::cout<< "xx "<<ti_key<<std::endl;
         As[ti_key][sign_sector_][0][0].add_values({0, i + 1}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
        As[ti_key][sign_sector_][0][0].add_values({i + 1, 0}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
        As[ti_key][sign_sector_][0][0].add_values({dim_0, i + 1 + dim_0}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
        As[ti_key][sign_sector_][0][0].add_values({i + 1 + dim_0, dim_0}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
      }
      assert(std::abs(coeff.imag()) < 1e-9);

      i++;
    }

    return;
  }
  void initialize_blocks_general()
  {

    int dim_x = lattice_.states_[sign_sector_].size(); // operators_.size(); // dimension of other blocks

    for (int j = 0; j < lattice_.Ly_; j++)
    {
      block_shifts.push_back({});

      for (int i = 0; i < lattice_.Lx_; i++)
      {
        block_shifts[j].push_back(dim_x);
      }
    }

    return;
  }
  void initialize_blocks(std::map<std::string, symmetry_sector> &As)
  {
    if (sign_sector_ == 0)
    {
      initialize_blocks_zero(As);
    }
    else
    {
      initialize_blocks_general();
    }
    return;
  }

  void generate_block(std::map<std::string, symmetry_sector> &As)
  {
    //     const auto start{std::chrono::steady_clock::now()};

    int i = 0;

    for (auto it1 = lattice_.states_[sign_sector_].begin(); it1 != lattice_.states_[sign_sector_].end(); ++it1)
    {
      int j = i;
      for (auto it2 = it1; it2 != lattice_.states_[sign_sector_].end(); ++it2)
      {

        for (int mat_pos_y = 0; mat_pos_y < lattice_.Ly_; mat_pos_y++)
        {
          for (int mat_pos_x = 0; mat_pos_x < lattice_.Lx_; mat_pos_x++)
          {

            // 			      // determines if first block of zeroth moment blocks
            int shift = block_shifts[mat_pos_y][mat_pos_x] % lattice_.states_[sign_sector_].size();

            // 			      // gives the shift between real and complex components
            int dim = block_shifts[mat_pos_y][mat_pos_x];

            for (int pos_y = 0; pos_y < lattice_.Ly_; pos_y++)
            {
              std::complex<double> FT_factor_y = FTy_(pos_y, mat_pos_y);

              for (int pos_x = 0; pos_x < lattice_.Lx_; pos_x++)
              {
                std::complex<double> FT_factor_x = FTx_(pos_x, mat_pos_x);

                //              // to do, correct so that all terms appearing here appear in map

                auto construct = lattice_.generate_G_element_sos(*it1, *it2, pos_y, pos_x);
                // if (i == j)
                // {

                //   std::cout << construct.prefac_ << "  " << construct.op_ << std::endl;
                // }
                std::complex<double>
                    total_prefactor = construct.prefac_ * FT_factor_x * FT_factor_y;
                 
                // assert(std::abs(total_prefactor)<1e-9); maybe not include values  that are zero

                if (std::abs(total_prefactor.real()) > 1e-9)
                {

                  As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values({i + shift, j + shift}, 1. / 2 * total_prefactor.real());
                  As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values({i + shift + dim, j + shift + dim}, 1. / 2 * total_prefactor.real());
                  if (i != j)
                  {
                    As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values({j + shift, i + shift}, 1. / 2 * total_prefactor.real());
                    As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values({j + shift + dim, i + shift + dim}, 1. / 2 * total_prefactor.real());
                  }
                }
                if (std::abs(total_prefactor.imag()) > 1e-9)
                {

                  // assert(i != j);
                  //  X^T[0,1]-X[0,1]=-H[0,1]

                  As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values({i + shift, j + shift + dim}, -1. / 2 * total_prefactor.imag());
                  As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values({j + shift, i + shift + dim}, 1. / 2 * total_prefactor.imag());
                  if (i != j)
                  {
                    As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values({i + shift + dim, j + shift}, 1. / 2 * total_prefactor.imag());
                    As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values({j + shift + dim, i + shift}, -1. / 2 * total_prefactor.imag());
                  }
                }
              }
            }
          }
        }

        j += 1;
      }
      i += 1;
    }

    return;
  }
};
template <typename Lattice>
class momentum_basis
{
  // note, the first sector must contain the unit element
  // solves min(by), with sum_i y_i A_i <<C
public:
  Model::t M_;
  std::map<int, momentum_block<Lattice>> sectors_;
  std::string sector_;

  Eigen::MatrixXcd FTx_;
  Eigen::MatrixXcd FTy_;
  Parameter::t b_;
  Parameter::t energy_vec_; // if b is observable then the energy is stored here
  bool bounding_observable_{false};
  std::map<std::string, Parameter::t> energy_bounds_; // contains two elements, upper bound lower bound
  Lattice &lattice_;
  // enforcing constarans ye nergy_vec_<=E_upper
  // contains the matrices As, for each sign symmetrye we have LxL blocks
  std::map<std::string, symmetry_sector> As_;
  // for the reduced density matrix
  std::map<rdm_operator, std::map<std::string, Matrix::t>> sigmas_;
  std::vector<Parameter::t>  linear_constraints_coefficients_;
   
  momentum_basis(Lattice &lattice, Model::t M, rdms_struct rdms) : lattice_(lattice), M_(M)
  {
    std::cout << "start" << std::endl;
    FTx_ = Eigen::MatrixXcd(lattice_.Lx_, lattice_.Lx_);
    for (int i = 0; i < lattice_.Lx_; i++)
    {
      for (int j = 0; j < lattice_.Lx_; j++)
      {
        std::complex<double> phase(0., -2. * i * j * pi / lattice_.Lx_);

        FTx_(i, j) = std::exp(phase);///std::sqrt(lattice_.Ly_);
      }
    }
    FTy_ = Eigen::MatrixXcd(lattice_.Ly_, lattice_.Ly_);
    for (int i = 0; i < lattice_.Ly_; i++)
    {
      for (int j = 0; j < lattice_.Ly_; j++)
      {
        std::complex<double> phase(0., -2. * i * j * pi / lattice_.Ly_);

        FTy_(i, j) = std::exp(phase);///std::sqrt(lattice_.Ly_);
      }
    }
    std::cout << "start initializeing blocks" << std::endl;
    for (auto it = lattice_.states_.begin(); it != lattice_.states_.end(); ++it)
    {

      auto Block = momentum_block(lattice_, M_, it->first, FTy_, FTx_, std::to_string(it->first));
      sectors_.insert({it->first, Block});
    }
    std::cout << "start initializeing maps" << std::endl;
    initialize_all_maps(rdms);

    std::cout << "size TI map " << lattice_.TI_map_.size() << std::endl;
    // for (auto a : lattice_.TI_map_)
    // {
    //   std::cout << a.first << "-> " << a.second.first << "   " << a.second.second << std::endl;
    // }
    std::cout << "size total refs " << lattice.variable_map_.size() << std::endl;

    for (auto it = lattice.variable_map_.begin(); it != lattice.variable_map_.end(); it++)
    {
      As_.insert({it->first, symmetry_sector()});
      for (auto it_sign_sector = sectors_.begin(); it_sign_sector != sectors_.end(); ++it_sign_sector)
      {
        As_[it->first][it_sign_sector->first] = {};

        for (int i = 0; i < lattice_.Ly_; i++)
        {
          As_[it->first][it_sign_sector->first].push_back({});
          for (int j = 0; j < lattice_.Lx_; j++)
          {
            As_[it->first][it_sign_sector->first][i].push_back(matrix_organizer());
          }
        }
      }
    }
     
    std::cout << "initialze blocks " << std::endl;
    for (auto &sector : sectors_)
    {
    
      sector.second.initialize_blocks(As_);
   
    }
    std::cout << "finished initializeing blocks" << std::endl;
    for (auto it_2 = sectors_.begin(); it_2 != sectors_.end(); ++it_2)
    {
      it_2->second.generate_block(As_);
  
    }
    std::cout << "finished making the As matrices" << std::endl;

    return;
  };
  void initialize_all_maps(rdms_struct rdms)
  {
    this->lattice_.generate_TI_map();
    if (rdms.size() > 0)
    {

      generate_rdms(rdms);
    }
    this->lattice_.make_map();

    b_ = M_->parameter("b", lattice_.variable_map_.size());
  }
  void set_b(std::vector<double> b)
  {

    auto a = monty::new_array_ptr<double>(b);
    b_->setValue(a);
    return;
  }
  void set_energy_vec(std::vector<double> energy_vec, double E_upper, double E_lower)
  {
    if (energy_bounds_.size() < 2)
    {
      bounding_observable_ = true;
      energy_vec_ = M_->parameter("energy vec", lattice_.variable_map_.size());
      energy_bounds_["E_upper"] = M_->parameter("E_upper");
      energy_bounds_["E_lower"] = M_->parameter("E_lower");
    }

    auto a = monty::new_array_ptr<double>(energy_vec);
    energy_vec_->setValue(a);

    std::cout << "bounding " << bounding_observable_ << std::endl;
    energy_bounds_["E_upper"]->setValue(E_upper);
    energy_bounds_["E_lower"]->setValue(E_lower);
    std::cout << "bounds up " << *(energy_bounds_["E_upper"]->getValue()) << std::endl;
    std::cout << "bounds low " << *(energy_bounds_["E_lower"]->getValue()) << std::endl;
    std::cout << "diff " << (*(energy_bounds_["E_upper"]->getValue()))[0] - (*(energy_bounds_["E_lower"]->getValue()))[0] << std::endl;

    return;
  }
  void set_linear_constraints_vec(std::vector<std::vector<double>>  linear_constraints)
  {
    if(linear_constraints_coefficients_.size()<1)
    {
      for(int i=0; i<linear_constraints.size(); i++)
      {
        linear_constraints_coefficients_.push_back(M_->parameter("linear_constraint_"+std::to_string(i), lattice_.variable_map_.size()));
      }
    }
    for(int i=0; i< linear_constraints.size(); i++)
    {
     
      auto a = monty::new_array_ptr<double>(linear_constraints[i]);
      linear_constraints_coefficients_[i]->setValue(a);
    }
    return;
  }
  void generate_rdms(rdms_struct rdms)
  {

    auto offset = lattice_.states_[1][0][0].offset_; // change this to be derived from baso

    int i = 0;
    std::cout << "rdms size " << rdms.rdms.size() << std::endl;
    for (auto site : rdms.rdms)
    {

      i++;
      auto sigmas_temp = lattice_.generate_rdms_primal_cp(site, offset);
      sigmas_.insert({site, sigmas_temp});
    }
    return;
  }
};
template <typename Lattice>
class momentum_symmetry_solver_dual : public momentum_basis<Lattice>
{
public:
  Variable::t y_;
  momentum_symmetry_solver_dual(Lattice &lattice, Model::t M, rdms_struct rdms) : momentum_basis<Lattice>(lattice, M, rdms)
  {
    std::cout << "start " << std::endl;
    y_ = this->M_->variable("T", this->lattice_.variable_map_.size());
    // fix 1
    auto el = this->lattice_.variable_map_.at("1");
    this->M_->constraint(y_->index(el), Domain::equalsTo(1.0));
    std::cout << "start " << std::endl;
    // fix zero
    el = this->lattice_.variable_map_.at("0");
    this->M_->constraint(y_->index(el), Domain::equalsTo(0.0));
    std::cout << "start " << std::endl;
  }
  void fix_constrains()
  {

    // iterate over sign sector
    for (auto sign_symm_sector : this->sectors_)
    {

      for (int i = 0; i < this->lattice_.Ly_; i++)
      {
        for (int j = 0; j < this->lattice_.Lx_; j++)
        {
          std::vector<Expression::t> matrices;

          for (auto op : this->lattice_.variable_map_)
          {
            if (op.first == "0")
            {
              continue;
            }
            if (this->As_[op.first][sign_symm_sector.first][i][j].has_elements_)
            {
              int matrix_dimension = 2 * sign_symm_sector.second.block_shifts[i][j];

              matrices.push_back(Expr::mul(y_->index(op.second), this->As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension, matrix_dimension)));
            }
          }

          if (matrices.size() > 0)
          {

            Expression::t ee = matrices[0];
            for (int n = 1; n < matrices.size(); n++)
            {
              ee = Expr::add(ee, matrices[n]);
            }

            this->M_->constraint(ee, Domain::inPSDCone());
          }
        }
      }
    }
    std::cout << "Finished generating the PSD constraints" << std::endl;
    for (auto state : this->sigmas_)
    {
      Expression::t ee = Expr::constTerm(state.second["1"]);
      // matrices[0]
      for (auto op_string : state.second)
      {
        if (op_string.first != "1")
        {
          ee = Expr::add(ee, Expr::mul(y_->index(this->lattice_.variable_map_[op_string.first]), op_string.second));
        }
      }
      this->M_->constraint(ee, Domain::inPSDCone());
    }
    std::cout << "Finished density matrices " << std::endl;
    // bounding energy
    if (this->bounding_observable_)
    {
      // std::cout << "introduing bounds" << std::endl;
      // std::cout << " upper " << this->energy_bounds_["E_upper"]->index(0) << std::endl;
      // std::cout << " lower " << this->energy_bounds_["E_lower"]->index(0) << std::endl;
      // this->M_->constraint(Expr::dot(this->energy_vec_, y_), Domain::lessThan(this->energy_bounds_["E_upper"]->index(0)));
      // this->M_->constraint(Expr::dot(this->energy_vec_, y_), Domain::greaterThan(this->energy_bounds_["E_lower"]->index(0)));
    }
    return;
  }
  Expression::t get_costfunction()
  {

    return Expr::dot(this->b_, y_);
  }
};
// struct linear_constraint{
// std::string key;
// double value;
// linear_constraint(std::string key, double value):key_(key), value_(value){};

// };
template <typename Lattice>
class momentum_symmetry_solver_sos : public momentum_basis<Lattice>
{
public:
  std::map<int, std::vector<std::vector<Expression::t>>> Xs_;
  // Lambdas are the Lagrangian stemmeing from psd density matrices
  std::map<rdm_operator, Expression::t> Lambdas_;

  // here we store the C matrices (the constants)
  // std::map<int, std::vector<std::vector<Matrix::t>>> Cs_;
  std::map<int, std::vector<std::vector<Matrix::t>>> zeros_;
  // variables introduced to bound the energy
  std::vector<Variable::t> energy_bouding_variables_;
  bool maximize_{true}; // if cost function is a maximization problem
  // to enforce 0=0 and 1=1
  Variable::t eta;
  Variable::t epsilon;

  // dummy variable used to constrain magnetization
  // todo add parameter that you can set and make a certain constraint fulfilled

  Variable::t delta;
  // a vector where each element is a constraint 
  std::vector<Variable::t> linear_constraints_variable_; 
  momentum_symmetry_solver_sos(Lattice &lattice, Model::t M, rdms_struct rdms, bool maximize = true) : maximize_(maximize), momentum_basis<Lattice>(lattice, M, rdms)
  {
    if (maximize_)
    {
      eta = this->M_->variable("eta", Domain::greaterThan(0.));
      epsilon = this->M_->variable("epsilon", Domain::greaterThan(0.));
    }
    else
    {
      eta = this->M_->variable("eta", Domain::lessThan(0.));
      epsilon = this->M_->variable("epsilon", Domain::lessThan(0.));
    }
    for (auto sign_symm_sector : this->sectors_)
    {
      Xs_[sign_symm_sector.first] = {};
      // Cs_[sign_symm_sector.first] = {};
      zeros_[sign_symm_sector.first] = {};

      for (int i = 0; i < this->lattice_.Ly_; i++)
      {
        Xs_[sign_symm_sector.first].push_back({});
        // Cs_[sign_symm_sector.first].push_back({});
        zeros_[sign_symm_sector.first].push_back({});

        for (int j = 0; j < this->lattice_.Lx_; j++)
        {

          int matrix_dimension = 2 * sign_symm_sector.second.block_shifts[i][j];

          auto X = this->M_->variable("X_" + std::to_string(sign_symm_sector.first) + "_" + std::to_string(i) + std::to_string(j), Domain::inPSDCone(matrix_dimension));

          // This is "minus" x, thus, we must replace all x with neg(x)
          if (maximize_)
          {
            Xs_[sign_symm_sector.first][i].push_back(Expr::neg(X));
          }
          else
          {
            Xs_[sign_symm_sector.first][i].push_back((X));
          }
        }
      }
    }

    int i = 0;
    for (auto op : rdms.rdms)
    {
      auto dm_dim = std::pow(2, op.size());

      auto beta = this->M_->variable("betas_" + std::to_string(i), Domain::inPSDCone(2 * dm_dim));
      i++;
      if (maximize_)
      {
        Lambdas_[op] = Expr::neg(beta);
      }
      else
      {
        Lambdas_[op] = (beta);
      }
    }
  }

  void fix_constrains()
  {
    // wanting to solve the sdp
    // Tr<X,C>, s.t. for all i, Tr<X,A_i>=b_i
    if (this->bounding_observable_)
    {
      std::cout << "true bounding observable " << std::endl;
      if (maximize_)
      {
        energy_bouding_variables_.push_back(this->M_->variable("upper energy", Domain::greaterThan(0.)));
        energy_bouding_variables_.push_back(this->M_->variable("lower energy", Domain::lessThan(0.)));
      }
      else
      {
        energy_bouding_variables_.push_back(this->M_->variable("upper energy", Domain::lessThan(0.)));
        energy_bouding_variables_.push_back(this->M_->variable("lower energy", Domain::greaterThan(0.)));
      }
    }
    
      for(int i=0; i<this->linear_constraints_coefficients_.size(); i++)
      {
        linear_constraints_variable_.push_back(this->M_->variable("LC"+std::to_string(i)));
      }
    

    std::vector<Expression::t> expressions_(this->lattice_.variable_map_.size(), Expr::constTerm(0));

    for (auto sign_symm_sector : this->sectors_)
    {

      for (int i = 0; i < this->lattice_.Ly_; i++)
      {
        for (int j = 0; j < this->lattice_.Lx_; j++)
        {

          for (auto op : this->lattice_.variable_map_)
          {
            int matrix_dimension = 2 * sign_symm_sector.second.block_shifts[i][j];
            if (op.first == "0")
            {
              continue;
            }
            else
            {
              if (op.first == "1")
              {
                if (this->As_[op.first][sign_symm_sector.first][i][j].has_elements_)
                {
                  int el = this->lattice_.variable_map_.at(op.first);
                  // generatin the C matrix blocks (the one that will be used om the cost function min(C,X))
                  // auto C = this->As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension, matrix_dimension);
                  // Cs_[sign_symm_sector.first][i].push_back(C);
                  expressions_[el] = Expr::add(expressions_[el], Expr::dot(this->As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension, matrix_dimension), (Xs_[sign_symm_sector.first][i][j])));
                }
              }
              else
              {
                if (this->As_[op.first][sign_symm_sector.first][i][j].has_elements_)
                {

                  int el = this->lattice_.variable_map_.at(op.first);
                  // making the constrains Tr<A_,X> which we will assign to b_i later
                  expressions_[el] = Expr::add(expressions_[el], Expr::dot(this->As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension, matrix_dimension), (Xs_[sign_symm_sector.first][i][j])));
                }
              }
            }
          }
        }
      }
    }
    // Enforcing reduced density matrices. Todo, see if this also can be simplified by removing moving "1" into this loop
    std::cout << "start generating constarins for rdms " << std::endl;
    for (auto lambda_ : Lambdas_)
    {
      for (auto string_and_matrix : this->sigmas_[lambda_.first])
      {
        if (string_and_matrix.first != "1")
        {

          int el = this->lattice_.variable_map_.at(string_and_matrix.first);
          auto a = lambda_.second;

          expressions_[el] = Expr::add(expressions_[el], (Expr::dot(lambda_.second, string_and_matrix.second)));
        }
      }
    }
    // adding the constrains enforcing energy </> to lower/upper bound
    if (this->bounding_observable_)
    {
      auto exp_temporary = Expr::mul(Expr::add(energy_bouding_variables_[0], energy_bouding_variables_[1]), this->energy_vec_);
      for (int i = 0; i < this->energy_vec_->getSize(); i++)
      {
        // energy_vec_->index(i)

        expressions_[i] = Expr::add(expressions_[i], (exp_temporary->index(i)));
      }
    }
    // adding lunear constarins
    if(linear_constraints_variable_.size()>0)
    {
      auto exp_temporary = Expr::mul(linear_constraints_variable_[0], this->linear_constraints_coefficients_[0]);
     
    
    for(int j=1; j<linear_constraints_variable_.size(); j++)
    {
      exp_temporary =Expr::add(exp_temporary,Expr::mul(linear_constraints_variable_[j], this->linear_constraints_coefficients_[j]));
      for (int i = 0; i < linear_constraints_variable_[0]->getSize(); i++)
      {
        // energy_vec_->index(i)

        expressions_[i] = Expr::add(expressions_[i], (exp_temporary->index(i)));
      }
      
    }
  }



    // adding a constant term for the 1:
    int el = this->lattice_.variable_map_.at("1");
    expressions_[el] = Expr::add(expressions_[el], epsilon);

    for (auto a : this->lattice_.variable_map_)
    {
      if (a.first != "0")
      {

        int el = this->lattice_.variable_map_.at(a.first);
        //=-1*b[el]
        this->M_->constraint(Expr::add(expressions_[el], this->b_->index(el)), Domain::equalsTo(0.));
      }
      // fixing that zero is zero
      // normally \eta=0 but we just eliminate eta
      if (a.first == "0")
      {
        int el = this->lattice_.variable_map_.at(a.first);
        this->M_->constraint(expressions_[el], Domain::equalsTo(0.));
      }
    }

    std::cout << "Finished generating the PSD constraints ones " << std::endl;
    return;
  }
  Expression::t get_costfunction()
  {
    Expression::t ee = Expr::constTerm(0.);
    ee = Expr::add(ee, Expr::neg(epsilon));

    // Adding matrices for the positive definite constrain of the RDMs
    for (auto lambda_ : Lambdas_)
    {

      ee = Expr::add(ee, Expr::dot(lambda_.second, this->sigmas_[lambda_.first]["1"]));
    }
    if (this->bounding_observable_)
    {

      ee = Expr::add(ee, Expr::neg(Expr::mul(this->energy_bounds_["E_upper"], this->energy_bouding_variables_[0])));
      ee = Expr::add(ee, Expr::neg(Expr::mul(this->energy_bounds_["E_lower"], this->energy_bouding_variables_[1])));
    }
    return ee;
  }
};