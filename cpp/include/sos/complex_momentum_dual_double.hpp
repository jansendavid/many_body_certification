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

template <typename Lattice>
class momentum_block_double
{
public:
  int sign_sector_{0};
  std::vector<std::vector<int>> block_shifts;
  Lattice &lattice_;

  Eigen::MatrixXcd &FTx_;
  Eigen::MatrixXcd &FTy_;

  momentum_block_double(Lattice &lattice, int sign_sector, Eigen::MatrixXcd &FTy, Eigen::MatrixXcd &FTx)
      : lattice_(lattice), sign_sector_(sign_sector), FTy_(FTy), FTx_(FTx)
  {
  }
  void initialize_blocks_zero(std::map<std::string, symmetry_sector> &As)
  {


    int dim_0 = lattice_.states_[sign_sector_][0].size()+lattice_.states_[sign_sector_][1].size() + 1; // dimension of 0th block
    int dim_x = lattice_.states_[sign_sector_][0].size()+lattice_.states_[sign_sector_][1].size();     // dimension of other blocks

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

    for (auto it = lattice_.states_[sign_sector_][0].begin(); it != lattice_.states_[sign_sector_][0].end(); ++it)
    {
      auto op = *it;
      // get normal form
      auto [coeff, nf] = get_normal_form(op);
      // get translation invariant representation

      auto ti_key = op_key_label(lattice_.TI_map_.at(key_dir_pos(nf)).first);

      //auto el = this->lattice_.variable_map_.at(ti_key);

      if (std::abs(coeff.real()) > 1e-9)
      {

        As[ti_key][sign_sector_][0][0].add_values({0, i + 1}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
        As[ti_key][sign_sector_][0][0].add_values({i + 1, 0}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
        As[ti_key][sign_sector_][0][0].add_values({dim_0, i + 1 + dim_0}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
        As[ti_key][sign_sector_][0][0].add_values({i + 1 + dim_0, dim_0}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
      }
      assert(std::abs(coeff.imag()) < 1e-9);

      i++;
    }
    i = 0;
    int shift=lattice_.states_[sign_sector_][0].size();
    for (auto it = lattice_.states_[sign_sector_][1].begin(); it != lattice_.states_[sign_sector_][1].end(); ++it)
    {
      auto op = *it;
      // get normal form
      auto [coeff, nf] = get_normal_form(op);
      // get translation invariant representation

      auto ti_key = op_key_label(lattice_.TI_map_.at(key_dir_pos(nf)).first);
      //auto el = this->lattice_.variable_map_.at(ti_key);

      if (std::abs(coeff.real()) > 1e-9)
      {

        As[ti_key][sign_sector_][0][0].add_values({0, i + 1+shift}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
        As[ti_key][sign_sector_][0][0].add_values({i + 1+shift, 0}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
        As[ti_key][sign_sector_][0][0].add_values({dim_0, i + 1 + dim_0+shift}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
        As[ti_key][sign_sector_][0][0].add_values({i + 1 + dim_0+shift, dim_0}, 1. / 2 * coeff.real() * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_));
      }
      assert(std::abs(coeff.imag()) < 1e-9);

      i++;
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
  void initialize_blocks_general()
  {


    int dim_x = lattice_.states_[sign_sector_][0].size()+lattice_.states_[sign_sector_][1].size(); // operators_.size(); // dimension of other blocks

    for (int j = 0; j < lattice_.Ly_; j++)
    {
      block_shifts.push_back({});

      for (int i = 0; i < lattice_.Lx_; i++)
      {
        block_shifts[j].push_back(dim_x);
      }
    }
  //}
    return;
  }

 void run_loop(std::vector<op_vec>& operator_1,std::vector<op_vec>& operator_2,std::map<std::string, symmetry_sector> &As, std::complex<double> fac_orig, std::pair<int,int> shift)
 {
  const int Ly = lattice_.Ly_;
  const int Lx = lattice_.Lx_;
  const int n_states = static_cast<int>(lattice_.states_[sign_sector_][0].size() +
                                        lattice_.states_[sign_sector_][1].size());

  int i = 0;
  for (auto it1 = operator_1.begin(); it1 != operator_1.end(); ++it1)
  {
    int j = 0;
    for (auto it2 = operator_2.begin(); it2 != operator_2.end(); ++it2)
    {
      // generate_G_element_sos depends only on (it1, it2, pos_y, pos_x), not on mat_pos.
      for (int pos_y = 0; pos_y < Ly; ++pos_y)
      {
        for (int pos_x = 0; pos_x < Lx; ++pos_x)
        {
          auto construct = lattice_.generate_G_element_sos(*it1, *it2, pos_y, pos_x);
          if (construct.op_ == "0")
            continue;

          assert(std::abs((construct.prefac_ * fac_orig).imag()) < 1e-9);

          for (int mat_pos_y = 0; mat_pos_y < Ly; ++mat_pos_y)
          {
            const std::complex<double> FT_factor_y = FTy_(pos_y, mat_pos_y);

            for (int mat_pos_x = 0; mat_pos_x < Lx; ++mat_pos_x)
            {
              const int shift_initial = block_shifts[mat_pos_y][mat_pos_x] % n_states;
              const int dim = block_shifts[mat_pos_y][mat_pos_x];
              const std::complex<double> FT_factor_x = FTx_(pos_x, mat_pos_x);
              const std::complex<double> total_prefactor =
                  construct.prefac_ * fac_orig * FT_factor_x * FT_factor_y;

              if (std::abs(total_prefactor.real()) > 1e-9)
              {
                As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values(
                    {i + shift.first + shift_initial, j + shift.second + shift_initial},
                    1. / 2 * total_prefactor.real());
                As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values(
                    {i + shift.first + dim + shift_initial, j + shift.second + dim + shift_initial},
                    1. / 2 * total_prefactor.real());
              }
              if (std::abs(total_prefactor.imag()) > 1e-9)
              {
                As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values(
                    {i + shift.first + shift_initial, j + shift.second + dim + shift_initial},
                    -1. / 2 * total_prefactor.imag());
                As[construct.op_][sign_sector_][mat_pos_y][mat_pos_x].add_values(
                    {i + shift.first + dim + shift_initial, j + shift.second + shift_initial},
                    1. / 2 * total_prefactor.imag());
              }
            }
          }
        }
      }

      j += 1;
    }
    i += 1;
  }
 }
  void generate_block(std::map<std::string, symmetry_sector> &As)
  {
std::complex<double> prefac(1.,0.);
int dim=lattice_.states_[sign_sector_][0].size();
std::pair<int,int> shift={0,0};
run_loop(lattice_.states_[sign_sector_][0],lattice_.states_[sign_sector_][0],As, prefac, shift);

prefac={1.,0.};

shift={dim,dim};
 run_loop(lattice_.states_[sign_sector_][1],lattice_.states_[sign_sector_][1],As, prefac, shift);

 prefac={0.,1.};

shift={0,dim};
 run_loop(lattice_.states_[sign_sector_][0],lattice_.states_[sign_sector_][1],As, prefac, shift);

 prefac={0.,-1.};

shift={dim,0};
 run_loop(lattice_.states_[sign_sector_][1],lattice_.states_[sign_sector_][0],As, prefac, shift);


    return;
  }
};
template <typename Lattice>
class momentum_basis_double
{
  // note, the first sector must contain the unit element
  // solves min(by), with sum_i y_i A_i <<C
public:
  Model::t M_;
  std::map<int, momentum_block_double<Lattice>> sectors_;

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
  std::map<rdm_operator, std::vector<std::map<std::string, Matrix::t>>> sigmas_;
  Matrix::t Psp;
  int nr_of_linear_constraints{0};
  bool U1=false;

  momentum_basis_double(Lattice &lattice, Model::t M, rdms_struct rdms, bool U1) : lattice_(lattice), M_(M), U1(U1)
  {
    FTx_ = Eigen::MatrixXcd(lattice_.Lx_, lattice_.Lx_);
    for (int i = 0; i < lattice_.Lx_; i++)
    {
      for (int j = 0; j < lattice_.Lx_; j++)
      {
        std::complex<double> phase(0., -2. * i * j * pi / lattice_.Lx_);

        FTx_(i, j) = std::exp(phase);
      }
    }
    FTy_ = Eigen::MatrixXcd(lattice_.Ly_, lattice_.Ly_);
    for (int i = 0; i < lattice_.Ly_; i++)
    {
      for (int j = 0; j < lattice_.Ly_; j++)
      {
        std::complex<double> phase(0., -2. * i * j * pi / lattice_.Ly_);

        FTy_(i, j) = std::exp(phase);
      }
    }
    for (auto it = lattice_.states_.begin(); it != lattice_.states_.end(); ++it)
    {
      auto Block = momentum_block_double(lattice_, it->first, FTy_, FTx_);
      sectors_.insert({it->first, Block});
    }
    initialize_all_maps(rdms);

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

    for (auto &sector : sectors_)
      sector.second.initialize_blocks(As_);

    for (auto it_2 = sectors_.begin(); it_2 != sectors_.end(); ++it_2)
      it_2->second.generate_block(As_);

    return;
  };
  void initialize_all_maps(rdms_struct rdms)
  {
    this->lattice_.generate_TI_map_double();
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
    energy_bounds_["E_upper"]->setValue(E_upper);
    energy_bounds_["E_lower"]->setValue(E_lower);

    return;
  }
  void set_linear_constraints_vec(std::vector<std::vector<double>> linear_constraints)
  {
    if (nr_of_linear_constraints < 1)
      nr_of_linear_constraints = static_cast<int>(linear_constraints.size());

    const int m = static_cast<int>(lattice_.variable_map_.size());
    const int n = nr_of_linear_constraints;

    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> vals;
    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < m; ++j)
      {
        double v = linear_constraints[i][j];
        if (std::abs(v) > 1e-15)
        {
          rows.push_back(j);
          cols.push_back(i);
          vals.push_back(v);
        }
      }
    }
    Psp = Matrix::t(Matrix::sparse(
        m, n, monty::new_array_ptr<int>(rows), monty::new_array_ptr<int>(cols),
        monty::new_array_ptr<double>(vals)));
    return;
  }
  void generate_rdms(rdms_struct rdms)
  {

    auto offset = lattice_.states_[1][0][0][0].offset_;
    std::cout << "rdms size " << rdms.rdms.size() << std::endl;
    for (auto site : rdms.rdms)
    {
      if(U1)
      {
        auto sigmas_temp = lattice_.generate_rdms_primal_U1(site, offset); 
        sigmas_.insert({site, sigmas_temp});
      }
      else{
        auto sigmas_temp = lattice_.generate_rdms_primal_cp(site, offset); 
      
        sigmas_.insert({site, sigmas_temp});
      }
    
    }
    return;
  }
};
template <typename Lattice>
class momentum_symmetry_solver_dual_double : public momentum_basis_double<Lattice>
{
public:
  Variable::t y_;
  momentum_symmetry_solver_dual_double(Lattice &lattice, Model::t M, rdms_struct rdms, bool U1=false) : momentum_basis_double<Lattice>(lattice, M, rdms, U1)
  {
    y_ = this->M_->variable("T", this->lattice_.variable_map_.size());
    this->M_->constraint(y_, Domain::lessThan(1.0));
    this->M_->constraint(y_, Domain::greaterThan(-1.0));

    auto el = this->lattice_.variable_map_.at("1");
    this->M_->constraint(y_->index(el), Domain::equalsTo(1.0));

    auto it = this->lattice_.variable_map_.find("0");
    if (it != this->lattice_.variable_map_.end())
      this->M_->constraint(y_->index(it->second), Domain::equalsTo(0.0));
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
    for(auto& psd_mat:this->sigmas_ )
    {
      for(auto& elements: psd_mat.second)
      {
    for (auto& state : elements)
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
  }
}
    std::cout << "Finished density matrices " << std::endl;

    if (this->nr_of_linear_constraints > 0)
{
    auto vals = this->Psp->getValue();
    auto shape = this->Psp->getShape();
    const int m = (*shape)[0];
    const int n = (*shape)[1];
    std::cout << "adding linear constrains " << n << std::endl;

    for (int i = 0; i < n; ++i)
    {
        // extract column i as a std::vector
        std::vector<double> col(m);
        for (int j = 0; j < m; ++j)
            col[j] = (*vals)[j * n + i];

        this->M_->constraint(
            Expr::dot(monty::new_array_ptr<double>(col), y_),
            Domain::equalsTo(0.0));
    }
    }

    return;
  }
  Expression::t get_costfunction()
  {

    return Expr::dot(this->b_, y_);
  }
};
template <typename Lattice>
class momentum_symmetry_solver_sos_double : public momentum_basis_double<Lattice>
{
public:
  std::map<int, std::vector<std::vector<Expression::t>>> Xs_;
  std::map<rdm_operator, std::vector<Expression::t>> Lambdas_;

  std::vector<Variable::t> energy_bouding_variables_;
  bool maximize_{true};
  Variable::t epsilon;
  Variable::t linear_constraints_variable2_;
  std::vector<std::vector<Variable::t>> linear_constraints_for_block_equality_variable_;
  Constraint::t final_constraint_;
  Expression::t A_vector = nullptr;
  Expression::t Lamba_vector = nullptr;
  Expression::t LC_vector = nullptr;
  Expression::t epsilon_vec_flat = nullptr;

  static int nrblocks_y_for(const Lattice &lattice)
  {
    return (lattice.Ly_ % 2 == 0) ? (2 + lattice.Ly_ / 2 - 1) : (1 + lattice.Ly_ / 2);
  }

  static int nrblocks_x_for(const Lattice &lattice)
  {
    return (lattice.Lx_ % 2 == 0) ? (2 + lattice.Lx_ / 2 - 1) : (1 + lattice.Lx_ / 2);
  }

  static double matrix_trace(const Matrix::t &mat, int block_size)
  {
    auto data = mat->getDataAsArray();
    double tr = 0.0;
    for (int r = 0; r < block_size; ++r)
      tr += (*data)[r * block_size + r];
    return tr;
  }

  momentum_symmetry_solver_sos_double(Lattice &lattice, Model::t M, rdms_struct rdms, bool maximize = true, bool U1=false)
      : maximize_(maximize), momentum_basis_double<Lattice>(lattice, M, rdms, U1)
  {
    epsilon = this->M_->variable("epsilon");
    for(int i=0; i<int(this->lattice_.Ly_/2); i++)
    {
      linear_constraints_for_block_equality_variable_.push_back({});
      for(int j=0; j<int(this->lattice_.Lx_/2); j++)
      {
        linear_constraints_for_block_equality_variable_[i].push_back(this->M_->variable("block_equality_variable_"+std::to_string(i)+"_"+std::to_string(j)));
      }
    }
    const int nrblocks_y = nrblocks_y_for(this->lattice_);
    const int nrblocks_x = nrblocks_x_for(this->lattice_);
    for (auto sign_symm_sector : this->sectors_)
    {
      Xs_[sign_symm_sector.first] = {};
      for (int i = 0; i < nrblocks_y; i++)
      {
        Xs_[sign_symm_sector.first].push_back({});
        for (int j = 0; j < nrblocks_x; j++)
        {
          int matrix_dimension = 2 * sign_symm_sector.second.block_shifts[i][j];
          auto X = this->M_->variable("X_" + std::to_string(sign_symm_sector.first) + "_" +
                                           std::to_string(i) + std::to_string(j),
                                       Domain::inPSDCone(matrix_dimension));
          if (maximize_)
            Xs_[sign_symm_sector.first][i].push_back(Expr::neg(X));
          else
            Xs_[sign_symm_sector.first][i].push_back(X);
        }
      }
    }


    int i = 0;
    for (auto  psd_mat : this->sigmas_ )
    {
    //  std::cout<< "psd sites "<<psd_mat.first.op_.size()<<std::endl;
      Lambdas_.insert({psd_mat.first, {}});
      std::cout<<"first "<<std::endl;
      for(auto& elements: psd_mat.second)
      {
//     // for (auto& state : elements)
//     // {
     auto it=elements.begin();
      auto rows=it->second->numRows();
      std::cout<<"rows "<< rows<<std::endl;


     
// {
      auto beta = this->M_->variable("betas_" + std::to_string(i), Domain::inPSDCone(rows));
      i++;
      if (maximize_)
      {
        Lambdas_[psd_mat.first].push_back(Expr::neg(beta));
      }
      else
      {
        Lambdas_[psd_mat.first].push_back(beta);
      }
     
//     //}
  


//   }
 }
  }
  }

  void update_constrains()
  {
    auto totalvec = Expr::add(A_vector, Lamba_vector);
    if (this->nr_of_linear_constraints > 0)
    {
      LC_vector = Expr::mul(this->Psp, linear_constraints_variable2_);
      totalvec = Expr::add(totalvec, LC_vector);
    }
    if (this->bounding_observable_)
    {
      auto exp_temporary =
          Expr::mul(Expr::add(energy_bouding_variables_[0], energy_bouding_variables_[1]), this->energy_vec_);
      totalvec = Expr::add(totalvec, exp_temporary);
    }
    final_constraint_->update(Expr::add(Expr::add(totalvec, epsilon_vec_flat), this->b_));
  }

  void fix_constrains()
  {
    if (this->bounding_observable_)
    {
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

    if (this->nr_of_linear_constraints > 0)
      linear_constraints_variable2_ = this->M_->variable(this->nr_of_linear_constraints);

    const int n_constraints = static_cast<int>(this->lattice_.variable_map_.size());
    const int nrblocks_y = nrblocks_y_for(this->lattice_);
    const int nrblocks_x = nrblocks_x_for(this->lattice_);

    for (auto &sign_symm_sector : this->sectors_)
    {
      for (int i = 0; i < nrblocks_y; i++)
      {
        for (int j = 0; j < nrblocks_x; j++)
        {
          int block_size = 2 * sign_symm_sector.second.block_shifts[i][j];
          if (block_size == 0)
            continue;

          int n_vars_block = block_size * block_size;
          auto x_block = Expr::reshape(Xs_[sign_symm_sector.first][i][j], n_vars_block);

          std::vector<int> rows_b, cols_b;
          std::vector<double> vals_b;

          for (auto &op : this->lattice_.variable_map_)
          {
            if (op.first == "0")
              continue;
            auto &A_block = this->As_[op.first][sign_symm_sector.first][i][j];
            if (!A_block.has_elements_)
              continue;

            int el = op.second;
            auto mat = A_block.make_matrix(block_size, block_size);
            auto data = mat->getDataAsArray();

            for (int r = 0; r < block_size; r++)
              for (int c = 0; c < block_size; c++)
              {
                double v = (*data)[r * block_size + c];
                if (std::abs(v) > 1e-15)
                {
                  rows_b.push_back(el);
                  cols_b.push_back(r * block_size + c);
                  vals_b.push_back(v);
                }
              }
          }

          if (vals_b.empty())
            continue;

          auto A_block_sparse = Matrix::sparse(
              n_constraints, n_vars_block, monty::new_array_ptr(rows_b),
              monty::new_array_ptr(cols_b), monty::new_array_ptr(vals_b));

          auto contribution = Expr::mul(A_block_sparse, x_block);
          if (A_vector == nullptr)
            A_vector = contribution;
          else
            A_vector = Expr::add(A_vector, contribution);
        }
      }
    }

    for (auto &sign_symm_sector : this->sectors_)
    {
      for (int i = 1; i < int(this->lattice_.Ly_ / 2); i++)
      {
        for (int j = 1; j < int(this->lattice_.Lx_ / 2); j++)
        {
          int block_size = 2 * sign_symm_sector.second.block_shifts[i][j];
          if (block_size == 0)
            continue;

          const int mir_y = this->lattice_.Ly_ - i;
          const int mir_x = this->lattice_.Lx_ - j;
          auto &eta_ij = linear_constraints_for_block_equality_variable_[i][j];

          for (auto &op : this->lattice_.variable_map_)
          {
            if (op.first == "0")
              continue;
            auto &A_ij = this->As_[op.first][sign_symm_sector.first][i][j];
            if (!A_ij.has_elements_)
              continue;

            int el = op.second;
            auto mat_ij = A_ij.make_matrix(block_size, block_size);
            auto mat_mir = this->As_[op.first][sign_symm_sector.first][mir_y][mir_x].make_matrix(
                block_size, block_size);
            double tr = matrix_trace(mat_ij, block_size) - matrix_trace(mat_mir, block_size);
            if (std::abs(tr) < 1e-15)
              continue;

            auto e_vec = Matrix::sparse(
                n_constraints, 1, monty::new_array_ptr(std::vector<int>{el}),
                monty::new_array_ptr(std::vector<int>{0}),
                monty::new_array_ptr(std::vector<double>{tr}));
            auto term = Expr::mul(e_vec, eta_ij);
            if (A_vector == nullptr)
              A_vector = term;
            else
              A_vector = Expr::add(A_vector, term);
          }
        }
      }
    }

    int i=0;

    for (auto &[key, lambda_vec] : Lambdas_)
{
  int ll=0;
for(auto& lambda_expr: lambda_vec)
{
    int block_size = (int)std::round(std::sqrt(lambda_expr->getSize()));
    int n_vars_block = block_size * block_size;
    auto l_block = Expr::reshape(lambda_expr, n_vars_block);
  //std::cout<< "matrix size "<<block_size <<std::endl;
    std::vector<int>    rows_b, cols_b;
    std::vector<double> vals_b;
//for(auto& elements: this->sigmas_[key])
auto elements= this->sigmas_[key][ll];

    for (auto& [op_string, mat] : elements)
    {
        if (op_string == "1") continue;

        int el    = this->lattice_.variable_map_.at(op_string);
        auto data = mat->getDataAsArray();

        for (int r = 0; r < block_size; r++)
            for (int c = 0; c < block_size; c++)
            {
                double v_entry = (*data)[r * block_size + c];
                if (std::abs(v_entry) > 1e-15)
                {
                    rows_b.push_back(el);
                    cols_b.push_back(r * block_size + c);
                    vals_b.push_back(v_entry);
                }
            }
    }
  
i++;
    if (vals_b.empty()) continue;

    auto S_block_sparse = Matrix::sparse(
        n_constraints, n_vars_block,
        monty::new_array_ptr(rows_b),
        monty::new_array_ptr(cols_b),
        monty::new_array_ptr(vals_b)
    );

    auto contribution = Expr::mul(S_block_sparse, l_block);
    if (Lamba_vector == nullptr)
      Lamba_vector = contribution;
    else
      Lamba_vector = Expr::add(Lamba_vector, contribution);

ll++;
    }
}
    int el = this->lattice_.variable_map_.at("1");
    auto e_vec = Matrix::sparse(
        n_constraints, 1, monty::new_array_ptr(std::vector<int>{el}),
        monty::new_array_ptr(std::vector<int>{0}),
        monty::new_array_ptr(std::vector<double>{1.0}));
    auto epsilon_vec = Expr::mul(e_vec, epsilon);
    epsilon_vec_flat = Expr::reshape(epsilon_vec, n_constraints);

    auto totalvec = Expr::add(A_vector, Lamba_vector);
    if (this->nr_of_linear_constraints > 0)
    {
      LC_vector = Expr::mul(this->Psp, linear_constraints_variable2_);
      totalvec = Expr::add(totalvec, LC_vector);
    }
    if (this->bounding_observable_)
    {
      auto exp_temporary =
          Expr::mul(Expr::add(energy_bouding_variables_[0], energy_bouding_variables_[1]), this->energy_vec_);
      totalvec = Expr::add(totalvec, exp_temporary);
    }

    final_constraint_ = this->M_->constraint(
        Expr::add(Expr::add(totalvec, epsilon_vec_flat), this->b_), Domain::equalsTo(0.));

    return;
  }
  Expression::t get_costfunction()
  {
    Expression::t ee = Expr::constTerm(0.);
    ee = Expr::add(ee, Expr::neg(epsilon));

    // Adding matrices for the positive definite constrain of the RDMs
    for (auto lambda_ : Lambdas_)
    {
      for (int i=0; i<lambda_.second.size(); i++)
      {
  
      ee = Expr::add(ee, Expr::dot(lambda_.second[i], this->sigmas_[lambda_.first][i]["1"]));
    }
  }
    if (this->bounding_observable_)
    {

      ee = Expr::add(ee, Expr::neg(Expr::mul(this->energy_bounds_["E_upper"], this->energy_bouding_variables_[0])));
      ee = Expr::add(ee, Expr::neg(Expr::mul(this->energy_bounds_["E_lower"], this->energy_bouding_variables_[1])));
    }
    return ee;
  }
};