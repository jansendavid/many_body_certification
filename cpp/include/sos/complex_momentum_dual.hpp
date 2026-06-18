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
#include <chrono>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
using namespace mosek::fusion;
using namespace monty;

using symmetry_sector = std::map<int, std::vector<std::vector<matrix_organizer>>>;

template <typename Lattice>
class momentum_block
{
public:
  int sign_sector_{0};
  std::vector<std::vector<int>> block_shifts;
  Lattice &lattice_;

  Eigen::MatrixXcd &FTx_;
  Eigen::MatrixXcd &FTy_;

  momentum_block(Lattice &lattice, int sign_sector, Eigen::MatrixXcd &FTx, Eigen::MatrixXcd &FTy)
      : lattice_(lattice), sign_sector_(sign_sector), FTy_(FTy), FTx_(FTx)
  {
  }
  void initialize_blocks_zero(std::map<std::string, symmetry_sector> &As)
  {

    int dim_0 = lattice_.states_[sign_sector_].size() + 1;
    int dim_x = lattice_.states_[sign_sector_].size();
    block_shifts.push_back({});

    // initializing block shifts
    block_shifts[0].push_back(dim_0);
    for (int i = 1; i < lattice_.Ly_; i++)
    {

      block_shifts[0].push_back(dim_x);
    }
   
    for (int i = 1; i < lattice_.Lx_; i++)
    {
      block_shifts.push_back({});
      for (int j = 0; j < lattice_.Ly_; j++)
      {
        block_shifts[i].push_back(dim_x);
      }
    }
 
    As["1"][sign_sector_][0][0].add_values({0, 0}, 1. / 2);
    As["1"][sign_sector_][0][0].add_values({dim_0, dim_0}, 1. / 2);

    int i = 0;
   
    for (auto it = lattice_.states_[sign_sector_].begin(); it != lattice_.states_[sign_sector_].end(); ++it)
    {
      auto op = *it;

      // get normal form
      auto [coeff, nf] = get_normal_form(op);
      // get translation invariant representation

      auto ti_key = op_key_label(lattice_.TI_map_.at(key_dir_pos(nf)).first);

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
 
    return;
  }
  void initialize_blocks_general()
  {

    int dim_x = lattice_.states_[sign_sector_].size();

    for (int j = 0; j < lattice_.Lx_; j++)
    {
      block_shifts.push_back({});

      for (int i = 0; i < lattice_.Ly_; i++)
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
    const int Ly = lattice_.Ly_;
    const int Lx = lattice_.Lx_;
    std::cout<<"start generate block"<<std::endl;
    const int n_states = static_cast<int>(lattice_.states_[sign_sector_].size());
    int i = 0;
    for (auto it1 = lattice_.states_[sign_sector_].begin(); it1 != lattice_.states_[sign_sector_].end(); ++it1)
    {
      int j = i;
      for (auto it2 = it1; it2 != lattice_.states_[sign_sector_].end(); ++it2)
      {
        // generate_G_element_sos depends only on (it1,it2,pos_y,pos_x), not on mat_pos.
        for (int pos_y = 0; pos_y < Ly; ++pos_y)
        {
          for (int pos_x = 0; pos_x < Lx; ++pos_x)
          {
            const auto construct = lattice_.generate_G_element_sos(*it1, *it2, pos_x, pos_y);
            if (construct.op_ == "0")
              continue;

            for (int mat_pos_y = 0; mat_pos_y < Ly; ++mat_pos_y)
            {
              const std::complex<double> ft_y = FTy_(pos_y, mat_pos_y);

              for (int mat_pos_x = 0; mat_pos_x < Lx; ++mat_pos_x)
              {
                const int shift = block_shifts[mat_pos_x][mat_pos_y] % n_states;
                const int dim = block_shifts[mat_pos_x][mat_pos_y];
                const int ii = i + shift;
                const int jj = j + shift;
                const int ii_dim = ii + dim;
                const int jj_dim = jj + dim;

                const std::complex<double> total_prefactor =
                    construct.prefac_ * FTx_(pos_x, mat_pos_x) * ft_y;

                auto &cell = As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y];

                const double re = 0.5 * total_prefactor.real();
                if (std::abs(re) > 1e-9)
                {
                  cell.add_values({ii, jj}, re);
                  cell.add_values({ii_dim, jj_dim}, re);
                  if (i != j)
                  {
                    cell.add_values({jj, ii}, re);
                    cell.add_values({jj_dim, ii_dim}, re);
                  }
                }

                const double im = 0.5 * total_prefactor.imag();
                if (std::abs(im) > 1e-9)
                {
                  cell.add_values({ii, jj + dim}, -im);
                  cell.add_values({jj, ii + dim}, im);
                  if (i != j)
                  {
                    cell.add_values({ii_dim, jj}, im);
                    cell.add_values({jj_dim, ii}, -im);
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
std::cout<< "out generating block "<<std::endl;
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
 
  //Parameter::t P;
  Matrix::t Psp;
  int nr_of_linear_constraints{0};
  bool U1;
  momentum_basis(Lattice &lattice, Model::t M, rdms_struct rdms, bool U1=false) : lattice_(lattice), M_(M), U1(U1)
  {
    std::cout << "start" << std::endl;
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
    std::cout << "start initializeing blocks" << std::endl;
    for (auto it = lattice_.states_.begin(); it != lattice_.states_.end(); ++it)
    {

      auto Block = momentum_block(lattice_, it->first, FTx_, FTy_);
      sectors_.insert({it->first, Block});
    }
    std::cout << "start initializeing maps" << std::endl;
    initialize_all_maps(rdms);

    std::cout << "size TI map " << lattice_.TI_map_.size() << std::endl;
    std::cout << "size total refs " << lattice.variable_map_.size() << std::endl;

    for (auto it = lattice.variable_map_.begin(); it != lattice.variable_map_.end(); it++)
    {
      As_.insert({it->first, symmetry_sector()});
      for (auto it_sign_sector = sectors_.begin(); it_sign_sector != sectors_.end(); ++it_sign_sector)
      {
        As_[it->first][it_sign_sector->first] = {};

        for (int i = 0; i < lattice_.Lx_; i++)
        {
          As_[it->first][it_sign_sector->first].push_back({});
          for (int j = 0; j < lattice_.Ly_; j++)
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
    std::cout<< "started generating rdm map"<<std::endl;
    if (rdms.size() > 0)
    {

      generate_rdms(rdms);
    }
    std::cout<< "finished generating rdm map"<<std::endl;
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

  void set_linear_constraints_vec(std::set<std::vector<double>> linear_constraints)
  {
    std::cout<< "start "<<nr_of_linear_constraints<<std::endl;
    if(nr_of_linear_constraints<1)
    {
      nr_of_linear_constraints=linear_constraints.size();
      auto shape = monty::new_array_ptr<int>({
        static_cast<int>(lattice_.variable_map_.size()),
        static_cast<int>(linear_constraints.size())
    });
    
   // P = M_->parameter(shape);
    }
{
  
  const int m = static_cast<int>(lattice_.variable_map_.size());
const int n = nr_of_linear_constraints;
std::cout<< "m,n="<<m<<","<<n<<std::endl;
// std::vector<double> flat(m * n, 0.0);
  std::vector<int> rows;
  std::vector<int> cols;
  std::vector<double> vals;
int i=0;
for(auto it=linear_constraints.begin(); it!=linear_constraints.end(); ++it)
  {  for (int j = 0; j < m; ++j)
  {
        // flat[j * n + i] =(*it)[j];
        double v =(*it)[j];

        if (std::abs(v) > 1e-15)
          {
              rows.push_back(j);
              cols.push_back(i);
              vals.push_back(v);
          }
  }
  i++;
}
std::cout<<"end "<<std::endl;
Psp=Matrix::t(Matrix::sparse(
        m,
        n,
        monty::new_array_ptr<int>(rows),
        monty::new_array_ptr<int>(cols),
        monty::new_array_ptr<double>(vals)));
        std::cout<<"mape P "<<std::endl;
//P->setValue(monty::new_array_ptr<double>(flat));
// Eigen::SparseMatrix<double> A(m, n);

// std::vector<Eigen::Triplet<double>> triplets;
// for (size_t k = 0; k < vals.size(); ++k)
// {
//     triplets.emplace_back(rows[k], cols[k], vals[k]);
// }

// A.setFromTriplets(triplets.begin(), triplets.end());

// Eigen::SparseQR<
//     Eigen::SparseMatrix<double>,
//     Eigen::COLAMDOrdering<int>
// > qr;

// qr.compute(A);

// if (qr.info() != Eigen::Success)
// {
//     std::cerr << "Factorization failed\n";
// }

// int r = qr.rank();
// int full = std::min(m, n);
// std::cout<< "size "<<m << " and "<<n <<std::endl;
// std::cout << "rank = " << r << "\n";
// std::cout << "full rank = " << full << "\n";
// std::cout << "is full rank = "
//           << (r == full ? "true" : "false")
//           << std::endl;
  }
    return;
  }
  void generate_rdms(rdms_struct rdms)
  {
    auto offset = lattice_.states_[0][0][0].offset_;

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
class momentum_symmetry_solver_dual : public momentum_basis<Lattice>
{
public:
  Variable::t y_;
  momentum_symmetry_solver_dual(Lattice &lattice, Model::t M, rdms_struct rdms, bool U1=false) : momentum_basis<Lattice>(lattice, M, rdms, U1)
  {
    std::cout << "start " << std::endl;
    y_ = this->M_->variable("T", this->lattice_.variable_map_.size());
    this->M_->constraint(y_,  Domain::lessThan(1.0));
    this->M_->constraint(y_, Domain::greaterThan(-1.0));
  
    // fix 1
    std::cout << "start fix one" << std::endl;
    auto el = this->lattice_.variable_map_.at("1");
    this->M_->constraint(y_->index(el), Domain::equalsTo(1.0));
    std::cout << "start fix zero " << std::endl;
    // fix zero
    auto it=this->lattice_.variable_map_.find("0");
    if(it!=this->lattice_.variable_map_.end())
    {
    el = this->lattice_.variable_map_.at("0");
    this->M_->constraint(y_->index(el), Domain::equalsTo(0.0));
    std::cout << "start " << std::endl;
    }
  }
  void fix_constrains()
  {

    // iterate over sign sector
    for (auto sign_symm_sector : this->sectors_)
    {

      for (int i = 0; i < this->lattice_.Lx_; i++)
      {
        for (int j = 0; j < this->lattice_.Ly_; j++)
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
    // if (this->bounding_observable_)
    // {
    //   // std::cout << "introduing bounds" << std::endl;
    //   // std::cout << " upper " << this->energy_bounds_["E_upper"]->index(0) << std::endl;
    //   // std::cout << " lower " << this->energy_bounds_["E_lower"]->index(0) << std::endl;
    //   // this->M_->constraint(Expr::dot(this->energy_vec_, y_), Domain::lessThan(this->energy_bounds_["E_upper"]->index(0)));
    //   // this->M_->constraint(Expr::dot(this->energy_vec_, y_), Domain::greaterThan(this->energy_bounds_["E_lower"]->index(0)));
    // }
    return;
  }
  Expression::t get_costfunction()
  {

    return Expr::dot(this->b_, y_);
  }
};
template <typename Lattice>
class momentum_symmetry_solver_sos : public momentum_basis<Lattice>
{
public:
  std::map<int, std::vector<std::vector<Expression::t>>> Xs_;
  // Lambdas are the Lagrangian stemmeing from psd density matrices
  std::map<rdm_operator, std::vector<Expression::t>> Lambdas_;

  // variables introduced to bound the energy
  std::vector<Variable::t> energy_bouding_variables_;
  bool maximize_{true}; // if cost function is a maximization problem
  Variable::t epsilon;
  Variable::t linear_constraints_variable2_;
  Constraint::t final_constraint_;
  Expression::t A_vector=nullptr;
  Expression::t Lamba_vector=nullptr;
  Expression::t LC_vector=nullptr;
  Expression::t epsilon_vec_flat=nullptr;

  momentum_symmetry_solver_sos(Lattice &lattice, Model::t M, rdms_struct rdms, bool maximize = true, bool U1=false) : maximize_(maximize), momentum_basis<Lattice>(lattice, M, rdms, U1)
  {
    epsilon = this->M_->variable("epsilon");
    for (auto sign_symm_sector : this->sectors_)
    {
      Xs_[sign_symm_sector.first] = {};

      for (int i = 0; i < this->lattice_.Lx_; i++)
      {
        Xs_[sign_symm_sector.first].push_back({});

        for (int j = 0; j < this->lattice_.Ly_; j++)
        {

          int matrix_dimension = 2 * sign_symm_sector.second.block_shifts[i][j];

          auto X = this->M_->variable("X_" + std::to_string(sign_symm_sector.first) + "_" + std::to_string(i) + std::to_string(j), Domain::inPSDCone(matrix_dimension));

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
  auto totalvec=Expr::add(A_vector,Lamba_vector);
  if(this->nr_of_linear_constraints > 0)
  {
    LC_vector=Expr::mul(this->Psp, linear_constraints_variable2_);
    totalvec=Expr::add(totalvec, LC_vector);
  }
  final_constraint_->update(Expr::add(Expr::add(totalvec, epsilon_vec_flat), this->b_));
 
 }
  void fix_constrains()
  {
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

    if (this->nr_of_linear_constraints > 0)
    {
      linear_constraints_variable2_ = this->M_->variable(this->nr_of_linear_constraints);
    }

    int n_constraints = this->lattice_.variable_map_.size();

    for (auto &sign_symm_sector : this->sectors_)
{
    for (int i = 0; i < this->lattice_.Lx_; i++)
    {
        for (int j = 0; j < this->lattice_.Ly_; j++)
        {
            int block_size = 2 * sign_symm_sector.second.block_shifts[i][j];
            if (block_size == 0) continue;

            int n_vars_block = block_size * block_size;
            auto x_block = Expr::reshape(Xs_[sign_symm_sector.first][i][j], n_vars_block);

            std::vector<int>    rows_b, cols_b;
            std::vector<double> vals_b;

            for (auto& op : this->lattice_.variable_map_)
            {
                if (op.first == "0") continue;
                auto& A_block = this->As_[op.first][sign_symm_sector.first][i][j];
                if (!A_block.has_elements_) continue;

                int el    = op.second;
                auto mat  = A_block.make_matrix(block_size, block_size);
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

            if (vals_b.empty()) continue;

            auto A_block_sparse = Matrix::sparse(
                n_constraints, n_vars_block,
                monty::new_array_ptr(rows_b),
                monty::new_array_ptr(cols_b),
                monty::new_array_ptr(vals_b)
            );

            auto contribution = Expr::mul(A_block_sparse, x_block);
            if (A_vector == nullptr)
              A_vector = contribution;
            else
              A_vector = Expr::add(A_vector, contribution);
        }
    }
}

    std::cout << "start generating constarins for rdms " << std::endl;
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

// sparse vector with 1.0 at position el, 0 elsewhere
auto e_vec = Matrix::sparse(
    n_constraints, 1,
    monty::new_array_ptr(std::vector<int>{el}),
    monty::new_array_ptr(std::vector<int>{0}),
    monty::new_array_ptr(std::vector<double>{1.0})
);

auto epsilon_vec = Expr::mul(e_vec, epsilon);  // shape [n_constraints, 1]
epsilon_vec_flat = Expr::reshape(epsilon_vec, n_constraints);  // shape [n_constraints]

  const std::size_t vm_total = this->lattice_.variable_map_.size();
  std::cout << "A_vector     = " << (A_vector != nullptr) << std::endl;
std::cout << "Lamba_vector = " << (Lamba_vector != nullptr) << std::endl;
  std::cout << "equality constraints: " << vm_total << " variables" << std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  auto totalvec=A_vector;
  if(Lamba_vector!=nullptr)
  { totalvec = Expr::add(A_vector, Lamba_vector);}
  if (this->nr_of_linear_constraints > 0)
  {
    std::cout<< "apply linear constarints "<<std::endl;
    LC_vector = Expr::mul(this->Psp, linear_constraints_variable2_);
    totalvec = Expr::add(totalvec, LC_vector);
  }
  final_constraint_=this->M_->constraint(Expr::add(Expr::add(totalvec, epsilon_vec_flat), this->b_),
           Domain::equalsTo(0.));
        
  if (vm_total > 0)
    std::cout << std::endl;
    std::cout << "Finished generating the PSD constraints ones " << std::endl;
    auto end = std::chrono::high_resolution_clock::now();

auto duration =
    std::chrono::duration_cast<std::chrono::milliseconds>(
        end - start);

std::cout << "Time: "
          << duration.count()
          << " ms"
          << std::endl;
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
