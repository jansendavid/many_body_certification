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
  std::vector<op_vec> operators_;
  Eigen::MatrixXcd &FTx_;
  Eigen::MatrixXcd &FTy_;

  std::map<std::string, int> &total_refs_;
  momentum_block(Lattice &lattice, std::vector<op_vec> operators, Model::t M, int sign_sector, std::map<std::string, int> &total_refs, Eigen::MatrixXcd &FTy, Eigen::MatrixXcd &FTx, std::string sector_label = "") : lattice_(lattice), operators_(operators), sign_sector_(sign_sector), total_refs_(total_refs), FTy_(FTy), FTx_(FTx)
  {
    // std::cout << FTx_ << std::endl;
    // std::cout << FTy_ << std::endl;
  }

  void initialize_blocks_general()
  {

    int dim_x = operators_.size(); // dimension of other blocks

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
  void initialize_blocks_zero(std::map<std::string, symmetry_sector> &As)
  {

    int dim_0 = operators_.size() + 1; // dimension of 0th block
    int dim_x = operators_.size();     // dimension of other blocks

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

    for (auto it = operators_.begin(); it != operators_.end(); ++it)
    {
      auto op = *it;
      // get normal form
      auto [coeff, nf] = get_normal_form(op);
      // get translation invariant representation

      auto ti_key = lattice_.TI_map_.at(print_op(nf)).first;

      auto el = total_refs_.at(ti_key);

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
  void generate_block(std::map<std::string, symmetry_sector> &As)
  {
    //     const auto start{std::chrono::steady_clock::now()};

    int i = 0;
    for (auto it1 = operators_.begin(); it1 != operators_.end(); ++it1)
    {
      int j = i;
      for (auto it2 = it1; it2 != operators_.end(); ++it2)
      {

        for (int mat_pos_y = 0; mat_pos_y < lattice_.Ly_; mat_pos_y++)
        {
          for (int mat_pos_x = 0; mat_pos_x < lattice_.Lx_; mat_pos_x++)
          {

            // 			      // determines if first block of zeroth moment blocks
            int shift = block_shifts[mat_pos_y][mat_pos_x] % operators_.size();

            // 			      // gives the shift between real and complex components
            int dim = block_shifts[mat_pos_y][mat_pos_x];

            for (int pos_y = 0; pos_y < lattice_.Ly_; pos_y++)
            {
              std::complex<double> FT_factor_y = FTy_(pos_y, mat_pos_y);

              for (int pos_x = 0; pos_x < lattice_.Lx_; pos_x++)
              {
                std::complex<double> FT_factor_x = FTx_(pos_x, mat_pos_x);

                //              // to do, correct so that all terms appearing here appear in map
                auto construct = lattice_.generate_G_element_sos(*it1, *it2, pos_y, pos_x, lattice_.TI_map_);
                std::complex<double> total_prefactor = construct.prefac_ * FT_factor_x * FT_factor_y;
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
  basis_structure operators_;
  Model::t M_;
  std::map<int, momentum_block<Lattice>> sectors_;
  std::string sector_;

  std::map<std::string, int> total_refs_;
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

  momentum_basis(Lattice &lattice, basis_structure operators, Model::t M, rdms_struct rdms) : lattice_(lattice), operators_(operators), M_(M)
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

    for (auto it = operators.begin(); it != operators.end(); ++it)
    {

      auto Block = momentum_block(lattice_, it->second, M_, it->first, total_refs_, FTy_, FTx_, std::to_string(it->first));
      sectors_.insert({it->first, Block});
    }

    initialize_all_maps(rdms);

    std::cout << "size TI map " << lattice_.TI_map_.size() << std::endl;
    std::cout << "size total refs " << total_refs_.size() << std::endl;

    for (auto it = total_refs_.begin(); it != total_refs_.end(); it++)
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
    {
      sector.second.initialize_blocks(As_);
    }

    for (auto it_2 = sectors_.begin(); it_2 != sectors_.end(); ++it_2)
    {
      it_2->second.generate_block(As_);
    }
    std::cout << "finished making the As matrices" << std::endl;

    return;
  };
  void initialize_all_maps(rdms_struct rdms)
  {
    std::map<std::string, op_vec> mat_terms;
    for (auto &b : sectors_)
    {
      lattice_.generate_TI_map(mat_terms, b.second.operators_, b.first);

      std::cout << "sizes " << mat_terms.size() << " and " << lattice_.TI_map_.size() << std::endl;
    }
    // for (auto a : lattice_.TI_map_)
    // {
    //   if (a.first.size() == 11 * 2)
    //   {
    //     std::cout << a.first << std::endl;
    //   }
    // }

    auto size_without_rdms = mat_terms.size();
    if (rdms.size() > 0)
    {

      generate_rdms(rdms, mat_terms);
    }
    assert(size_without_rdms == mat_terms.size());
    std::cout << "sizes after rdm " << mat_terms.size() << " and " << lattice_.TI_map_.size() << std::endl;

    int new_index = 0;

    for (auto a : mat_terms)
    {

      total_refs_.insert({a.first, new_index});

      new_index += 1;
    }

    b_ = M_->parameter("b", total_refs_.size());
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
      energy_vec_ = M_->parameter("energy vec", total_refs_.size());
      energy_bounds_["E_upper"] = M_->parameter("E_upper");
      energy_bounds_["E_lower"] = M_->parameter("E_lower");
    }

    auto a = monty::new_array_ptr<double>(energy_vec);
    energy_vec_->setValue(a);

    std::cout << "bounding " << bounding_observable_ << std::endl;
    energy_bounds_["E_upper"]->setValue(E_upper);
    energy_bounds_["E_lower"]->setValue(E_lower);
    return;
  }

  std::map<std::string, Matrix::t> generate_rdms_primal_cp(rdm_operator sites, std::vector<int> offset, std::map<std::string, op_vec> &mat_terms) //,std::map<std::string, int> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, Variable::t var,int Lx)
  {

    std::map<std::string, Matrix::t> sigmas_temp_;
    std::map<std::string, mat_type> rdms_eigen_;
    mat_type pauliI = mat_type::Zero(2, 2);
    pauliI(0, 0) = 1;
    pauliI(1, 1) = 1;

    // pauliI.makeCompressed();
    mat_type pauliZ = mat_type::Zero(2, 2);
    pauliZ(0, 0) = 1;
    pauliZ(1, 1) = -1;

    // pauliZ.makeCompressed();
    mat_type pauliX = mat_type::Zero(2, 2);
    pauliX(0, 1) = 1;
    pauliX(1, 0) = 1;

    mat_type pauliY = mat_type::Zero(2, 2);
    pauliY(0, 1) = std::complex<double>(0, -1);
    pauliY(1, 0) = std::complex<double>(0, 1);
    int degree = sites.size();
    std::vector<std::string> terms;
    std::map<std::string, mat_type> sigma_map;

    sigma_map.insert({"1", pauliI});
    sigma_map.insert({"x", pauliX});
    sigma_map.insert({"y", pauliY});
    sigma_map.insert({"z", pauliZ});

    auto dirs = std::vector<std::string>{"1", "x", "y", "z"};
    std::set<std::string> tots;

    for (auto d1 : dirs)
    {

      if (degree == 1)
      {
        tots.insert(d1);
        continue;
      }
      for (auto d2 : dirs)
      {
        if (degree == 2)
        {
          tots.insert(d1 + d2);
          continue;
        }
        for (auto d3 : dirs)
        {
          if (degree == 3)
          {
            tots.insert(d1 + d2 + d3);
            continue;
          }
          for (auto d4 : dirs)
          {
            if (degree == 4)
            {
              tots.insert(d1 + d2 + d3 + d4);
              continue;
            }
            for (auto d5 : dirs)
            {
              if (degree == 5)
              {
                tots.insert(d1 + d2 + d3 + d4 + d5);
                continue;
              }
              for (auto d6 : dirs)
              {
                if (degree == 6)
                {
                  tots.insert(d1 + d2 + d3 + d4 + d5 + d6);
                  continue;
                }
                for (auto d7 : dirs)
                {
                  if (degree == 7)
                  {
                    tots.insert(d1 + d2 + d3 + d4 + d5 + d6 + d7);
                    continue;
                  }
                  for (auto d8 : dirs)
                  {
                    if (degree == 8)
                    {
                      tots.insert(d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8);
                      continue;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    int dim = std::pow(2, degree);

    std::vector<Expression::t> matrices;

    double prefac = 1;

    for (auto t : tots)
    {

      mat_type mat;
      op_vec state;

      for (int i = 0; i < t.size(); i++)
      {

        std::string key = t.substr(i, 1);

        if (key != "1")
        {
          state.push_back(spin_op(key, sites.at(i), offset));
        }
        if (i == 0)
        {
          mat = sigma_map[key];
        }
        else
        {

          mat = Eigen::KroneckerProduct(mat, sigma_map[key]).eval();
        }
      }

      // if matrix element exists I only
      if (print_op(state) == "1")
      {

        mat = mat / std::
                        pow(2, degree);
        if (rdms_eigen_.find("1") != rdms_eigen_.end())
        {
          rdms_eigen_["1"] += mat;
        }
        else
        {
          rdms_eigen_.insert({"1", mat});
        }
      }
      else
      {

        auto [fac, nf] = get_normal_form(state);

        auto [found, op_string] = lattice_.check_if_operator_exists(nf, mat_terms, true);
        // call error if operator does not exist, we do not add new

        assert(found);
        lattice_.TI_map_.insert({print_op(nf), {op_string, 1}});
        auto [state_from_map, coeff] = lattice_.TI_map_.at(print_op(nf));
        if (state_from_map == "0")
        {
        }
        else
        {

          assert(std::abs((fac * coeff).imag()) < 1e-9);
          mat = mat * (fac * coeff).real() / std::pow(2, degree);

          if (rdms_eigen_.find(state_from_map) != rdms_eigen_.end())
          {

            rdms_eigen_[state_from_map] += mat;
          }
          else
          {

            rdms_eigen_.insert({state_from_map, mat});
          }
        }
      }
    }
    // convert to mosek format
    for (auto eigen_matrix : rdms_eigen_)
    {
      auto Alpha = get_sparse_from_eigen(eigen_matrix.second);

      sigmas_temp_.insert({eigen_matrix.first, Alpha});
    }
    return sigmas_temp_;
  }
  void generate_rdms(rdms_struct rdms, std::map<std::string, op_vec> &mat_terms)
  {

    auto offset = operators_[0][0][0].offset_; // change this to be derived from baso
    int i = 0;
    std::cout << "rdms size " << rdms.rdms.size() << std::endl;
    for (auto site : rdms.rdms)
    {

      i++;
      auto sigmas_temp = generate_rdms_primal_cp(site, offset, mat_terms);
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
  momentum_symmetry_solver_dual(Lattice &lattice, basis_structure operators, Model::t M, rdms_struct rdms) : momentum_basis<Lattice>(lattice, operators, M, rdms)
  {
    y_ = this->M_->variable("T", this->total_refs_.size());
    // fix 1
    auto el = this->total_refs_.at("1");
    this->M_->constraint(y_->index(el), Domain::equalsTo(1.0));
    // fix zero
    el = this->total_refs_.at("0");
    this->M_->constraint(y_->index(el), Domain::equalsTo(0.0));
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

          for (auto op : this->total_refs_)
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
          ee = Expr::add(ee, Expr::mul(y_->index(this->total_refs_[op_string.first]), op_string.second));
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
template <typename Lattice>
class momentum_symmetry_solver_sos : public momentum_basis<Lattice>
{
public:
  std::map<int, std::vector<std::vector<Expression::t>>> Xs_;
  // Lambdas are the Lagrangian stemmeing from psd density matrices
  std::map<rdm_operator, Expression::t> Lambdas_;

  // here we store the C matrices (the constants)
  std::map<int, std::vector<std::vector<Matrix::t>>> Cs_;
  std::map<int, std::vector<std::vector<Matrix::t>>> zeros_;
  // variables introduced to bound the energy
  std::vector<Variable::t> energy_bouding_variables_;
  bool maximize_{true}; // if cost function is a maximization problem
  momentum_symmetry_solver_sos(Lattice &lattice, basis_structure operators, Model::t M, rdms_struct rdms, bool maximize = true) : maximize_(maximize), momentum_basis<Lattice>(lattice, operators, M, rdms)
  {
    for (auto sign_symm_sector : this->sectors_)
    {
      Xs_[sign_symm_sector.first] = {};
      Cs_[sign_symm_sector.first] = {};
      zeros_[sign_symm_sector.first] = {};

      for (int i = 0; i < this->lattice_.Ly_; i++)
      {
        Xs_[sign_symm_sector.first].push_back({});
        Cs_[sign_symm_sector.first].push_back({});
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

    std::vector<Expression::t> expressions_(this->total_refs_.size(), Expr::constTerm(0));

    for (auto sign_symm_sector : this->sectors_)
    {

      for (int i = 0; i < this->lattice_.Ly_; i++)
      {
        for (int j = 0; j < this->lattice_.Lx_; j++)
        {

          for (auto op : this->total_refs_)
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
                  auto C = this->As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension, matrix_dimension);
                  Cs_[sign_symm_sector.first][i].push_back(C);
                }
              }
              else
              {
                if (this->As_[op.first][sign_symm_sector.first][i][j].has_elements_)
                {

                  int el = this->total_refs_.at(op.first);

                  expressions_[el] = Expr::add(expressions_[el], Expr::dot(this->As_[op.first][sign_symm_sector.first][i][j].make_matrix(matrix_dimension, matrix_dimension), (Xs_[sign_symm_sector.first][i][j])));
                }
              }
            }
          }
        }
      }
    }

    for (auto lambda_ : Lambdas_)
    {
      for (auto string_and_matrix : this->sigmas_[lambda_.first])
      {
        if (string_and_matrix.first != "1")
        {
          int el = this->total_refs_.at(string_and_matrix.first);
          auto a = lambda_.second;
          // why not neg (Expr::dot(lambda_.second,string_and_matrix.second )) ?
          expressions_[el] = Expr::add(expressions_[el], (Expr::dot(lambda_.second, string_and_matrix.second)));
        }
      }
    }
    if (this->bounding_observable_)
    {
      auto exp_temporary = Expr::mul(Expr::add(energy_bouding_variables_[0], energy_bouding_variables_[1]), this->energy_vec_);
      for (int i = 0; i < this->energy_vec_->getSize(); i++)
      {
        // energy_vec_->index(i)

        expressions_[i] = Expr::add(expressions_[i], (exp_temporary->index(i)));
      }
    }

    for (auto a : this->total_refs_)
    {
      if (a.first != "1" && a.first != "0")
      {

        int el = this->total_refs_.at(a.first);
        //=-1*b[el]
        this->M_->constraint(Expr::add(expressions_[el], this->b_->index(el)), Domain::equalsTo(0.));
      }
    }

    std::cout << "Finished generating the PSD constraints ones " << std::endl;
    return;
  }
  Expression::t get_costfunction()
  {
    Expression::t ee = Expr::constTerm(0.);
    for (auto sign_symm_sector : this->sectors_)
    {

      for (int i = 0; i < this->lattice_.Ly_; i++)
      {
        for (int j = 0; j < this->lattice_.Lx_; j++)
        {

          ee = Expr::add(ee, Expr::dot(Cs_[sign_symm_sector.first][i][j], (Xs_[sign_symm_sector.first][i][j])));
        }
      }
    }

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