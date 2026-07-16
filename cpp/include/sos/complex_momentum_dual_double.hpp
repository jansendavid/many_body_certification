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
class momentum_block_double
{
public:
  int sign_sector_{0};
  std::vector<std::vector<int>> block_shifts;
	std::vector<std::vector<int>> state_optimality_block_shifts;
  Lattice &lattice_;

  Eigen::MatrixXcd &FTx_;
  Eigen::MatrixXcd &FTy_;

  momentum_block_double(Lattice &lattice, int sign_sector, Eigen::MatrixXcd &FTx, Eigen::MatrixXcd &FTy)
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

    const double zero_block_scale =
        0.5 * std::sqrt(lattice_.Lx_) * std::sqrt(lattice_.Ly_);
    const auto add_zero_block_state =
        [&](const std::string &ti_key, int state_index,
            std::complex<double> coeff) {
          auto &cell = As[ti_key][sign_sector_][0][0];
          const double re = coeff.real() * zero_block_scale;
          if (std::abs(re) > 1e-9)
          {
            cell.add_values({0, state_index}, re);
            cell.add_values({state_index, 0}, re);
            cell.add_values({dim_0, state_index + dim_0}, re);
            cell.add_values({state_index + dim_0, dim_0}, re);
          }

          const double im = coeff.imag() * zero_block_scale;
          if (std::abs(im) > 1e-9)
          {
            cell.add_values({0, state_index + dim_0}, -im);
            cell.add_values({state_index + dim_0, 0}, -im);
            cell.add_values({dim_0, state_index}, im);
            cell.add_values({state_index, dim_0}, im);
          }
        };

    //   //     // The "c" terms first row and column in block 0
    int i = 0;

    for (auto it = lattice_.states_[sign_sector_][0].begin(); it != lattice_.states_[sign_sector_][0].end(); ++it)
    {
      auto op = *it;
      // get normal form
      auto [coeff, nf] = lattice_.get_form_of_TI_map(op);
      // get translation invariant representation

      auto ti_key = op_key_label(lattice_.TI_map_.at(key_dir_pos(nf)).first);

      //auto el = this->lattice_.variable_map_.at(ti_key);

      add_zero_block_state(ti_key, i + 1, coeff);

      i++;
    }
    i = 0;
    int shift=lattice_.states_[sign_sector_][0].size();
    for (auto it = lattice_.states_[sign_sector_][1].begin(); it != lattice_.states_[sign_sector_][1].end(); ++it)
    {
      auto op = *it;
      // get normal form
      auto [coeff, nf] = lattice_.get_form_of_TI_map(op);
      // get translation invariant representation

      auto ti_key = op_key_label(lattice_.TI_map_.at(key_dir_pos(nf)).first);
      //auto el = this->lattice_.variable_map_.at(ti_key);

      add_zero_block_state(ti_key, i + 1 + shift, coeff);

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

	void initialize_state_optimality_block_shifts()
	{
	  const int dimension = static_cast<int>(
		  lattice_.state_optimality_states_[sign_sector_][0].size() +
		  lattice_.state_optimality_states_[sign_sector_][1].size());
	  state_optimality_block_shifts.assign(
		  lattice_.Lx_, std::vector<int>(lattice_.Ly_, dimension));
	}
  void initialize_blocks_general()
  {


    int dim_x = lattice_.states_[sign_sector_][0].size()+lattice_.states_[sign_sector_][1].size(); // operators_.size(); // dimension of other blocks

    for (int j = 0; j < lattice_.Lx_; j++)
    {
      block_shifts.push_back({});

      for (int i = 0; i < lattice_.Ly_; i++)
      {
        block_shifts[j].push_back(dim_x);
      }
    }
  //}
    return;
  }

 template <typename OperatorVector>
 void run_loop(OperatorVector &operator_1, OperatorVector &operator_2,
               std::map<std::string, symmetry_sector> &As,
               std::complex<double> fac_orig, std::pair<int,int> shift)
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
          auto construct = lattice_.generate_G_element_sos_double(*it1, *it2, pos_x, pos_y);
          if (construct.op_ == "0")
            continue;

          for (int mat_pos_y = 0; mat_pos_y < Ly; ++mat_pos_y)
          {
            const std::complex<double> FT_factor_y = FTy_(pos_y, mat_pos_y);

            for (int mat_pos_x = 0; mat_pos_x < Lx; ++mat_pos_x)
            {
              const int shift_initial = block_shifts[mat_pos_x][mat_pos_y] % n_states;
              const int dim = block_shifts[mat_pos_x][mat_pos_y];
              const std::complex<double> FT_factor_x = FTx_(pos_x, mat_pos_x);
              const std::complex<double> total_prefactor =
                  construct.prefac_ * fac_orig * FT_factor_x * FT_factor_y;

              if (std::abs(total_prefactor.real()) > 1e-9)
              {
                As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values(
                    {i + shift.first + shift_initial, j + shift.second + shift_initial},
                    1. / 2 * total_prefactor.real());
                As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values(
                    {i + shift.first + dim + shift_initial, j + shift.second + dim + shift_initial},
                    1. / 2 * total_prefactor.real());
              }
              if (std::abs(total_prefactor.imag()) > 1e-9)
              {
                As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values(
                    {i + shift.first + shift_initial, j + shift.second + dim + shift_initial},
                    -1. / 2 * total_prefactor.imag());
                As[construct.op_][sign_sector_][mat_pos_x][mat_pos_y].add_values(
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

  void generate_state_optimality_block(
      std::map<std::string, symmetry_sector> &state_optimality_As)
  {
    if (!lattice_.state_optimality_hamiltonian_)
      return;

    const int Ly = lattice_.Ly_;
    const int Lx = lattice_.Lx_;
    const int even_dimension =
		static_cast<int>(
			lattice_.state_optimality_states_[sign_sector_][0].size());
    auto operators = lattice_.state_optimality_states_[sign_sector_][0];
    operators.insert(
        operators.end(),
        lattice_.state_optimality_states_[sign_sector_][1].begin(),
        lattice_.state_optimality_states_[sign_sector_][1].end());
    const auto phase_start = std::chrono::steady_clock::now();
    const std::size_t progress_interval =
        std::max<std::size_t>(1, operators.size() / 20);
    std::cout << "State-optimality coefficients sector " << sign_sector_
              << " start: basis=" << operators.size()
              << ", upper-triangle pairs="
              << operators.size() * (operators.size() + 1) / 2 << std::endl;

    for (std::size_t row = 0; row < operators.size(); ++row)
    {
      if (row % progress_interval == 0)
      {
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - phase_start).count();
        std::cout << "State-optimality coefficients sector " << sign_sector_
                  << ": row " << row << '/' << operators.size()
                  << ", elapsed=" << elapsed << " s, entry_cache="
                  << lattice_.state_optimality_entry_cache_.size()
                  << std::endl;
      }
      const std::complex<double> row_phase =
          row < static_cast<std::size_t>(even_dimension)
              ? std::complex<double>{1., 0.}
              : std::complex<double>{0., 1.};
      for (std::size_t col = row; col < operators.size(); ++col)
      {
        const std::complex<double> col_phase =
            col < static_cast<std::size_t>(even_dimension)
                ? std::complex<double>{1., 0.}
                : std::complex<double>{0., 1.};
        const auto sector_prefactor = std::conj(row_phase) * col_phase;

        for (int pos_y = 0; pos_y < Ly; ++pos_y)
        {
          for (int pos_x = 0; pos_x < Lx; ++pos_x)
          {
            std::map<std::string, std::complex<double>> entry;
            const auto &reduced_entry = lattice_.get_state_optimality_entry(
                operators[row], operators[col], pos_x, pos_y);
            for (const auto &[label, reduced_term] : reduced_entry.get_terms())
            {
              (void)label;
              const auto normal_key = key_dir_pos(reduced_term.get_op());
              const auto relation = lattice_.TI_map_.find(normal_key);
              if (relation == lattice_.TI_map_.end())
                throw std::logic_error(
                    "missing state-optimality TI-map entry in sector " +
                    std::to_string(sign_sector_) + ": " +
                    op_key_label(normal_key));
              const std::string moment =
                  op_key_label(relation->second.first);
              if (moment != "0")
                entry[moment] += reduced_term.get_coeff() *
                                 relation->second.second;
            }

            for (int momentum_y = 0; momentum_y < Ly; ++momentum_y)
            {
              const auto fourier_y = FTy_(pos_y, momentum_y);
              for (int momentum_x = 0; momentum_x < Lx; ++momentum_x)
              {
                const int dim = state_optimality_block_shifts
                    [momentum_x][momentum_y];
                const auto fourier =
                    FTx_(pos_x, momentum_x) * fourier_y;
                for (const auto &[moment, moment_coefficient] : entry)
                {
                  const auto total = moment_coefficient * sector_prefactor *
                                     fourier;
                  auto &cell = state_optimality_As[moment][sign_sector_]
                                                     [momentum_x][momentum_y];
                  const double re = 0.5 * total.real();
                  if (std::abs(re) > 1e-9)
                  {
                    cell.add_values(
                        {static_cast<int>(row), static_cast<int>(col)}, re);
                    cell.add_values(
                        {static_cast<int>(row) + dim,
                         static_cast<int>(col) + dim}, re);
                    if (row != col)
                    {
                      cell.add_values(
                          {static_cast<int>(col), static_cast<int>(row)}, re);
                      cell.add_values(
                          {static_cast<int>(col) + dim,
                           static_cast<int>(row) + dim}, re);
                    }
                  }

                  const double im = 0.5 * total.imag();
                  if (std::abs(im) > 1e-9)
                  {
                    cell.add_values(
                        {static_cast<int>(row),
                         static_cast<int>(col) + dim}, -im);
                    cell.add_values(
                        {static_cast<int>(col),
                         static_cast<int>(row) + dim}, im);
                    if (row != col)
                    {
                      cell.add_values(
                          {static_cast<int>(row) + dim,
                           static_cast<int>(col)}, im);
                      cell.add_values(
                          {static_cast<int>(col) + dim,
                           static_cast<int>(row)}, -im);
                    }
                  }
                }
              }
            }
          }
        }
      }
      lattice_.state_optimality_entry_cache_.clear();
      lattice_.clear_caches();
    }
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - phase_start).count();
    std::cout << "State-optimality coefficients sector " << sign_sector_
              << " complete: elapsed=" << elapsed
              << " s, entry_cache="
              << lattice_.state_optimality_entry_cache_.size() << std::endl;
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
	std::map<std::string, symmetry_sector> state_optimality_As_;
  // for the reduced density matrix
  std::map<rdm_operator, std::vector<std::map<std::string, Matrix::t>>> sigmas_;
  Matrix::t Psp;
  int nr_of_linear_constraints{0};
  bool U1=false;
	bool enable_state_optimality_conditions_{false};

  momentum_basis_double(Lattice &lattice, Model::t M, rdms_struct rdms,
						bool U1, bool enable_state_optimality_conditions)
		: lattice_(lattice), M_(M), U1(U1),
		  enable_state_optimality_conditions_(enable_state_optimality_conditions)
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
      auto Block = momentum_block_double(lattice_, it->first, FTx_, FTy_);
      sectors_.insert({it->first, Block});
    }
    initialize_all_maps(rdms);
    std::cout << "size TI map " << lattice_.TI_map_.size() << std::endl;
    std::cout << "size total refs " << lattice.variable_map_.size() << std::endl;

    for (auto it = lattice.variable_map_.begin(); it != lattice.variable_map_.end(); it++)
    {
      As_.insert({it->first, symmetry_sector()});
		if (enable_state_optimality_conditions_)
		  state_optimality_As_.insert({it->first, symmetry_sector()});
      for (auto it_sign_sector = sectors_.begin(); it_sign_sector != sectors_.end(); ++it_sign_sector)
      {
        As_[it->first][it_sign_sector->first] = {};
		if (enable_state_optimality_conditions_)
		  state_optimality_As_[it->first][it_sign_sector->first] = {};

        for (int i = 0; i < lattice_.Lx_; i++)
        {
          As_[it->first][it_sign_sector->first].push_back({});
		  if (enable_state_optimality_conditions_)
			state_optimality_As_[it->first][it_sign_sector->first].push_back({});
          for (int j = 0; j < lattice_.Ly_; j++)
          {
            As_[it->first][it_sign_sector->first][i].push_back(matrix_organizer());
			if (enable_state_optimality_conditions_)
			  state_optimality_As_[it->first][it_sign_sector->first][i]
				  .push_back(matrix_organizer());
          }
        }
      }
    }

    for (auto &sector : sectors_)
      sector.second.initialize_blocks(As_);
	if (enable_state_optimality_conditions_)
	  for (auto &sector : sectors_)
		sector.second.initialize_state_optimality_block_shifts();

    std::cout << "Coefficient phase 1/2: ordinary moment blocks" << std::endl;
    for (auto it_2 = sectors_.begin(); it_2 != sectors_.end(); ++it_2)
    {
      const auto sector_start = std::chrono::steady_clock::now();
      std::cout << "Ordinary coefficient sector " << it_2->first
                << " start" << std::endl;
      it_2->second.generate_block(As_);
      std::cout << "Ordinary coefficient sector " << it_2->first
                << " complete: elapsed="
                << std::chrono::duration<double>(
                       std::chrono::steady_clock::now() - sector_start).count()
                << " s" << std::endl;
    }

	if (enable_state_optimality_conditions_)
	{
	  std::cout << "Coefficient phase 2/2: state-optimality blocks"
				<< std::endl;
	  std::cout << "State-optimality basis maximum degree: "
				<< lattice_.state_optimality_basis_degree_ << std::endl;
	  for (const auto &[sector, block] : sectors_)
		std::cout << "State-optimality sector " << sector
				  << " basis dimension: "
				  << block.state_optimality_block_shifts[0][0] << std::endl;
	  for (auto &sector : sectors_)
		sector.second.generate_state_optimality_block(state_optimality_As_);
	  report_state_optimality_matrix_diagnostics();
	}

    return;
  };

	void report_state_optimality_matrix_diagnostics() const
	{
	  std::size_t referenced_moments = 0;
	  std::size_t nonzero_blocks = 0;
	  std::size_t stored_entries = 0;
	  double largest_symmetry_discrepancy = 0.;

	  for (const auto &[moment, by_sector] : state_optimality_As_)
	  {
		(void)moment;
		bool moment_referenced = false;
		for (const auto &[sector, by_x] : by_sector)
		{
		  (void)sector;
		  for (const auto &by_y : by_x)
			for (const auto &matrix : by_y)
			{
			  if (!matrix.has_elements_)
				continue;
			  moment_referenced = true;
			  ++nonzero_blocks;
			  stored_entries += matrix.matrix_positions.size();

			  std::map<int_pair, double> values;
			  for (std::size_t i = 0; i < matrix.matrix_positions.size(); ++i)
				values[matrix.matrix_positions[i]] += matrix.matrix_values[i];
			  for (const auto &[position, value] : values)
			  {
				const auto transposed = values.find(
					{position.second, position.first});
				const double transposed_value =
					transposed == values.end() ? 0. : transposed->second;
				largest_symmetry_discrepancy = std::max(
					largest_symmetry_discrepancy,
					std::abs(value - transposed_value));
			  }
			}
		}
		if (moment_referenced)
		  ++referenced_moments;
	  }

	  std::cout << "State-optimality coefficient matrices: moments="
				<< referenced_moments << ", nonzero blocks=" << nonzero_blocks
				<< ", stored entries=" << stored_entries
				<< ", largest realified symmetry discrepancy="
				<< largest_symmetry_discrepancy << std::endl;
	  if (largest_symmetry_discrepancy > 1e-7)
		throw std::logic_error(
			"state-optimality coefficient matrix is not Hermitian");
	}
  void initialize_all_maps(rdms_struct rdms)
  {
		this->lattice_.generate_TI_map_double(
			enable_state_optimality_conditions_);
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

    auto offset = lattice_.states_[1][0][0][0].offset_;
    std::cout << "rdms size " << rdms.rdms.size() << std::endl;
    for (auto site : rdms.rdms)
    {
      // Match the standard formulation: U1 selects fixed-number fermionic
      // blocks, while false imposes the corrected full (CP) density matrix.
      auto sigmas_temp = U1
          ? lattice_.generate_rdms_primal_U1(site, offset)
          : lattice_.generate_rdms_primal_cp(site, offset);
      sigmas_.insert({site, std::move(sigmas_temp)});
    
    }
    return;
  }
};
template <typename Lattice>
class momentum_symmetry_solver_dual_double : public momentum_basis_double<Lattice>
{
public:
  Variable::t y_;
  momentum_symmetry_solver_dual_double(
		Lattice &lattice, Model::t M, rdms_struct rdms, bool U1 = false,
		bool enable_state_optimality_conditions = false)
		: momentum_basis_double<Lattice>(
			  lattice, M, rdms, U1, enable_state_optimality_conditions)
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

	if (this->enable_state_optimality_conditions_)
	{
	  std::size_t state_optimality_psd_blocks = 0;
	  for (auto &sign_symm_sector : this->sectors_)
	  {
		for (int i = 0; i < this->lattice_.Lx_; ++i)
		{
		  for (int j = 0; j < this->lattice_.Ly_; ++j)
		  {
			const int matrix_dimension = 2 *
				sign_symm_sector.second.state_optimality_block_shifts[i][j];
			std::vector<Expression::t> matrices;
			for (const auto &op : this->lattice_.variable_map_)
			{
			  if (op.first == "0")
				continue;
			  auto &coefficient_matrix =
				  this->state_optimality_As_[op.first]
					  [sign_symm_sector.first][i][j];
			  if (coefficient_matrix.has_elements_)
				matrices.push_back(Expr::mul(
					y_->index(op.second),
					coefficient_matrix.make_matrix(
						matrix_dimension, matrix_dimension)));
			}
			if (matrices.empty())
			  continue;

			Expression::t state_optimality_matrix = matrices.front();
			for (std::size_t n = 1; n < matrices.size(); ++n)
			  state_optimality_matrix =
				  Expr::add(state_optimality_matrix, matrices[n]);
			this->M_->constraint(
				state_optimality_matrix, Domain::inPSDCone());
			++state_optimality_psd_blocks;
		  }
		}
	  }
	  std::cout << "Finished generating " << state_optimality_psd_blocks
				<< " state-optimality PSD constraints" << std::endl;
	}

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

    if (this->bounding_observable_)
    {
      auto energy = Expr::dot(this->energy_vec_, y_);
      this->M_->constraint(
          Expr::sub(energy, this->energy_bounds_["E_upper"]->index(0)),
          Domain::lessThan(0.0));
      this->M_->constraint(
          Expr::sub(energy, this->energy_bounds_["E_lower"]->index(0)),
          Domain::greaterThan(0.0));
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
	std::map<int, std::vector<std::vector<Expression::t>>>
		state_optimality_Xs_;
  std::map<rdm_operator, std::vector<Expression::t>> Lambdas_;

  std::vector<Variable::t> energy_bouding_variables_;
  bool maximize_{true};
  Variable::t epsilon;
  Variable::t upper_box_multiplier_;
  Variable::t lower_box_multiplier_;
  Variable::t zero_moment_multiplier_ = nullptr;
  Variable::t linear_constraints_variable2_;
  std::vector<std::vector<Variable::t>> linear_constraints_for_block_equality_variable_;
  Constraint::t final_constraint_;
  Expression::t A_vector = nullptr;
	Expression::t state_optimality_A_vector = nullptr;
  Expression::t Lamba_vector = nullptr;
  Expression::t LC_vector = nullptr;
  Expression::t epsilon_vec_flat = nullptr;

  static int nrblocks_y_for(const Lattice &lattice)
  {
    return lattice.Ly_;
  }

  static int nrblocks_x_for(const Lattice &lattice)
  {
    return lattice.Lx_;
  }

  static std::pair<int, int> conjugate_momentum(const Lattice &lattice,
                                                int kx, int ky)
  {
    // Time reversal/complex conjugation maps a translation momentum block
    // k=(kx,ky) to -k modulo the finite lattice periods.
    return {(lattice.Lx_ - kx) % lattice.Lx_,
            (lattice.Ly_ - ky) % lattice.Ly_};
  }

  static bool is_momentum_representative(const Lattice &lattice,
                                         int kx, int ky)
  {
    const auto [cx, cy] = conjugate_momentum(lattice, kx, ky);
    // Keep one lexicographic representative of each {k,-k} orbit.  For a
    // 4x4 lattice this gives 10 momenta per sector: 4 self-conjugate points
    // plus 6 two-point conjugate orbits.  The old rectangular 3x3 selection
    // missed the (1,3)<->(3,1) orbit.
    return std::tie(kx, ky) <= std::tie(cx, cy);
  }

  momentum_symmetry_solver_sos_double(
		Lattice &lattice, Model::t M, rdms_struct rdms,
		bool maximize = true, bool U1 = false,
		bool enable_state_optimality_conditions = false)
      : maximize_(maximize), momentum_basis_double<Lattice>(
			lattice, M, rdms, U1, enable_state_optimality_conditions)
  {
    epsilon = this->M_->variable("epsilon");
    const int number_of_moments =
        static_cast<int>(this->lattice_.variable_map_.size());
    upper_box_multiplier_ = this->M_->variable(
        "upper box multipliers", number_of_moments,
        Domain::greaterThan(0.0));
    lower_box_multiplier_ = this->M_->variable(
        "lower box multipliers", number_of_moments,
        Domain::greaterThan(0.0));
    const int nrblocks_x = nrblocks_x_for(this->lattice_);
    const int nrblocks_y = nrblocks_y_for(this->lattice_);
    for (auto sign_symm_sector : this->sectors_)
    {
      Xs_[sign_symm_sector.first] = {};
	  if (this->enable_state_optimality_conditions_)
		state_optimality_Xs_[sign_symm_sector.first] = {};
      for (int i = 0; i < nrblocks_x; i++)
      {
        Xs_[sign_symm_sector.first].push_back({});
		if (this->enable_state_optimality_conditions_)
		  state_optimality_Xs_[sign_symm_sector.first].push_back({});
        for (int j = 0; j < nrblocks_y; j++)
        {
          if (!is_momentum_representative(this->lattice_, i, j))
          {
            // Non-representative conjugate blocks are encoded by the
            // representative block for their {k,-k} orbit.
            Xs_[sign_symm_sector.first][i].push_back(nullptr);
			if (this->enable_state_optimality_conditions_)
			  state_optimality_Xs_[sign_symm_sector.first][i]
				  .push_back(nullptr);
            continue;
          }
          int matrix_dimension = 2 * sign_symm_sector.second.block_shifts[i][j];
          auto X = this->M_->variable("X_" + std::to_string(sign_symm_sector.first) + "_" +
                                           std::to_string(i) + std::to_string(j),
                                       Domain::inPSDCone(matrix_dimension));
          if (maximize_)
            Xs_[sign_symm_sector.first][i].push_back(Expr::neg(X));
          else
            Xs_[sign_symm_sector.first][i].push_back(X);

		  if (this->enable_state_optimality_conditions_)
		  {
			const int state_optimality_matrix_dimension = 2 *
				sign_symm_sector.second.state_optimality_block_shifts[i][j];
			auto state_optimality_X = this->M_->variable(
				"state_optimality_X_" +
					std::to_string(sign_symm_sector.first) + "_" +
					std::to_string(i) + "_" + std::to_string(j),
				Domain::inPSDCone(state_optimality_matrix_dimension));
			if (maximize_)
			  state_optimality_Xs_[sign_symm_sector.first][i].push_back(
				  Expr::neg(state_optimality_X));
			else
			  state_optimality_Xs_[sign_symm_sector.first][i].push_back(
				  state_optimality_X);
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
	Expression::t totalvec = A_vector;
	if (state_optimality_A_vector.get() != nullptr)
	  totalvec = Expr::add(totalvec, state_optimality_A_vector);
	if (Lamba_vector.get() != nullptr)
	  totalvec = Expr::add(totalvec, Lamba_vector);
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
    auto box_contribution = maximize_
        ? Expr::sub(upper_box_multiplier_, lower_box_multiplier_)
        : Expr::sub(lower_box_multiplier_, upper_box_multiplier_);
    totalvec = Expr::add(totalvec, box_contribution);
    if (zero_moment_multiplier_.get() != nullptr)
    {
      const int zero_index = this->lattice_.variable_map_.at("0");
      auto zero_selector = Matrix::sparse(
          static_cast<int>(this->lattice_.variable_map_.size()), 1,
          monty::new_array_ptr(std::vector<int>{zero_index}),
          monty::new_array_ptr(std::vector<int>{0}),
          monty::new_array_ptr(std::vector<double>{1.0}));
      totalvec = Expr::add(
          totalvec, Expr::mul(zero_selector, zero_moment_multiplier_));
    }
    auto updated_expression =
        Expr::add(Expr::add(totalvec, epsilon_vec_flat), this->b_);
    if (this->lattice_.uses_1d_reflection())
    {
      final_constraint_->remove();
      final_constraint_ = this->M_->constraint(
          updated_expression, Domain::equalsTo(0.0));
    }
    else
    {
      final_constraint_->update(updated_expression);
    }
  }

  void fix_constrains()
  {
    if (this->lattice_.variable_map_.find("0") !=
        this->lattice_.variable_map_.end())
    {
      zero_moment_multiplier_ = this->M_->variable("zero moment multiplier");
    }

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
    const int nrblocks_x = nrblocks_x_for(this->lattice_);
    const int nrblocks_y = nrblocks_y_for(this->lattice_);

    for (auto &sign_symm_sector : this->sectors_)
    {
      for (int i = 0; i < nrblocks_x; i++)
      {
        for (int j = 0; j < nrblocks_y; j++)
        {
          if (!is_momentum_representative(this->lattice_, i, j))
            continue;
          int block_size = 2 * sign_symm_sector.second.block_shifts[i][j];
          if (block_size == 0)
            continue;

          int n_vars_block = block_size * block_size;
          auto x_block = Expr::reshape(Xs_[sign_symm_sector.first][i][j], n_vars_block);

          std::vector<int> rows_b, cols_b;
          std::vector<double> vals_b;

          const auto add_block_entries =
              [&](const Matrix::t &mat, int moment_index) {
                auto data = mat->getDataAsArray();

                for (int r = 0; r < block_size; r++)
                  for (int c = 0; c < block_size; c++)
                  {
                    double v = (*data)[r * block_size + c];
                    if (std::abs(v) > 1e-15)
                    {
                      rows_b.push_back(moment_index);
                      cols_b.push_back(r * block_size + c);
                      vals_b.push_back(v);
                    }
                  }
              };

          for (auto &op : this->lattice_.variable_map_)
          {
            if (op.first == "0")
              continue;
            int el = op.second;
            auto &A_block = this->As_[op.first][sign_symm_sector.first][i][j];
            if (A_block.has_elements_)
            {
              // Use only the representative momentum block.  The realified
              // construction already packages the conjugate block into this
              // real PSD variable; adding the mirror block a second time
              // over-constrains the relaxation.
              auto mat = A_block.make_matrix(block_size, block_size);
              add_block_entries(mat, el);
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

	if (this->enable_state_optimality_conditions_)
	{
	  for (auto &sign_symm_sector : this->sectors_)
	  {
		for (int i = 0; i < nrblocks_x; ++i)
		{
		  for (int j = 0; j < nrblocks_y; ++j)
		  {
			if (!is_momentum_representative(this->lattice_, i, j))
			  continue;
			const int block_size = 2 *
				sign_symm_sector.second.state_optimality_block_shifts[i][j];
			if (block_size == 0)
			  continue;

			const int n_vars_block = block_size * block_size;
			auto x_block = Expr::reshape(
				state_optimality_Xs_[sign_symm_sector.first][i][j],
				n_vars_block);
			std::vector<int> rows, cols;
			std::vector<double> values;

			for (const auto &op : this->lattice_.variable_map_)
			{
			  if (op.first == "0")
				continue;
			  auto &coefficient_matrix =
				  this->state_optimality_As_[op.first]
					  [sign_symm_sector.first][i][j];
			  if (!coefficient_matrix.has_elements_)
				continue;

			  auto matrix = coefficient_matrix.make_matrix(
				  block_size, block_size);
			  auto data = matrix->getDataAsArray();
			  for (int row = 0; row < block_size; ++row)
				for (int col = 0; col < block_size; ++col)
				{
				  const double value = (*data)[row * block_size + col];
				  if (std::abs(value) > 1e-15)
				  {
					rows.push_back(op.second);
					cols.push_back(row * block_size + col);
					values.push_back(value);
				  }
				}
			}

			if (values.empty())
			  continue;
			auto coefficient_map = Matrix::sparse(
				n_constraints, n_vars_block,
				monty::new_array_ptr(rows), monty::new_array_ptr(cols),
				monty::new_array_ptr(values));
			auto contribution = Expr::mul(coefficient_map, x_block);
			if (state_optimality_A_vector.get() == nullptr)
			  state_optimality_A_vector = contribution;
			else
			  state_optimality_A_vector =
				  Expr::add(state_optimality_A_vector, contribution);
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
    const int sigma_index = ll++;
    int block_size = (int)std::round(std::sqrt(lambda_expr->getSize()));
    int n_vars_block = block_size * block_size;
    auto l_block = Expr::reshape(lambda_expr, n_vars_block);
  //std::cout<< "matrix size "<<block_size <<std::endl;
    std::vector<int>    rows_b, cols_b;
    std::vector<double> vals_b;
//for(auto& elements: this->sigmas_[key])
auto elements= this->sigmas_[key][sigma_index];

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
    }
}
    int el = this->lattice_.variable_map_.at("1");
    auto e_vec = Matrix::sparse(
        n_constraints, 1, monty::new_array_ptr(std::vector<int>{el}),
        monty::new_array_ptr(std::vector<int>{0}),
        monty::new_array_ptr(std::vector<double>{1.0}));
    auto epsilon_vec = Expr::mul(e_vec, epsilon);
    epsilon_vec_flat = Expr::reshape(epsilon_vec, n_constraints);

	Expression::t totalvec = A_vector;
	if (state_optimality_A_vector.get() != nullptr)
	  totalvec = Expr::add(totalvec, state_optimality_A_vector);
    if (Lamba_vector.get() != nullptr)
	{ totalvec = Expr::add(totalvec, Lamba_vector);}
    if (this->nr_of_linear_constraints > 0)
    {
      std::cout<< "apply linear constarints "<<std::endl;
      LC_vector = Expr::mul(this->Psp, linear_constraints_variable2_);
      totalvec = Expr::add(totalvec, LC_vector);
    }
    if (this->bounding_observable_)
    {
      auto exp_temporary =
          Expr::mul(Expr::add(energy_bouding_variables_[0], energy_bouding_variables_[1]), this->energy_vec_);
      totalvec = Expr::add(totalvec, exp_temporary);
    }
    auto box_contribution = maximize_
        ? Expr::sub(upper_box_multiplier_, lower_box_multiplier_)
        : Expr::sub(lower_box_multiplier_, upper_box_multiplier_);
    totalvec = Expr::add(totalvec, box_contribution);
    if (zero_moment_multiplier_.get() != nullptr)
    {
      const int zero_index = this->lattice_.variable_map_.at("0");
      auto zero_selector = Matrix::sparse(
          n_constraints, 1,
          monty::new_array_ptr(std::vector<int>{zero_index}),
          monty::new_array_ptr(std::vector<int>{0}),
          monty::new_array_ptr(std::vector<double>{1.0}));
      totalvec = Expr::add(
          totalvec, Expr::mul(zero_selector, zero_moment_multiplier_));
    }

    final_constraint_ = this->M_->constraint(
        Expr::add(Expr::add(totalvec, epsilon_vec_flat), this->b_), Domain::equalsTo(0.));

    return;
  }
  Expression::t get_costfunction()
  {
    Expression::t ee = Expr::constTerm(0.);
    ee = Expr::add(ee, Expr::neg(epsilon));

    auto box_constant = Expr::sum(
        Expr::add(upper_box_multiplier_, lower_box_multiplier_));
    if (maximize_)
      ee = Expr::sub(ee, box_constant);
    else
      ee = Expr::add(ee, box_constant);

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
