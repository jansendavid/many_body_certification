#pragma once
#include "fusion.h"
#include "spins.hpp"
#include <unordered_map>
#include <memory>
#include "symmetries.hpp"
using namespace mosek::fusion;
using namespace monty;
using int_pair = std::pair<int, int>;
using TI_map_type = std::map<std::string, std::pair<std::string, std::complex<double>>>;
class LatticeBase
{
public:
	LatticeBase(int Ly, int Lx) : Ly_(Ly), Lx_(Lx) {};
	int Ly_;
	int Lx_;
	TI_map_type TI_map_;
	void generate_TI_map(std::map<std::string, op_vec> &mat_terms, std::vector<op_vec> &operators_, int sign_sector_) {};
	struct G_el
	{
		double prefac_;
		int pos_;
		G_el(double prefac, int pos) : prefac_(prefac), pos_(pos) {};
	};

	struct G_op
	{
		std::complex<double> prefac_;
		std::string op_;
		G_op(std::complex<double> prefac, std::string op) : prefac_(prefac), op_(op) {};
	};

	G_op generate_G_element_sos(op_vec op1, op_vec op2, int j, int i, TI_map_type &TI_map_)
	{
		// Generate all elements of the first row with translation in y direction. j go in y direction

		auto op_dagg_first = dagger_operator(op1);

		auto [fac_dagg, op_dagger] = get_normal_form(op_dagg_first);

		assert(std::abs(fac_dagg.imag()) < 1e-9);
		op_vec new_op_y;
		op_vec new_op;
		if (j > 0)
		{
			new_op_y = translation_y(op2, j, Ly_);
		}
		else
		{

			new_op_y = op2;
		}
		if (i > 0)
		{
			new_op = translation(new_op_y, i, Lx_);
		}
		else
		{
			new_op = new_op_y;
		}
		auto v_x = op_dagger;

		v_x.insert(v_x.end(), new_op.begin(), new_op.end());
		auto [fac, vec] = get_normal_form(v_x);

		auto [ti_key, ti_val] = TI_map_.at(print_op(vec));

		cpx total_fac = fac; // fac*ti_val;
		return G_op(total_fac, ti_key);
	}
};

class SquareLattice : public LatticeBase
{
public:
	std::string permuts_;
	// sign symmetry of the Hamiltonian
	std::string signsym_;

	bool bilayer_;
	bool square_;
	std::vector<int> get_offset_vec()
	{
		if (bilayer_)
		{
			return {2, Ly_, Lx_};
		}
		else
		{
			return {Ly_, Lx_};
		}
	}

	SquareLattice(int Ly, int Lx, bool square, bool bilayer, std::string permuts = "xyz", std::string signsym = "xyz") : LatticeBase(Ly, Lx), bilayer_(bilayer), square_(square), permuts_(permuts), signsym_(signsym)
	{
		// assert(Lx == Ly);
		if (permuts != "xyz" and permuts != "xy" and permuts != "None")
		{
			std::cout << "permutation error" << std::endl;
		}
		if (signsym != "xyz" and signsym != "y" and signsym != "None")
		{
			std::cout << "sign symmetrie error" << std::endl;
		}
	};
	template <typename container>
	std::pair<bool, std::string> check_if_operator_exists(op_vec op, container &mat_terms, bool rdm_check = false)
	{
		// rdm_check check if exist, if not print something
		//  check if the operator is contained in functions
		// returs if_found, and where
		if (op.size() < 1)
		{

			auto it = mat_terms.find(print_op(op));
			if (it != mat_terms.end())
			{
				return {true, print_op(op)};
			}
			else
			{
				return {false, print_op(op)};
			}
		}
		if (signsym_ == "xyz")
		{
			if (is_zero_signsym_xyz(op))
			{

				if (mat_terms.find("0") != mat_terms.end())
				{

					return {true, "0"};
				}
				else
				{

					return {false, "0"};
				}
			}
		}
		else if (signsym_ == "y")
		{
			if (is_zero_signsym_y(op))
			{

				if (mat_terms.find("0") != mat_terms.end())
				{

					return {true, "0"};
				}
				else
				{

					return {false, "0"};
				}
			}
		}
		auto all_ty = generate_all_translations_y(op, Ly_, 1);

		bool found = false;

		for (auto op_ty : all_ty)
		{
			auto all_t = generate_all_translations(op_ty, Lx_);

			for (auto op_t : all_t)
			{
				auto it = mat_terms.find(print_op(op_t));
				if (it != mat_terms.end())
				{

					// TI_map_.insert({print_op(op), {it->first, 1}});
					return {true, it->first};
				}
				std::vector<op_vec> all_p;
				if (permuts_ == "xyz" or permuts_ == "yxz" or permuts_ == "zxy" or permuts_ == "zyx")
				{
					all_p = generate_all_permutations_xyz(op_t);
				}
				else if (permuts_ == "xy")
				{
					all_p = generate_all_permutations_xy(op_t);
				}
				else if (permuts_ == "None")
				{
					all_p.push_back(op_t);
				}
				for (auto op_p : all_p)
				{
					if (!square_)
					{
						auto it = mat_terms.find(print_op(op_p));
						if (it != mat_terms.end())
						{

							return {true, it->first};
						}
					}
					else
					{
						auto all_d8sym = generate_all_d8(op_p, Lx_);

						for (auto d8s : all_d8sym)
						{
							auto it = mat_terms.find(print_op(d8s));
							if (it != mat_terms.end())
							{

								return {true, it->first};
							}

							auto op_mirror = mirror(d8s);

							it = mat_terms.find(print_op(op_mirror));
							if (it != mat_terms.end())
							{

								return {true, it->first};
							}
							if (bilayer_)
							{
								auto op_flip_layer = flip_layer(d8s);
								it = mat_terms.find(print_op(op_flip_layer));
								if (it != mat_terms.end())
								{

									return {true, it->first};
								}
								op_flip_layer = flip_layer(op_mirror);
								it = mat_terms.find(print_op(op_flip_layer));
								if (it != mat_terms.end())
								{

									return {true, it->first};
								}
							}
						}
					}
				}
			}
		}
		if (rdm_check)
		{

			std::cout << print_op(op) << " in rdm did not exist" << std::endl;
			assert(false);
		}
		return {false, print_op(op)};
	}
	void generate_TI_map(std::map<std::string, op_vec> &mat_terms, std::vector<op_vec> &operators_, int sign_sector)
	{
		//     // generate all elemenets with translation symmetrie in x and y direction

		int index = 0;
		for (auto it1 = operators_.begin(); it1 != operators_.end(); ++it1)
		{
			auto [fac, vec] = get_normal_form(*it1);

			auto [found, op_string] = check_if_operator_exists(vec, mat_terms);
			if (found)
			{
				TI_map_.insert({print_op(vec), {op_string, 1}});
			}
			else
			{
				// std::cout << "called" << std::endl;
				mat_terms.insert({print_op(vec), vec});
				TI_map_.insert({print_op(vec), {print_op(vec), 1}});
			}

			// diagonal elements

			for (auto it2 = it1; it2 != operators_.end(); ++it2)
			{

				// if(it2!=it1)
				{
					// std::cout<<"xx"<<std::endl;
					auto op_cp = *it2;

					for (int n = 0; n < Ly_; n++)
					{
						for (int m = 0; m < Lx_; m++)
						{
							op_vec new_op_y;
							op_vec new_op;
							if (n > 0)
							{

								new_op_y = translation_y(op_cp, n, Ly_);
							}
							else
							{
								new_op_y = op_cp;
							}

							if (m > 0)
							{

								new_op = translation(new_op_y, m, Lx_);
							}
							else
							{
								new_op = new_op_y;
							}

							auto op_dagger = dagger_operator(*it1);
							auto [fac, v_x] = get_normal_form(op_dagger);
							assert(std::abs(fac.imag() < 1e-9));

							v_x.insert(v_x.end(), new_op.begin(), new_op.end());

							auto [fac_tot, vec_tot] = get_normal_form(v_x);

							auto [found, op_string] = check_if_operator_exists(vec_tot, mat_terms);
							if (found)
							{
								TI_map_.insert({print_op(vec_tot), {op_string, 1}});
							}
							else
							{
								mat_terms.insert({op_string, vec_tot});
								TI_map_.insert({op_string, {op_string, 1}});
							}
						}
					}
				}
			}
		}

		return;
	}
};
