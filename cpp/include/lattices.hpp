#pragma once
#include "fusion.h"
#include "spins.hpp"
#include <unordered_map>
#include <memory>
#include "symmetries.hpp"
#include "reduced_dms.hpp"
using namespace mosek::fusion;
using namespace monty;
using int_pair = std::pair<int, int>;

struct op_key_hash
{
	std::size_t operator()(const op_key &k) const noexcept
	{
		// 64-bit-ish mix; deterministic, fast, good enough for small vectors.
		std::size_t h = 1469598103934665603ull;
		for (std::uint32_t x : k)
		{
			h ^= static_cast<std::size_t>(x) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
		}
		return h;
	}
};

using TI_map_type = std::unordered_map<op_key, std::pair<op_key, std::complex<double>>, op_key_hash>;

class LatticeBase
{
public:
	LatticeBase(int Ly, int Lx) : Ly_(Ly), Lx_(Lx) {};
	int Ly_;
	int Lx_;
	TI_map_type TI_map_;
	std::map<std::string, int> variable_map_;
	std::unordered_map<op_key, std::pair<std::complex<double>, op_vec>, op_key_hash> nf_cache;
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
		G_op(std::complex<double> prefac, std::string op) : prefac_(prefac), op_(std::move(op)) {};
	};
	auto get_nf_cached(const op_vec &op)
    {
        op_key key = key_dir_pos(op);

        auto it = nf_cache.find(key);
        if (it != nf_cache.end())
            return it->second;

        auto res = get_normal_form(op);
        nf_cache.emplace(key, res);
        return res;
    }
	G_op generate_G_element_sos_double(op_vec op1, op_vec op2, int j, int i)
	{
		// Generate all elements of the first row with translation in y direction. j go in y direction

		auto op_dagg_first = dagger_operator(op1);

		// auto [fac_dagg, op_dagger] = get_normal_form(op_dagg_first);

		// assert(std::abs(fac_dagg.imag()) < 1e-9);
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
		auto v_x = op_dagg_first;

		v_x.insert(v_x.end(), new_op.begin(), new_op.end());
		// auto [fac, vec] = get_normal_form(v_x);
		// std::cout << print_op(v_x) << std::endl;
		auto [fac, nf] = get_nf_cached(v_x);
		if(nf.size()%2==0)
		{
			
		}

		auto [ti_key, ti_val] = TI_map_.at(key_dir_pos(nf));
		// std::cout << "end" << std::endl;
		//  assert(fac == ti_val);std::cout<<"start"<<std::endl;

		cpx total_fac = fac * ti_val;
		return G_op(total_fac, op_key_label(ti_key));
	}
	G_op generate_G_element_sos(op_vec op1, op_vec op2, int j, int i)
	{
		// Generate all elements of the first row with translation in y direction. j go in y direction

		auto op_dagg_first = dagger_operator(op1);

		// auto [fac_dagg, op_dagger] = get_normal_form(op_dagg_first);

		// assert(std::abs(fac_dagg.imag()) < 1e-9);
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
		auto v_x = op_dagg_first;

		v_x.insert(v_x.end(), new_op.begin(), new_op.end());

		auto [fac, nf] = get_normal_form(v_x);

		auto [ti_key, ti_val] = TI_map_.at(key_dir_pos(nf));
		
		cpx total_fac = fac * ti_val;
		return G_op(total_fac, op_key_label(ti_key));
	}
};
template<typename Basis>
class SquareLattice : public LatticeBase
{
public:
	std::string permuts_;
	// sign symmetry of the Hamiltonian
	std::string signsym_;
	Basis states_;
	// vector in which all elements found while looking for translation invariance will be added and flushed (reset ) once an element is added
	std::vector<op_vec> flush_vector;

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

	SquareLattice(Basis& states, int Ly, int Lx, bool square, bool bilayer, std::string permuts = "xyz", std::string signsym = "xyz") : LatticeBase(Ly, Lx), states_(states),bilayer_(bilayer), square_(square), permuts_(permuts), signsym_(signsym)
	{
		// assert(Lx == Ly);
		if (permuts != "xyz" and permuts != "yxz" and permuts != "zxy" and permuts != "xy" and permuts != "None")
		{
			std::cout << "permutation error" << std::endl;
		}
		if (signsym != "xyz" and signsym != "y" and signsym != "None")
		{
			std::cout << "sign symmetrie error" << std::endl;
		}
	};
	
	void flush(op_vec op_in)
	{

		auto [fac_key, nf_key] = get_nf_cached(op_in);
		const op_key nf_key_k = key_dir_pos(nf_key);

		for (auto &op : flush_vector)
		{
			// std::cout << print_op(op) << std::endl;
			auto [fac, nf] = get_nf_cached(op);

			TI_map_.insert({key_dir_pos(nf),
							{nf_key_k, std::conj(fac) * fac_key}});
		}

		flush_vector.clear();
	};
	bool check_operator_translation(op_vec op)
	{
		if(op.size()<1)
		{
			auto it = TI_map_.find(key_dir_pos(op));
		
			if (it != TI_map_.end())
			{

				TI_map_.insert({key_dir_pos(op),
								{it->second.first, 1. * it->second.second}});
				// flush(it->second.first, fac);

				return true;
			}
			else{
				return false;
			}
		
		}
		bool found = false;

		auto all_t = all_translations(op, Lx_, Ly_);
		auto [fac_op, nf_op] = get_nf_cached(op);

		for (const auto &op_t : all_t)
		{

			auto [fac, nf] = get_nf_cached(op_t);
			auto it = TI_map_.find(key_dir_pos(nf));
			flush_vector.push_back(op_t);
			if (it != TI_map_.end())
			{

				TI_map_.insert({key_dir_pos(nf_op),
								{it->second.first, std::conj(fac_op) * fac * it->second.second}});
				// flush(it->second.first, fac);

				return true;
			}
			else
			{

				found = check_additional_symmetries(op, op_t);
			}
		}

		return found;
	}
	bool check_permutation_symm(op_vec op_org, op_vec op)
	{
	
		std::set<op_vec> all_p;
		auto [fac_org, nf_org] = get_nf_cached(op_org);
		auto [fac, nf] = get_nf_cached(op);
		if (permuts_ == "xyz" or permuts_ == "yxz" or permuts_ == "zxy" or permuts_ == "zyx")
		{
			
			all_p = generate_all_permutations_xyz(nf);
		}
		else if (permuts_ == "xy")
		{
			all_p = generate_all_permutations_xy(nf);
		}
		else if (permuts_ == "None")
		{
			all_p.insert(nf);
		}

		for (auto op_p : all_p)
		{
			//auto [fac, nf] = get_normal_form(op_p);
			auto it = TI_map_.find(key_dir_pos(op_p));

			if (it != TI_map_.end())
			{

				TI_map_.insert({key_dir_pos(nf_org),
								{it->second.first, std::conj(fac_org) * fac * it->second.second}});
				// flush(it->second.first, fac);
				return true;
			}
			else
			{
				//flush_vector.push_back(op_p);
			}
		}
		return false;
	}
	bool check_additional_symmetries(op_vec op_org, op_vec op)
	{
		auto [fac_org, nf_org] = get_nf_cached(op_org);
		bool found = false;
		if (square_)
		{
			auto dsvec = generate_all_d8(op, Lx_);
			for (auto &d8s : dsvec)
			{
				auto [fac, nf] = get_nf_cached(d8s);
				auto it = TI_map_.find(key_dir_pos(nf));

				auto mirrored_ds8 = mirror(d8s);

				auto [fac_mir, nf_mir] = get_nf_cached(mirrored_ds8);

				if (it != TI_map_.end())
				{

					TI_map_.insert({key_dir_pos(nf_org),
									{it->second.first, std::conj(fac_org) * fac * it->second.second}});

					return true;
				}
				else
				{
					flush_vector.push_back(d8s);
				}
				auto it_mirrored = TI_map_.find(key_dir_pos(nf_mir));
				if (it_mirrored != TI_map_.end())
				{

					TI_map_.insert({key_dir_pos(nf_org),
									{it_mirrored->second.first, std::conj(fac_org) * fac_mir * it_mirrored->second.second}});

					return true;
				}
				else
				{
					flush_vector.push_back(mirrored_ds8);
				}

				found = check_permutation_symm(op_org, d8s);

				if (found)
				{
					return found;
				}

				found = check_permutation_symm(op_org, mirrored_ds8);

				if (found)
				{
					return found;
				}

				if (bilayer_)
				{
					auto op_flip_layer = flip_layer((d8s));
					auto [fac_flip, nf_flip] = get_nf_cached(op_flip_layer);
					auto it_flip = TI_map_.find(key_dir_pos(nf_flip));
					if (it_flip != TI_map_.end())
					{

						TI_map_.insert({key_dir_pos(nf_org),
										{it->second.first, std::conj(fac_org) * fac_flip * it->second.second}});

						return true;
					}
					else
					{
						flush_vector.push_back(op_flip_layer);
					}
					found = check_permutation_symm(op_org, op_flip_layer);

					if (found)
					{
						return found;
					}

					auto op_flip_layer_mirr = flip_layer(mirrored_ds8);
					auto [fac_flip_mirr, nf_flip_mirr] = get_nf_cached(op_flip_layer_mirr);
					auto it_flip_mirr = TI_map_.find(key_dir_pos(nf_flip_mirr));
					if (it_flip_mirr != TI_map_.end())
					{

						TI_map_.insert({key_dir_pos(nf_org),
										{it->second.first, std::conj(fac_org) * fac_flip_mirr * it->second.second}});

						return true;
					}
					else
					{
						flush_vector.push_back(op_flip_layer_mirr);
					}
					found = check_permutation_symm(op_org, op_flip_layer_mirr);

					if (found)
					{
						return found;
					}
				}
			}
		}
		else
		{
			found = check_permutation_symm(op_org, op);
			if (found)
			{
				return found;
			}
			found = check_permutation_symm(op_org, mirror(op));
			if (found)
			{
				return found;
			}
		}

		return false;
	}

	std::pair<op_key, std::complex<double>> get_key(op_vec spin_op)
	{

		auto [fac, nf] = get_nf_cached(spin_op);

		op_key key = key_dir_pos(nf);

		if (signsym_ == "xyz")
		{
			if (is_zero_signsym_xyz(nf))
			{
				key = op_key_zero();
			}
		}
		else if (signsym_ == "xy")
		{
			if (is_zero_signsym_xy(nf))
			{
				key = op_key_zero();
			}
		}
		else if (signsym_ == "y")
		{
			if (is_zero_signsym_y(nf))
			{
				key = op_key_zero();
			}
		}
		else
		{
		}
		return std::pair<op_key, std::complex<double>>(key, fac);
	}
	void clear_caches() {
        nf_cache.clear();

    }
	bool see_if_state_exists(op_vec spin_op)
	{
		flush_vector.clear();
		bool found = false;
		found = check_operator_translation(spin_op);
		return found;
	}
	void generate_TI_map()
	{

		for (auto sector : states_)
		{
			std::cout << "sector " << sector.first << std::endl;
			auto operators = sector.second;
			for (auto it1 = operators.begin(); it1 != operators.end(); ++it1)
			{
				auto op = *it1;
			
				for (auto it2 = it1; it2 != operators.end(); ++it2)
				{
					//std::cout << " op 2: " << print_op(*it1) << std::endl;
					auto op_dagg_first = dagger_operator(op);
					std::vector<op_vec> all_t;
					if(it2->size()>0)
					{
						 all_t = all_translations(*it2, Lx_, Ly_);
					}
					else{
						all_t.push_back(*it2);
					}
					for (auto &op_right : all_t)
					{
						flush_vector.clear();
						auto v_x = op_dagg_first;

						v_x.insert(v_x.end(), op_right.begin(), op_right.end());
						auto [key, fac] = get_key(v_x);

						auto [fac_, nf] = get_nf_cached(v_x);
						bool found = false;
						if (is_zero_key(key))
						{
						}
						else
						{
							found = check_operator_translation(v_x);
						}
						if (found == false)
						{

							TI_map_.insert({key_dir_pos(nf),
											{key, 1}});

							flush(v_x);
						}
					}
				}
			}
			clear_caches();
		}
	

		return;
	}
	void operator_run(std::vector<op_vec>& operators_1, std::vector<op_vec>& operators_2)
	{
		for (auto it1 = operators_1.begin(); it1 != operators_1.end(); ++it1)
			{
				auto op=*it1;
			for (auto it2 = operators_2.begin(); it2 != operators_2.end(); ++it2)
				{
					auto op_dagg_first = dagger_operator(op);
					std::vector<op_vec> all_t;
					if(it2->size()>0)
					{
						 all_t = all_translations(*it2, Lx_, Ly_);
					}
					else{
						all_t.push_back(*it2);
					}
							for (auto &op_right : all_t)
					{
						flush_vector.clear();
						auto v_x = op_dagg_first;

						v_x.insert(v_x.end(), op_right.begin(), op_right.end());
						auto [key, fac] = get_key(v_x);

						auto [fac_, nf] = get_nf_cached(v_x);
						if(nf.size()%2!=0)
						{key = op_key_zero();}
						bool found = false;
						if (is_zero_key(key))
						{
						}
						else
						{
							found = check_operator_translation(v_x);
						}
						if (found == false)
						{

							TI_map_.insert({key_dir_pos(nf),
											{key, 1}});

							flush(v_x);
						}
					}
				}
	
				clear_caches();
			}

		return;
	}
	void generate_TI_map_double()
	{

		for (auto &sector : states_)
		{
			operator_run(sector.second.at(0), sector.second.at(0));
			operator_run(sector.second.at(1), sector.second.at(1));
			operator_run(sector.second.at(0), sector.second.at(1));
			operator_run(sector.second.at(1), sector.second.at(0));

		}
	
	

		return;
	}
	void make_map()
	{

		std::set<std::string> unique_values;

		for (const auto &[k, v] : TI_map_)
		{

			unique_values.insert(op_key_label(v.first));
		}

		int i = 0;
		for (auto a : unique_values)
		{

			variable_map_.insert({a, i});
			i += 1;
		}

		return;
	}

	std::map<std::string, Matrix::t> generate_rdms_primal_cp(rdm_operator sites, std::vector<int> offset)
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

				auto [key, fac] = get_key(state);
				if (is_zero_key(key))
				{
				}
				else
				{
					auto [fac, nf] = get_nf_cached(state);
					bool found = see_if_state_exists(nf);

					if (!found)
					{
						std::cout << "adding rdm operator" << std::endl;
						TI_map_.insert({key_dir_pos(nf),
										{key, 1}});
					}
					auto it = TI_map_.find(key_dir_pos(nf));

					const std::string rep_label = op_key_label(it->second.first);

					assert(std::abs((fac * it->second.second).imag()) < 1e-9);
					mat = mat * (fac * it->second.second).real() / std::pow(2, degree);

					if (rdms_eigen_.find(rep_label) != rdms_eigen_.end())
					{

						rdms_eigen_[rep_label] += mat;
					}
					else
					{

						rdms_eigen_.insert({rep_label, mat});
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
	std::pair<std::complex<double>, op_vec> get_form_of_TI_map(const op_vec &op)
	{

		auto [coeff_, nf] = get_nf_cached(op);

		return std::pair<std::complex<double>, op_vec>(coeff_, nf);
	}
};
