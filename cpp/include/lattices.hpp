#pragma once
#include <algorithm>
#include "fusion.h"
#include "spins.hpp"
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include "operator_operations.hpp"
#include "symmetries.hpp"
#include "reduced_dms.hpp"
#include <cstdlib>
using namespace mosek::fusion;
using namespace monty;
using int_pair = std::pair<int, int>;

template <typename T>
struct is_std_pair : std::false_type
{
};

template <typename First, typename Second>
struct is_std_pair<std::pair<First, Second>> : std::true_type
{
};

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
	LatticeBase(int Lx, int Ly) : Ly_(Ly), Lx_(Lx) {};
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
		// j is x translation (pos_x), i is y translation (pos_y)

		auto op_dagg_first = dagger_operator(op1);

		op_vec new_op_x;
		op_vec new_op;
		auto number_of_indices=op1[0].get_site().size();
		if (j > 0)
		{
			new_op_x = translation(op2, j, Lx_, number_of_indices-2);
		}
		else
		{
			new_op_x = op2;
		}
		if (i > 0)
		{
			new_op = translation(new_op_x, i, Ly_, number_of_indices-1);
		}
		else
		{
			new_op = new_op_x;
		}
		auto v_x = op_dagg_first;

		v_x.insert(v_x.end(), new_op.begin(), new_op.end());

		auto [fac, nf] = get_nf_cached(v_x);

		auto [ti_key, ti_val] = TI_map_.at(key_dir_pos(nf));
		
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
		auto number_of_indices=op1[0].get_site().size();
		if (j > 0)
		{
			new_op_y = translation(op2, j, Lx_,number_of_indices-2);
		}
		else
		{

			new_op_y = op2;
		}
		if (i > 0)
		{
			new_op = translation(new_op_y, i, Ly_,number_of_indices-1);
		}
		else
		{
			new_op = new_op_y;
		}
		auto v_x = op_dagg_first;

		v_x.insert(v_x.end(), new_op.begin(), new_op.end());

		auto [fac, nf] = get_nf_cached(v_x);

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
	Basis state_optimality_states_;
	int state_optimality_basis_degree_{3};
	// vector in which all elements found while looking for translation invariance will be added and flushed (reset ) once an element is added
	std::vector<op_vec> flush_vector;

	bool bilayer_;
	bool square_;
	std::set<op_vec> extra_states_;
	std::optional<SumOfOperators> state_optimality_hamiltonian_;
	std::vector<int> offset_vec_;

	struct StateOptimalityEntryKey
	{
		op_key v;
		op_key w;
		int dx{0};
		int dy{0};

		bool operator==(const StateOptimalityEntryKey &other) const
		{
			return v == other.v && w == other.w && dx == other.dx &&
				   dy == other.dy;
		}
	};

	struct StateOptimalityEntryKeyHash
	{
		std::size_t operator()(const StateOptimalityEntryKey &key) const noexcept
		{
			op_key_hash hash_op;
			std::size_t hash = hash_op(key.v);
			hash ^= hash_op(key.w) + 0x9e3779b97f4a7c15ull + (hash << 6) +
					(hash >> 2);
			hash ^= static_cast<std::size_t>(key.dx) +
					0x9e3779b97f4a7c15ull + (hash << 6) + (hash >> 2);
			hash ^= static_cast<std::size_t>(key.dy) +
					0x9e3779b97f4a7c15ull + (hash << 6) + (hash >> 2);
			return hash;
		}
	};

	std::unordered_map<op_key, std::vector<std::string>, op_key_hash>
		state_optimality_anticommuting_terms_cache_;
	std::unordered_map<StateOptimalityEntryKey, SumOfOperators,
					   StateOptimalityEntryKeyHash>
		state_optimality_entry_cache_;

	// Fast TI construction.  The old check/flush routines below are retained as
	// reference helpers, but generate_TI_map[_double]() now uses these compact
	// key caches.  A stored relation means raw Pauli string = coefficient *
	// representative Pauli string.
	struct FastNormalizedImage
	{
		op_key raw_key;
		op_key symmetry_key;
		std::complex<double> phase{1.0, 0.0};
		bool zero{false};
	};
	std::unordered_map<op_key, std::pair<op_key, std::complex<double>>,
					   op_key_hash> fast_orbit_cache_;
	std::set<op_key> fast_forced_zero_representatives_;

	bool is_zero_by_sign_symmetry(const op_vec &normal_form) const
	{
		if (signsym_ == "xyz")
			return is_zero_signsym_xyz(normal_form);
		if (signsym_ == "xy")
			return is_zero_signsym_xy(normal_form);
		if (signsym_ == "y")
			return is_zero_signsym_y(normal_form);
		return false;
	}

	void insert_fast_relation(const op_key &raw_key,
						  const op_key &representative,
						  std::complex<double> coefficient,
						  const char *source)
	{
		constexpr double tolerance = 1e-10;
		auto [it, inserted] = TI_map_.insert(
			{raw_key, {representative, coefficient}});
		if (inserted)
			return;
		if (it->second.first == representative &&
			std::abs(it->second.second - coefficient) <= tolerance)
			return;

		if (is_zero_key(representative))
		{
			if (!is_zero_key(it->second.first))
				fast_forced_zero_representatives_.insert(it->second.first);
			it->second = {op_key_zero(), {1.0, 0.0}};
			return;
		}
		if (is_zero_key(it->second.first) ||
			fast_forced_zero_representatives_.count(representative) != 0)
		{
			it->second = {op_key_zero(), {1.0, 0.0}};
			return;
		}
		if (it->second.first == representative)
		{
			// Symmetry gives R=cR with c!=1, hence this moment must vanish.
			fast_forced_zero_representatives_.insert(representative);
			it->second = {op_key_zero(), {1.0, 0.0}};
			return;
		}

		std::cerr << "TI_map contradiction from " << source << '\n'
				  << op_key_label(raw_key) << " -> "
				  << it->second.second << " * "
				  << op_key_label(it->second.first) << " versus "
				  << coefficient << " * "
				  << op_key_label(representative) << '\n';
		throw std::logic_error("incompatible fast TI-map relation");
	}

	void apply_fast_forced_zeros()
	{
		if (fast_forced_zero_representatives_.empty())
			return;
		for (auto &[raw, relation] : TI_map_)
		{
			(void)raw;
			if (fast_forced_zero_representatives_.count(relation.first) != 0)
				relation = {op_key_zero(), {1.0, 0.0}};
		}
	}

	FastNormalizedImage normalize_fast_image(const op_vec &image,
										 bool force_odd_zero = false)
	{
		auto [phase, normal_form] = get_nf_cached(image);
		FastNormalizedImage result;
		result.raw_key = key_dir_pos(normal_form);
		result.phase = phase;
		result.zero = is_zero_by_sign_symmetry(normal_form) ||
			(force_odd_zero && ((normal_form.size() & 1u) != 0));
		result.symmetry_key = result.zero ? op_key_zero() : result.raw_key;
		return result;
	}

	template <typename Fn>
	void for_each_fast_symmetry_image(const op_vec &op, Fn &&fn)
	{
		const auto emit_with_permutations = [&](const op_vec &spatial_image) {
			fn(spatial_image, std::complex<double>{1.0, 0.0});
			auto [phase, normal_form] = get_nf_cached(spatial_image);
			std::set<op_vec> permutations;
			if (permuts_ == "xyz" || permuts_ == "yxz" ||
				permuts_ == "zxy" || permuts_ == "zyx")
				permutations = generate_all_permutations_xyz(normal_form);
			else if (permuts_ == "xy")
				permutations = generate_all_permutations_xy(normal_form);
			else
				permutations.insert(normal_form);
			for (const auto &permuted : permutations)
			{
				// Permutations are generated from the normalized spatial
				// image, so retain the phase extracted from that image.
				fn(permuted, phase);
			}
		};

		const auto emit_spatial = [&](const op_vec &translated) {
			const bool one_dimensional = (Lx_ == 1) != (Ly_ == 1);
			if (square_ && one_dimensional)
			{
				for (const auto &reflected :
					 generate_all_1d_reflections(translated, Lx_, Ly_))
					emit_with_permutations(reflected);
			}
			else if (square_)
			{
				for (const auto &d8 : generate_all_d8(translated, Lx_))
				{
					emit_with_permutations(d8);
					emit_with_permutations(mirror(d8));
					if (bilayer_)
					{
						emit_with_permutations(flip_layer(d8));
						emit_with_permutations(flip_layer(mirror(d8)));
					}
				}
			}
			else
			{
				emit_with_permutations(translated);
				emit_with_permutations(mirror(translated));
			}
		};

		if (op.empty())
		{
			emit_spatial(op);
			return;
		}
		for (int dx = 0; dx < Lx_; ++dx)
			for (int dy = 0; dy < Ly_; ++dy)
			{
				op_vec translated = op;
				for (auto &factor : translated)
				{
					auto site = factor.get_site();
					site[site.size() - 2] =
						(site[site.size() - 2] + dx) % Lx_;
					site[site.size() - 1] =
						(site[site.size() - 1] + dy) % Ly_;
					factor.set_site(site);
				}
				emit_spatial(translated);
			}
	}

	void process_fast_orbit_candidate(const op_vec &op, const char *source,
									 bool force_odd_zero = false)
	{
		const auto original = normalize_fast_image(op, force_odd_zero);
		const auto orbit_cached = fast_orbit_cache_.find(original.raw_key);
		if (orbit_cached != fast_orbit_cache_.end())
		{
			insert_fast_relation(original.raw_key, orbit_cached->second.first,
							 orbit_cached->second.second, source);
			return;
		}
		const auto ti_cached = TI_map_.find(original.raw_key);
		if (ti_cached != TI_map_.end())
		{
			fast_orbit_cache_[original.raw_key] = ti_cached->second;
			return;
		}
		if (original.zero)
		{
			insert_fast_relation(original.raw_key, op_key_zero(), {1.0, 0.0},
							 source);
			fast_orbit_cache_[original.raw_key] =
				{op_key_zero(), {1.0, 0.0}};
			return;
		}

		// Preserve the legacy representative convention: the first-discovered
		// (original) normal form names a new orbit.  Scan lazily for an already
		// known image before inserting the full new orbit; most later products
		// therefore avoid a complete set of TI-map writes.
		std::unordered_map<op_key, FastNormalizedImage, op_key_hash> pending;
		bool found_known_image = false;
		bool orbit_zero = false;
		for_each_fast_symmetry_image(op, [&](const op_vec &image,
											 std::complex<double> inherited_phase) {
			if (found_known_image || orbit_zero)
				return;
			auto normalized = normalize_fast_image(image, force_odd_zero);
			normalized.phase *= inherited_phase;

			if (normalized.raw_key == original.raw_key)
			{
				if (std::abs(normalized.phase - original.phase) > 1e-10)
					orbit_zero = true;
				return;
			}
			if (normalized.zero)
			{
				orbit_zero = true;
				return;
			}

			const auto cached_image = fast_orbit_cache_.find(normalized.raw_key);
			const auto mapped_image = TI_map_.find(normalized.raw_key);
			if (cached_image != fast_orbit_cache_.end() ||
				mapped_image != TI_map_.end())
			{
				const auto &known = cached_image != fast_orbit_cache_.end()
					? cached_image->second
					: mapped_image->second;
				if (is_zero_key(known.first) ||
					fast_forced_zero_representatives_.count(known.first) != 0)
				{
					insert_fast_relation(original.raw_key, op_key_zero(),
								 {1.0, 0.0}, source);
					fast_orbit_cache_[original.raw_key] =
						{op_key_zero(), {1.0, 0.0}};
				}
				else
				{
					const auto coefficient = std::conj(original.phase) *
						normalized.phase * known.second;
					insert_fast_relation(original.raw_key, known.first,
								 coefficient, source);
					fast_orbit_cache_[original.raw_key] =
						{known.first, coefficient};
				}
				found_known_image = true;
				return;
			}

			auto [it, inserted] = pending.emplace(normalized.raw_key, normalized);
			if (!inserted &&
				std::abs(it->second.phase - normalized.phase) > 1e-10)
				orbit_zero = true;
		});

		if (found_known_image)
			return;
		if (orbit_zero)
		{
			fast_forced_zero_representatives_.insert(original.raw_key);
			insert_fast_relation(original.raw_key, op_key_zero(), {1.0, 0.0},
							 source);
			fast_orbit_cache_[original.raw_key] =
				{op_key_zero(), {1.0, 0.0}};
			return;
		}

		insert_fast_relation(original.raw_key, original.raw_key, {1.0, 0.0},
						 source);
		fast_orbit_cache_[original.raw_key] =
			{original.raw_key, {1.0, 0.0}};
		for (const auto &[raw, image] : pending)
		{
			const std::complex<double> coefficient =
				std::conj(image.phase) * original.phase;
			insert_fast_relation(raw, original.raw_key, coefficient, source);
			const auto inserted = TI_map_.find(raw);
			if (inserted != TI_map_.end())
			{
				fast_orbit_cache_[raw] = inserted->second;
			}
		}
	}

	void append_translated_fast(op_vec &destination, const op_vec &source,
							 int dx, int dy) const
	{
		for (const auto &factor : source)
		{
			auto translated = factor;
			auto site = translated.get_site();
			site[site.size() - 2] = (site[site.size() - 2] + dx) % Lx_;
			site[site.size() - 1] = (site[site.size() - 1] + dy) % Ly_;
			translated.set_site(site);
			destination.push_back(std::move(translated));
		}
	}

	static std::vector<int> extract_offset_vec(const Basis &states)
	{
		std::vector<int> offset_vec;
		auto visit = [&offset_vec](auto &&self, const auto &node) -> void
		{
			using Node = std::decay_t<decltype(node)>;
			if constexpr (std::is_same_v<Node, spin_op>)
			{
				if (offset_vec.empty())
				{
					offset_vec = node.offset_;
				}
				else if (offset_vec != node.offset_)
				{
					throw std::invalid_argument(
						"SquareLattice states contain inconsistent offset vectors");
				}
			}
			else if constexpr (is_std_pair<Node>::value)
			{
				self(self, node.second);
			}
			else
			{
				for (const auto &child : node)
				{
					self(self, child);
				}
			}
		};

		visit(visit, states);
		if (offset_vec.empty())
		{
			throw std::logic_error(
				"SquareLattice cannot infer an offset vector from empty states");
		}

		return offset_vec;
	}

	static Basis filter_state_optimality_states(const Basis &states,
										   int maximum_degree)
	{
		Basis filtered = states;
		if constexpr (std::is_same_v<Basis, basis_structure_with_sub>)
		{
			for (auto &[sector, subsectors] : filtered)
			{
				(void)sector;
				for (auto &[subsector, operators] : subsectors)
				{
					(void)subsector;
					operators.erase(
						std::remove_if(
							operators.begin(), operators.end(),
							[maximum_degree](const op_vec &op) {
								return static_cast<int>(op.size()) > maximum_degree;
							}),
						operators.end());
				}
			}
		}
		else if constexpr (std::is_same_v<Basis, basis_structure>)
		{
			for (auto &[sector, operators] : filtered)
			{
				(void)sector;
				operators.erase(
					std::remove_if(
						operators.begin(), operators.end(),
						[maximum_degree](const op_vec &op) {
							return static_cast<int>(op.size()) > maximum_degree;
						}),
					operators.end());
			}
		}
		return filtered;
	}

	std::vector<int> get_offset_vec() const
	{
		return offset_vec_;
	}

	bool uses_1d_reflection() const
	{
		return square_ && ((Lx_ == 1) != (Ly_ == 1));
	}

	SquareLattice(Basis &states, int Lx, int Ly, bool square, bool bilayer,
			  std::string permuts = "xyz", std::string signsym = "xyz",
			  std::set<op_vec> extra_states = {},
			  std::optional<SumOfOperators> state_optimality_hamiltonian = std::nullopt,
			  int state_optimality_basis_degree = 3)
		: LatticeBase(Lx, Ly), permuts_(permuts), signsym_(signsym),
		  states_(states),
		  state_optimality_states_(filter_state_optimality_states(
			  states, state_optimality_basis_degree)),
		  state_optimality_basis_degree_(state_optimality_basis_degree),
		  bilayer_(bilayer), square_(square),
		  extra_states_(extra_states),
		  state_optimality_hamiltonian_(std::move(state_optimality_hamiltonian)),
		  offset_vec_(extract_offset_vec(states))
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
		std::cout << "square "<<square_ << std::endl;

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
		flush_vector.clear();
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
				if (found)
				{
					return true;
				}
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
		const bool one_dimensional = (Lx_ == 1) != (Ly_ == 1);
		if (square_ && one_dimensional)
		{
			auto reflections = generate_all_1d_reflections(op, Lx_, Ly_);
			for (auto &reflected : reflections)
			{
				auto [fac, nf] = get_nf_cached(reflected);
				auto it = TI_map_.find(key_dir_pos(nf));
				if (it != TI_map_.end())
				{
					TI_map_.insert({key_dir_pos(nf_org),
								{it->second.first, std::conj(fac_org) * fac * it->second.second}});
					return true;
				}

				flush_vector.push_back(reflected);
				found = check_permutation_symm(op_org, reflected);
				if (found)
				{
					return true;
				}
			}
		}
		else if (square_)
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
						const auto mapped = it_flip->second;
						TI_map_.insert({key_dir_pos(nf_org),
										{mapped.first, std::conj(fac_org) * fac_flip * mapped.second}});

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
						const auto mapped = it_flip_mirr->second;
						TI_map_.insert({key_dir_pos(nf_org),
										{mapped.first, std::conj(fac_org) * fac_flip_mirr * mapped.second}});

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
	const std::vector<std::string> &state_optimality_anticommuting_terms(
		const op_vec &op)
	{
		const auto key = key_dir_pos(op);
		const auto cached = state_optimality_anticommuting_terms_cache_.find(key);
		if (cached != state_optimality_anticommuting_terms_cache_.end())
			return cached->second;

		std::vector<std::string> labels;
		if (state_optimality_hamiltonian_)
		{
			for (const auto &[label, term] :
				 state_optimality_hamiltonian_->get_terms())
				if (pauli_strings_anticommute(op, term.get_op()))
					labels.push_back(label);
		}
		return state_optimality_anticommuting_terms_cache_
			.emplace(key, std::move(labels)).first->second;
	}

	const SumOfOperators &get_state_optimality_entry(
		const op_vec &v, const op_vec &w, int dx, int dy)
	{
		if (!state_optimality_hamiltonian_)
			throw std::logic_error(
				"state-optimality entry requested without a Hamiltonian");

		StateOptimalityEntryKey key{key_dir_pos(v), key_dir_pos(w), dx, dy};
		const auto cached = state_optimality_entry_cache_.find(key);
		if (cached != state_optimality_entry_cache_.end())
			return cached->second;

		op_vec translated_w;
		translated_w.reserve(w.size());
		append_translated_fast(translated_w, w, dx, dy);
		const auto translated_w_dagger = dagger_operator(translated_w);
		const auto &v_terms = state_optimality_anticommuting_terms(v);
		const auto &w_terms =
			state_optimality_anticommuting_terms(translated_w_dagger);
		auto entry = build_state_optimality_entry_from_anticommuting_terms(
			v, translated_w_dagger, *state_optimality_hamiltonian_, v_terms,
			w_terms);
		return state_optimality_entry_cache_
			.emplace(std::move(key), std::move(entry)).first->second;
	}
	bool see_if_state_exists(op_vec spin_op)
	{
		flush_vector.clear();
		bool found = false;
		found = check_operator_translation(spin_op);
		return found;
	}
	void generate_TI_map(bool enable_state_optimality_conditions = false)
	{
		fast_orbit_cache_.clear();
		fast_forced_zero_representatives_.clear();
		TI_map_.reserve(std::max<std::size_t>(TI_map_.size(), 1024));
		insert_fast_relation(op_key_identity(), op_key_identity(), 1.0,
						 "identity base relation");
		insert_fast_relation(op_key_zero(), op_key_zero(), 1.0,
						 "zero base relation");

		for (const auto &sector : states_)
		{
			std::cout << "sector " << sector.first << std::endl;
			const auto &operators = sector.second;
			std::size_t maximum_size = 0;
			for (const auto &op : operators)
				maximum_size = std::max(maximum_size, op.size());
			const std::size_t candidates = operators.size() +
				(operators.size() * (operators.size() + 1) / 2) *
				static_cast<std::size_t>(Lx_) * static_cast<std::size_t>(Ly_);
			TI_map_.reserve(TI_map_.size() + candidates * 8);
			fast_orbit_cache_.reserve(fast_orbit_cache_.size() + candidates * 8);

			op_vec product;
			product.reserve(2 * maximum_size);
			for (auto left = operators.begin(); left != operators.end(); ++left)
			{
				process_fast_orbit_candidate(*left, "single basis operator");
				const auto left_dagger = dagger_operator(*left);
				for (auto right = left; right != operators.end(); ++right)
					for (int dx = 0; dx < Lx_; ++dx)
						for (int dy = 0; dy < Ly_; ++dy)
						{
							product.clear();
							product.insert(product.end(), left_dagger.begin(),
										   left_dagger.end());
							append_translated_fast(product, *right, dx, dy);
							process_fast_orbit_candidate(product,
												 "moment-matrix product");
						}
			}
			clear_caches();
		}
		if (enable_state_optimality_conditions)
		{
			state_optimality_anticommuting_terms_cache_.clear();
			state_optimality_entry_cache_.clear();
			for (auto &sector : state_optimality_states_)
				generate_TI_map_state_optimality_conditions_double(
					sector.second, sector.second);
		}
		for (const auto &state : extra_states_)
			process_fast_orbit_candidate(state, "extra observable");
		apply_fast_forced_zeros();
		return;
	}
	void operator_run(std::vector<op_vec>& operators_1, std::vector<op_vec>& operators_2)
	{
		std::size_t maximum_size = 0;
		for (const auto &op : operators_1)
			maximum_size = std::max(maximum_size, op.size());
		for (const auto &op : operators_2)
			maximum_size = std::max(maximum_size, op.size());
		op_vec product;
		product.reserve(2 * maximum_size);
		for (const auto &left : operators_1)
		{
			const auto left_dagger = dagger_operator(left);
			for (const auto &right : operators_2)
				for (int dx = 0; dx < Lx_; ++dx)
					for (int dy = 0; dy < Ly_; ++dy)
					{
						product.clear();
						product.insert(product.end(), left_dagger.begin(),
									   left_dagger.end());
						append_translated_fast(product, right, dx, dy);
						process_fast_orbit_candidate(product,
											 "double moment product", true);
					}
			clear_caches();
		}
		return;
	}
	void generate_TI_map_state_optimality_conditions_double(
		const std::vector<op_vec> &operators_1,
		const std::vector<op_vec> &operators_2)
	{
		if (!state_optimality_hamiltonian_)
			return;

		for (const auto &v : operators_1)
		{
			for (const auto &w : operators_2)
			{
				for (int dx = 0; dx < Lx_; ++dx)
				{
					for (int dy = 0; dy < Ly_; ++dy)
					{
						const auto &reduced_entry =
							get_state_optimality_entry(v, w, dx, dy);
						for (const auto &[label, term] : reduced_entry.get_terms())
						{
							(void)label;
							process_fast_orbit_candidate(
								term.get_op(),
								"reduced state-optimality entry", true);
						}
					}
				}
			}
			clear_caches();
		}
		return;
	}
	void generate_TI_map_double(bool enable_state_optimality_conditions = false)
	{
		fast_orbit_cache_.clear();
		fast_forced_zero_representatives_.clear();
		state_optimality_anticommuting_terms_cache_.clear();
		state_optimality_entry_cache_.clear();
		insert_fast_relation(op_key_identity(), op_key_identity(), 1.0,
						 "identity base relation");
		insert_fast_relation(op_key_zero(), op_key_zero(), 1.0,
						 "zero base relation");
		for (auto &sector : states_)
		{
			operator_run(sector.second.at(0), sector.second.at(0));
			operator_run(sector.second.at(1), sector.second.at(1));
			operator_run(sector.second.at(0), sector.second.at(1));
			operator_run(sector.second.at(1), sector.second.at(0));

		}
		if (enable_state_optimality_conditions)
		{
			for (auto &sector : state_optimality_states_)
			{
				generate_TI_map_state_optimality_conditions_double(
					sector.second.at(0), sector.second.at(0));
				generate_TI_map_state_optimality_conditions_double(
					sector.second.at(1), sector.second.at(1));
				generate_TI_map_state_optimality_conditions_double(
					sector.second.at(0), sector.second.at(1));
				generate_TI_map_state_optimality_conditions_double(
					sector.second.at(1), sector.second.at(0));
			}
		}
		for (const auto &state : extra_states_)
			process_fast_orbit_candidate(state, "extra double observable", true);
		apply_fast_forced_zeros();
		return;
	}
	void make_map()
	{
		std::unordered_set<op_key, op_key_hash> seen;
		std::vector<std::string> labels;
		labels.reserve(TI_map_.size());
		for (const auto &[raw, relation] : TI_map_)
		{
			(void)raw;
			if (seen.insert(relation.first).second)
				labels.push_back(op_key_label(relation.first));
		}
		std::sort(labels.begin(), labels.end());
		for (std::size_t index = 0; index < labels.size(); ++index)
			variable_map_.insert({labels[index], static_cast<int>(index)});
		return;
	}

	std::vector<std::map<std::string, Matrix::t>> generate_rdms_primal_cp(rdm_operator sites, std::vector<int> offset)
	{

		std::vector<std::map<std::string, Matrix::t>> sigmas_temp_;
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
			auto [key, fac_i] = get_key(state);
			auto [fac, nf] = get_nf_cached(state);
			if (print_op(nf) == "1")
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

			
				if (is_zero_key(key))
				{
				}
				else
				{
				 
					bool found = see_if_state_exists(nf);

					if (!found)
					{
						//std::cout << "adding rdm operator" << std::endl;
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
		std::map<std::string, Matrix::t> sigmas_temp_el;
		for (auto eigen_matrix : rdms_eigen_)
		{
			auto Alpha = get_sparse_from_eigen(eigen_matrix.second);

			sigmas_temp_el.insert({eigen_matrix.first, Alpha});
		}
		sigmas_temp_.push_back(sigmas_temp_el);
		return sigmas_temp_;
	}
	std::vector<std::vector<int>> generate_binary_vectors(int L, int N) {
		std::vector<std::vector<int>> result;
	
		// initial vector: N ones, L-N zeros
		std::vector<int> v(L, 0);
		for (int i = 0; i < N; ++i)
			v[i] = 1;
	
		// generate all permutations
		do {
			result.push_back(v);
		} while (std::prev_permutation(v.begin(), v.end()));
	
		return result;
	}
	struct U1rdm_element
	{
		std::vector<std::string> operators;
		std::pair<int,int> indices;
		int dim{0};
	};
	std::vector<U1rdm_element> make_terms(
		std::vector<std::vector<int>> vecs)
	{
		std::map<std::pair<int,int>, std::string> stringmap;
		stringmap[{0,0}]="n";
		stringmap[{1,1}]="(1-n)";
		stringmap[{0,1}]="cdag";
		stringmap[{1,0}]="c";
		std::map<std::string, mat_type> rdms_eigen_;
		auto M=mat_type::Zero(vecs.size(), vecs.size());
		std::vector<U1rdm_element> results;
		
		for(int i=0; i<vecs.size(); i++)
		{
			for(int j=0; j<vecs.size(); j++)
			{
				// std::string s="";
				U1rdm_element result;
				result.dim=vecs.size();
				for(int l=0; l<vecs[0].size(); l++)
				{
					result.operators.push_back(stringmap[{vecs[i][l],vecs[j][l] }]);
						// s+=stringmap[{vecs[i][l],vecs[j][l] }]+"_"+std::to_string(l);
				}
				// int site1=1;
				// int site2=2;
				// int layers=0;
				// std::vector<int> offset={2,3,3};
				// std::cout<< "endtry ("<<i<<","<<j<< ") = "<<s<<std::endl;
				result.indices={i,j};
	results.push_back(result);
			}
	
		}
	
	 return results;}
	 std::vector<std::pair<std::complex<double>, op_vec>>get_res(std::string s, std::vector<int> indices, std::vector<int> offset)
{
    std::vector<std::pair<std::complex<double>, op_vec>> res;
    if(s=="n")
    {
        res.push_back({1./2, {}});
        res.push_back({1./2, {spin_op("z", indices, offset)}});
    }
    else if(s=="c")
    {
        res.push_back({1./2, {spin_op("x", indices, offset)}});
        res.push_back({std::complex<double>(0,-1.)*1./2., {spin_op("y", indices, offset)}});
    }
    else if(s=="cdag")
    {
        res.push_back({1./2, {spin_op("x", indices, offset)}});
        res.push_back({std::complex<double>(0,1.)*1./2., {spin_op("y", indices, offset)}});
    }
    else if(s=="(1-n)")
    {
        res.push_back({1./2, {}});
        res.push_back({-1./2, {spin_op("z", indices, offset)}});
    }
    else{
        std::cout<< "error: "<<s<<std::endl;
    }
    return res;
}
template<typename T>
std::map<std::string, Matrix::t>  get_temp_sig(T& rdms_eigen_){
	//  returns one matrix with U(1) symm
	std::map<std::string, Matrix::t> sigmas_temp_;

		for (auto eigen_matrix : rdms_eigen_)
		{
			auto Alpha = get_sparse_from_eigen(eigen_matrix.second);

			sigmas_temp_.insert({eigen_matrix.first, Alpha});
			//std::cout<< eigen_matrix.first <<" inserting "<< eigen_matrix.second.rows()<<std::endl;
		}
		
		//sigmas_temp_.push_back(sigmas_temp_el);
	
	return sigmas_temp_;
}
	std::vector<std::map<std::string, Matrix::t>> generate_rdms_primal_U1(rdm_operator sites, std::vector<int> offset)
	{
		
		std::vector<std::map<std::string, Matrix::t>> sigmas_temp_;
		std::vector<std::map<std::string, mat_type>> rdms_eigen_;
		std::map<std::string, mat_type> sigma_map;
		std::map<std::pair<int,int>, std::string> stringmap;
		// double check convention
		stringmap[{0,0}]="n";
		stringmap[{1,1}]="(1-n)";
		stringmap[{0,1}]="cdag";
		stringmap[{1,0}]="c";
		std::vector<std::vector<U1rdm_element>>  matrices;
		for(int i=0; i<=sites.size(); i++)
		{
			auto res=generate_binary_vectors(sites.size(), i);
			auto obj=make_terms(res);
			// for(auto a:obj)
			// {std::cout<<a.indices.first<< ";"<<a.indices.second<<std::endl;}
			matrices.push_back(obj);
			// std::cout<<"start"<<std::endl;
			// auto res= generate_binary_vectors(L,  i);
	
			// std::cout<< binom(L, i) << " and "<<res.size()<<std::endl;
			// make_terms(res);

	//		std::cout<<std::endl;

		}
	//	std::cout<<"mats "<<matrices.size()<<std::endl;
		for(auto& vect_of_op: matrices)
		{
			
			//std::cout<< "mat runcs "<<std::endl;
			std::map<std::string, mat_type> rdms_eigen_temp_;
			for(auto op: vect_of_op)
			{
			
		
				int n=0;
				//std::pair<std::complex<double>, op_vec> initial_pair={{1.0, {}}};
				std::vector<std::pair<std::complex<double>, op_vec>> total_ops={{1.0, {}}};
				for(auto op_string: op.operators)
				{
					std::vector<std::pair<std::complex<double>, op_vec>> next;
					//std::cout<< "string op "<<op_string<<std::endl;
					std::vector<std::pair<std::complex<double>, op_vec>> conv=get_res(op_string, sites.at(n), offset);
					//std::cout<<"conv "<<conv.size()<<std::endl;
					
					for(auto& obj: total_ops)
					{
								for(auto& final_op:conv )
				{
					auto new_state=obj;
					new_state.first*=final_op.first;
					new_state.second.insert(new_state.second.end(),final_op.second.begin(),final_op.second.end());
				//	append(final_op.second);
					next.push_back(std::move(new_state));
					
					}

					}
					total_ops=std::move(next);
					n++;
				}

				
			 
				for(auto final_op:total_ops)
				{
					
			 		auto [key, fac] = get_key(final_op.second);
				if (is_zero_key(key))
				{
				}
				else{

				auto [fac, nf] = get_nf_cached(final_op.second);
					bool found = see_if_state_exists(nf);
					mat_type mat=mat_type::Zero(op.dim, op.dim);
					mat(op.indices.first, op.indices.second)=1.;
					if (!found)
					{
						//std::cout << "adding rdm operator" << std::endl;
						TI_map_.insert({key_dir_pos(nf),
										{key, 1}});
					}
					auto it = TI_map_.find(key_dir_pos(nf));

					const std::string rep_label = op_key_label(it->second.first);
					if (rdms_eigen_temp_.find(rep_label) != rdms_eigen_temp_.end())
			{
				
				
				rdms_eigen_temp_[rep_label] += mat*final_op.first;
			}
			else
			{

				rdms_eigen_temp_.insert({rep_label, mat*final_op.first});
			}
			 	}
			// 	rdms_eigen_.push_back(rdms_eigen_temp_);
				} // end of iteratying over total obs
				
			 }
			 auto sig_output=  get_temp_sig(rdms_eigen_temp_);
				sigmas_temp_.push_back(sig_output);
			}
			
			
		
		
		//exit(1);
		// auto M=mat_type::Zero(vecs.size(), vecs.size());
		// sites.at(i)

		std::cout<< "sig size "<<sigmas_temp_.size()<<std::endl;
		return sigmas_temp_;
	}
	std::pair<std::complex<double>, op_vec> get_form_of_TI_map(const op_vec &op)
	{

		auto [coeff_, nf] = get_nf_cached(op);

		return std::pair<std::complex<double>, op_vec>(coeff_, nf);
	}
};
spin_op unpack_word(uint32_t word,
	const std::vector<int>& offset)
{
uint32_t dir_code = word & 0x3;
uint32_t pos      = word >> 2;

std::string dir;
switch (dir_code)
{
case 1: dir = "x"; break;
case 2: dir = "y"; break;
case 3: dir = "z"; break;
default:
throw std::runtime_error("invalid dir code");
}

std::vector<int> shifts(offset.size());
shifts[0] = 1;

for (size_t i = 1; i < offset.size(); ++i)
shifts[i] = shifts[i - 1] * offset[i - 1];

std::vector<int> site(offset.size());

for (int i = static_cast<int>(offset.size()) - 1; i >= 0; --i)
{
site[i] = pos / shifts[i];
pos %= shifts[i];
}

return spin_op(dir, site, offset);
}

op_vec unpack_key(const op_key& key,
	const std::vector<int>& offset)
{
op_vec result;
result.reserve(key.size());

for (auto word : key)
result.push_back(unpack_word(word, offset));

return result;
}
