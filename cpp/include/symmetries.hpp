#pragma once
#include "fusion.h"
#include "spins.hpp"
#include <unordered_map>
#include <memory>
#include <set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
using namespace mosek::fusion;
using namespace monty;
std::shared_ptr<ndarray<int, 1>> nint(const std::vector<int> &X) { return new_array_ptr<int>(X); }
std::shared_ptr<ndarray<double, 1>> ndou(const std::vector<double> &X) { return new_array_ptr<double>(X); }
template <typename T>
std::vector<std::vector<T>>
all_translations(const std::vector<T> &config, int Lx, int Ly)
{
  std::vector<std::vector<T>> result;
  int size_of_vec = config[0].get_site().size();
  result.reserve(Lx * Ly);

  for (int dx = 0; dx < Lx; ++dx)
  {
    for (int dy = 0; dy < Ly; ++dy)
    {

      std::vector<T> translated = config;

      for (auto &s : translated)
      {
        auto old_site = s.get_site();

        old_site[size_of_vec - 1] = (old_site[size_of_vec - 1] + dx) % Lx;
        old_site[size_of_vec - 2] = (old_site[size_of_vec - 2] + dy) % Ly;
        s.set_site(old_site);
      }

      result.push_back(std::move(translated));
    }
  }

  return result;
}
struct rdm_operator
{

  std::vector<std::vector<int>> op_;
  rdm_operator(std::vector<std::vector<int>> op) : op_(op) {};
  rdm_operator() {};
  unsigned int size()
  {
    return op_.size();
  }
  std::vector<int> at(int i)
  {
    return op_[i];
  }
  bool operator<(const rdm_operator &rhs) const
  {
    return op_ < rhs.op_;
  }
};
struct rdms_struct
{
  // struct to manage the reduced density matrices
  // store a vector with operator and each opertor has vector of its sites
  std::vector<rdm_operator> rdms;
  rdms_struct(std::vector<rdm_operator> rdms) : rdms(rdms) {};
  rdms_struct() {};
  void add_operator(rdm_operator new_operator_indices)
  {
    // appends operators, e.g. (1,2)(0,3) (2,2)
    rdms.push_back(new_operator_indices);
  }
  unsigned int size()
  {
    return rdms.size();
  }
};
rdms_struct get_rdms(int Lx, int dim)
{
  rdms_struct data;

  rdm_operator newstate({{0, 0}, {0, 1}});
  data.add_operator(newstate);
  if (dim >= 4)
  {
    rdm_operator newstate({{0, 0}, {0, 1}, {0, 2}});
    data.add_operator(newstate);
    // rdm_operator newstate_1({{0, 0}, {0, 1}, {0, 2}, {0, 3}});
    // data.add_operator(newstate_1);
    // rdm_operator newstate_2({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}});
    // data.add_operator(newstate_2);
    // rdm_operator newstate_3({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}});
    // data.add_operator(newstate_3);
    // rdm_operator newstate_4({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}});
    // data.add_operator(newstate_4);
    // rdm_operator newstate_5({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 2}});
    // data.add_operator(newstate_5);

    // rdm_operator newstate_55({{0, 0}, {2, 1}, {0, 2}, {2, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 2}});
    // data.add_operator(newstate_55);
    //  rdm_operator newstate_6({{0,0}, {0,1}, {0,2},{0,3},{1,1},{1,2},{1,3},{2,2},{3,3}});
    // data.add_operator(newstate_6);
    //  rdm_operator newstate_7({{0,0}, {0,1}, {0,2},{0,3},{1,1},{1,2},{1,3},{2,2},{3,3}});
    // data.add_operator(newstate_7);
  }
  if (dim >= 6)
  {
    rdm_operator newstate({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}});
    data.add_operator(newstate);
    rdm_operator newstate_1({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}});
    data.add_operator(newstate_1);
  }

  return data;
}

rdms_struct get_rdms_1d(int Lx, int dim)
{
  rdms_struct data;

  rdm_operator newstate({{0, 0}, {0, 1}});
  data.add_operator(newstate);
  if (dim >= 4)
  {
    rdm_operator newstate({{0, 0}, {0, 1}, {0, 2}});
    data.add_operator(newstate);
    // rdm_operator newstate_1({{0, 0}, {0, 1}, {0, 2}, {0, 3}});
    // data.add_operator(newstate_1);
    // rdm_operator newstate_2({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}});
    // data.add_operator(newstate_2);
    // rdm_operator newstate_3({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}});
    // data.add_operator(newstate_3);
    // rdm_operator newstate_4({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}});
    // data.add_operator(newstate_4);
    // rdm_operator newstate_5({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 2}});
    // data.add_operator(newstate_5);

    // rdm_operator newstate_55({{0, 0}, {2, 1}, {0, 2}, {2, 3}, {1, 1}, {1, 2}, {1, 3}, {2, 2}});
    // data.add_operator(newstate_55);
    //  rdm_operator newstate_6({{0,0}, {0,1}, {0,2},{0,3},{1,1},{1,2},{1,3},{2,2},{3,3}});
    // data.add_operator(newstate_6);
    //  rdm_operator newstate_7({{0,0}, {0,1}, {0,2},{0,3},{1,1},{1,2},{1,3},{2,2},{3,3}});
    // data.add_operator(newstate_7);
  }
  // if (dim >= 6)
  // {
  //   rdm_operator newstate({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}});
  //   data.add_operator(newstate);
  //   rdm_operator newstate_1({{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}});
  //   data.add_operator(newstate_1);
  // }

  return data;
}

rdms_struct translation_invariant_rdms_2nd(int Lx, int Ly)
{
  rdms_struct data;
  std::string s = "x";
  //   std::set<std::string> list_of_operators;
  //   for(int i=0; i<Lx; i++)
  //   {
  //     for(int j=0; j<Lx; j++)
  //     {
  //       if(i!=0 and j!=0)
  //       {
  //         op_vec v0={spin_op(s, {0,0}, {Ly,Lx}),spin_op(s, {i,j},{Ly,Lx})};
  //         auto [fac, vec] =get_normal_form(v0);
  //         std::vector<std::vector<int>> states;
  //         for(auto O: vec)
  //         {
  //           states.push_back(O.sites_);
  //         }
  //         rdm_operator newstate(states);

  // data.add_operator(newstate);
  //       }
  //     }
  //   }

  return data;
}

rdms_struct translation_invariant_rdms_3d(int Lx, int Ly)
{
  rdms_struct data;
  std::vector<std::vector<int>> combs;
  for (int i = 0; i < Lx; i++)
  {
    for (int j = 0; j < Lx; j++)
    {
      if (i != 0 and j != 0)
      {
        combs.push_back({i, j});
      }
    }
  }
  for (auto p : combs)
  {
    for (auto o : combs)
    {
      for (auto l : combs)
      {
        if (p != o and p != l and o != l)
        {
          rdm_operator newstate({{0, 0}, {o}, {p}, {l}});
          data.add_operator(newstate);
        }
      }
    }
  }

  return data;
}
rdms_struct translation_invariant_rdms_4th(int Lx, int Ly)
{
  rdms_struct data;
  std::vector<std::vector<int>> combs;
  for (int i = 0; i < Lx; i++)
  {
    for (int j = 0; j < Lx; j++)
    {
      if (i != 0 and j != 0)
      {
        combs.push_back({i, j});
      }
    }
  }
  for (auto p : combs)
  {
    for (auto o : combs)
    {
      for (auto l : combs)
      {
        for (auto m : combs)
        {
          std::set<std::vector<int>> set_ = {p, o, l, m};

          if (set_.size() == 4)
          {
            rdm_operator newstate({{0, 0}, {o}, {p}, {l}, {m}});
            data.add_operator(newstate);
          }
        }
      }
    }
  }

  return data;
}
template <typename T>
T apply_group_trafo(const T &op, const Eigen::Matrix2i &mat, int L)
{
  T vec;
  vec.reserve(op.size()); // preallocate

  for (const auto &o : op)
  {
    // Access last two coordinates directly
    auto site = o.get_site();
    int n = o.get_site().size();
    int x = site[n - 1];
    int y = site[n - 2];

    // Wrap around L/2
    if (x > L / 2)
      x -= L;
    if (y > L / 2)
      y -= L;

    // Apply linear transformation manually (avoid Eigen small vector overhead)
    int x_new = mat(0, 0) * x + mat(0, 1) * y;
    int y_new = mat(1, 0) * x + mat(1, 1) * y;

    // Map back to [0, L)
    x_new = (x_new + L) % L;
    y_new = (y_new + L) % L;

    // Create new sites vector (reuse existing storage)
    auto sites_new = o.get_site();
    sites_new[n - 1] = x_new;
    sites_new[n - 2] = y_new;

    vec.emplace_back(o.get_dir(), std::move(sites_new), o.offset_);
  }

  return vec;
}

// op_vec apply_group_trafo(op_vec op, Eigen::Matrix2i mat, int L)
// {

//   op_vec vec;
//   for (int i = 0; i < op.size(); i++)
//   {

//     int x_cor = op[i].site_[op[i].site_.size() - 1];
//     int y_cor = op[i].site_[op[i].site_.size() - 2];
//     if (x_cor > L / 2)
//     {
//       x_cor = -L + x_cor;
//     }
//     if (y_cor > L / 2)
//     {
//       y_cor = -L + y_cor;
//     }

//     Eigen::Vector2i location(x_cor, y_cor);

//     auto new_location = mat * location;
//     int x_cor_new = (L + new_location(0)) % L;
//     int y_cor_new = (L + new_location(1)) % L;

//     auto sites_new = op[i].site_;
//     sites_new[op[i].site_.size() - 1] = x_cor_new;
//     sites_new[op[i].site_.size() - 2] = y_cor_new;
//     vec.push_back(spin_op(op[i].dir_, sites_new, op[i].offset_));
//   }
//   return vec;
// }
template <class F>
inline void for_each_d8(const op_vec &op, int L, F &&f)
{
  const Eigen::Matrix2i a = (Eigen::Matrix2i() << 0, -1,
                             1, 0)
                                .finished();

  const Eigen::Matrix2i x = (Eigen::Matrix2i() << 1, 0,
                             0, -1)
                                .finished();

  Eigen::Matrix2i a_pow = Eigen::Matrix2i::Identity();

  for (int k = 0; k < 4; ++k)
  {
    f(apply_group_trafo(op, a_pow, L));
    a_pow = a * a_pow;
  }

  a_pow = Eigen::Matrix2i::Identity();
  for (int k = 0; k < 4; ++k)
  {
    f(apply_group_trafo(op, a_pow * x, L));
    a_pow = a * a_pow;
  }
}
template <typename T>
std::vector<T> generate_all_d8(const T &op, int L)
{
  std::vector<T> all_d8;
  all_d8.reserve(8);

  // Generators of D8
  const Eigen::Matrix2i a = (Eigen::Matrix2i() << 0, -1,
                             1, 0)
                                .finished(); // rotation by pi/2

  const Eigen::Matrix2i x = (Eigen::Matrix2i() << 1, 0,
                             0, -1)
                                .finished(); // reflection

  // auto [fac, vec] = get_normal_form(op);

  // Precompute powers of a
  const Eigen::Matrix2i a2 = a * a;
  const Eigen::Matrix2i a3 = a2 * a;

  // Identity
  all_d8.push_back(op);

  // Rotations
  all_d8.push_back(apply_group_trafo(op, a, L));  // a
  all_d8.push_back(apply_group_trafo(op, a2, L)); // a^2
  all_d8.push_back(apply_group_trafo(op, a3, L)); // a^3

  // Reflections
  all_d8.push_back(apply_group_trafo(op, x, L));      // x
  all_d8.push_back(apply_group_trafo(op, a * x, L));  // ax
  all_d8.push_back(apply_group_trafo(op, a2 * x, L)); // a^2 x
  all_d8.push_back(apply_group_trafo(op, a3 * x, L)); // a^3 x

  return all_d8;
}

// std::vector<op_vec> generate_all_d8(op_vec op, int L)
// {
//   std::vector<op_vec> all_d8;
//   Eigen::Matrix2i a(2, 2);
//   a(0, 0) = 0;
//   a(0, 1) = -1;
//   a(1, 0) = 1;
//   a(1, 1) = 0;

//   Eigen::Matrix2i x;
//   x(0, 0) = 1;
//   x(0, 1) = 0;
//   x(1, 0) = 0;
//   x(1, 1) = -1;

//   auto [fac, vec] = get_normal_form(op);
//   // unit element
//   all_d8.push_back(vec);
//   // a (rotation pi/2)
//   all_d8.push_back(apply_group_trafo(vec, a, L));
//   // a^2
//   all_d8.push_back(apply_group_trafo(vec, a * a, L));

//   // a^3
//   all_d8.push_back(apply_group_trafo(vec, a * a * a, L));

//   // // x
//   all_d8.push_back(apply_group_trafo(vec, x, L));

//   // ax
//   all_d8.push_back(apply_group_trafo(vec, a * x, L));

//   // a^2x

//   all_d8.push_back(apply_group_trafo(vec, a * a * x, L));

//   // a^3x
//   all_d8.push_back(apply_group_trafo(vec, a * a * a * x, L));

//   return all_d8;
// }
template <typename operator_class>
operator_class translation_y(operator_class op, int j, int L)
{
  operator_class vec;
  for (int i = 0; i < op.size(); i++)
  {
    vec.push_back(op[i].get_translated_y(j, L));
  }
  return vec;
}
template <typename operator_class>
std::vector<operator_class> generate_all_translations_y(operator_class op, int L, int inc = 1)
{

  std::vector<operator_class> all_T;
  if (print_op(op) == "1")
  {
    all_T.push_back(op);
    return all_T;
  }

  all_T.push_back(op);
  for (int i = 1; i < L; i += inc)
  {
    auto new_op = translation_y(op, i, L);

    all_T.push_back(new_op);
  }
  return all_T;
}
template <typename operator_class>
operator_class translation(operator_class op, int j, int L)
{
  operator_class vec;
  for (int i = 0; i < op.size(); i++)
  {
    vec.push_back(op[i].get_translated(j, L));
  }
  return vec;
}
template <typename operator_class>
std::vector<operator_class> generate_all_translations(operator_class op, int L)
{
  std::vector<operator_class> all_T;

  if (print_op(op) == "1")
  {
    all_T.push_back(op);
    return all_T;
  }

  all_T.push_back(op);

  for (int i = 1; i < L; i++)
  {
    auto new_op = translation(op, i, L);
    //    do we need auto
    all_T.push_back(new_op);
  }
  return all_T;
}
std::pair<int, int> get_sec(op_vec op)
{
  // computes the sector of a given vector of operators
  // NOTE: all are of the form vec{S_i}vec{S_j} in the Hamiltonian is needed. For e.g., TFI, small modifications must be made
  int sxy = 1;
  int syz = 1;
  for (auto a : op)
  {
    if (a.get_dir() == "x")
    {
      sxy *= -1;
    }
    if (a.get_dir() == "y")
    {
      syz *= -1;
      sxy *= -1;
    }
    if (a.get_dir() == "z")
    {
      syz *= -1;
    }
  }
  return std::pair<int, int>(sxy, syz);
}

void add_state(basis_structure &states, op_vec op, std::map<std::pair<int, int>, int> map_sec)
{
  // adds a state to a basis
  auto [fac, nf] = get_normal_form(op);
  auto sign = get_sec(nf);
  if (nf.size() > 0)
  {
    states.at(map_sec.at(sign)).push_back(nf);
  }

  return;
}

void add_state_with_symmetries(basis_structure &states, op_vec op, std::map<std::pair<int, int>, int> map_sec, int Ly, int Lx)
{
  // adds a state to a basis

  auto [fac, nf] = get_normal_form(op);
  bool print = false;

  auto sign = get_sec(nf);
  bool found = false;
  if (nf.size() > 0)
  {
    auto all_t = generate_all_translations(nf, Lx);
    bool found = false;

    for (auto op_t : all_t)
    {

      auto all_ty = generate_all_translations_y(op_t, Ly, 1);

      for (auto op_ty : all_ty)
      {
        if (print)
        {
          std::cout << print_op(op_ty) << std::endl;
        }
        auto it = std::find(states.at(map_sec.at(sign)).begin(), states.at(map_sec.at(sign)).end(), op_ty);
        if (it != states.at(map_sec.at(sign)).end())
        {
          // std::cout<<print_op(nf)<< " was "<< print_op(*it)<<std::endl;
          found = true;
          break;
        }
      }
    }

    if (!found)
    {
      states.at(map_sec.at(sign)).push_back(nf);
    }
  }

  return;
}

template <class T>
int getIndex(std::vector<T> v, T K)
{
  auto it = find(v.begin(), v.end(), K);

  // If element was found
  if (it != v.end())
  {

    // calculating the index
    // of K
    int index = it - v.begin();
    return index;
  }
  else
  {
    // If the element is not
    // present in the vector
    return -1;
  }
}
struct mom_ref
{
  mom_ref(Variable::t var, int i1, op_vec vec) : var_(var), i1_(i1), vec_(vec) {}
  Variable::t var_;
  int i1_;
  op_vec vec_;
  mom_ref() = default;
};
template <typename operator_vector>
bool is_zero_signsym_xyz(operator_vector op)
{

  std::vector<std::string> dirs = {"x", "y", "z"};
  for (auto dir_ : dirs)
  {
    int fac = 1;
    for (auto a : op)
    {

      if (a.get_dir() == dir_)
      {
        fac *= -1;
      }
    }
    if (fac < 0)
    {

      return true;
    }
  }

  return false;
}
template <typename operator_vector>
bool is_zero_signsym_xy(operator_vector op)
{

  std::vector<std::string> dirs = {"x", "y"};
  for (auto dir_ : dirs)
  {
    int fac = 1;
    for (auto a : op)
    {

      if (a.get_dir() == dir_)
      {
        fac *= -1;
      }
    }
    if (fac < 0)
    {

      return true;
    }
  }

  return false;
}
template <typename operator_vector>
bool is_zero_signsym_y(operator_vector op)
{

  std::vector<std::string> dirs = {"y"};
  for (auto dir_ : dirs)
  {
    int fac = 1;
    for (auto a : op)
    {

      if (a.get_dir() == dir_)
      {
        fac *= -1;
      }
    }
    if (fac < 0)
    {

      return true;
    }
  }

  return false;
}

void display(char a[], int n)
{
  for (int i = 0; i < n; i++)
  {
    std::cout << a[i] << " ";
  }
  std::cout << std::endl;
}
template <typename T>
T mirror(T op)
{
  T vec;
  for (int i = 0; i < op.size(); i++)
  {
    vec.push_back(op[i].get_mirror());
  }
  // auto [fac, nf] = get_normal_form(vec);

  // assert(fac.imag() < 1e-9);
  return vec;
}
template <typename T>
T flip_layer(T op)
{
  T vec;
  for (int i = 0; i < op.size(); i++)
  {
    vec.push_back(op[i].get_flipped_layer());
  }
  // auto [fac, nf] = get_normal_form(vec);

  // assert(fac.imag() < 1e-9);
  return vec;
}
void findPermutations(char a[], int n)
{

  // Sort the given array
  std::sort(a, a + n);

  // Find all possible permutations
  std::cout << "Possible permutations are:\n";
  do
  {
    display(a, n);
  } while (std::next_permutation(a, a + n));
}
std::set<op_vec> generate_all_permutations_xy(op_vec op)
{

  std::set<op_vec> all_P;
  all_P.insert(op);

  std::vector<std::map<std::string, std::string>> permutations(1);
  permutations[0].insert({"x", "y"});
  permutations[0].insert({"y", "x"});
  permutations[0].insert({"z", "z"});

  for (auto &a : permutations)
  {
    auto new_op = op;
    std::for_each(new_op.begin(), new_op.end(), [a](spin_op n)
                  {
       auto old_d=n.get_dir();
       n.set_dir(a.at(old_d)); });
    // auto [fac, vec] = get_normal_form(new_op);
    // assert(fac.imag() < 1e-9);
    all_P.insert(new_op);
  }

  return all_P;
}
std::string apply_permutation(
    const std::string &s,
    const std::map<std::string, std::string> &perm)
{
  std::string out = s;
  for (char &c : out)
  {
    auto it = perm.find(std::string(1, c));
    if (it != perm.end())
      c = it->second[0];
  }
  return out;
}
// std::vector<op_vec> generate_all_permutations_xyz_2(op_vec op)
// {

//   std::vector<op_vec> all_P;
//   // auto [fac, vec] = get_normal_form(op);
//   // assert(fac.imag() < 1e-9);

//   std::vector<std::map<std::string, std::string>> permutations(4);

//   //  trivial permutation
//   all_P.push_back(op);
//   // permutations[0].insert({"x", "x"});
//   // permutations[0].insert({"y", "z"});
//   // permutations[0].insert({"z", "y"});

//   permutations[0].insert({"x", "y"});
//   permutations[0].insert({"y", "x"});
//   permutations[0].insert({"z", "z"});

//   permutations[1].insert({"x", "y"});
//   permutations[1].insert({"y", "z"});
//   permutations[1].insert({"z", "x"});

//   permutations[2].insert({"x", "z"});
//   permutations[2].insert({"y", "x"});
//   permutations[2].insert({"z", "y"});

//   permutations[3].insert({"x", "z"});
//   permutations[3].insert({"y", "y"});
//   permutations[3].insert({"z", "x"});

//   for (auto &a : permutations)
//   {

//     auto new_op = op;
//     std::for_each(new_op.begin(), new_op.end(), [a](spin_op &n)
//                   {
//        auto old_d=n.get_dir();
//        n.set_dir(a.at(old_d)); });
//     // auto [fac, vec] = get_normal_form(new_op);
//     // assert(fac.imag() < 1e-9);
//     all_P.push_back(new_op);
//   }

//   return all_P;
// }

std::set<op_vec> generate_all_permutations_xyz(op_vec op)
{

  std::set<op_vec> all_P;
  // auto [fac, vec] = get_normal_form(op);
  // assert(fac.imag() < 1e-9);

  std::vector<std::map<std::string, std::string>> permutations(4);

  //  trivial permutation
  all_P.insert(op);
  // permutations[0].insert({"x", "x"});
  // permutations[0].insert({"y", "z"});
  // permutations[0].insert({"z", "y"});

  permutations[0].insert({"x", "y"});
  permutations[0].insert({"y", "x"});
  permutations[0].insert({"z", "z"});

  permutations[1].insert({"x", "y"});
  permutations[1].insert({"y", "z"});
  permutations[1].insert({"z", "x"});

  permutations[2].insert({"x", "z"});
  permutations[2].insert({"y", "x"});
  permutations[2].insert({"z", "y"});

  permutations[3].insert({"x", "z"});
  permutations[3].insert({"y", "y"});
  permutations[3].insert({"z", "x"});

  for (auto &a : permutations)
  {

    auto new_op = op;
    std::for_each(new_op.begin(), new_op.end(), [a](spin_op n)
                  {
       auto old_d=n.get_dir();
       n.set_dir(a.at(old_d)); });
    // auto [fac, vec] = get_normal_form(new_op);
    // assert(fac.imag() < 1e-9);
    all_P.insert(new_op);
  }

  return all_P;
}

std::map<std::pair<int, int>, int> get_sector_map()
{
  std::map<std::pair<int, int>, int> map_sec;
  map_sec.insert({std::pair<int, int>(1, 1), 0});
  map_sec.insert({std::pair<int, int>(1, -1), 1});
  map_sec.insert({std::pair<int, int>(-1, 1), 2});
  map_sec.insert({std::pair<int, int>(-1, -1), 3});
  return map_sec;
}
basis_structure get_states()
{
  basis_structure states;
  std::vector<op_vec> v_block_0;
  std::vector<op_vec> v_block_1;
  std::vector<op_vec> v_block_2;
  std::vector<op_vec> v_block_3;
  states.insert({0, v_block_0});
  states.insert({1, v_block_1});
  states.insert({2, v_block_2});
  states.insert({3, v_block_3});

  return states;
}
void get_order_one_monomials(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int Ly, int Lx, bool use_symm)
{

  std::vector<std::string> dirs = {"x", "y", "z"};

  for (auto s : dirs)
  {

    op_vec v0 = {spin_op(s, {0, 0}, {Ly, Lx})};

    auto sign = get_sec(v0);

    if (use_symm)
    {
      add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
    }
    else
    {
      add_state(states, v0, map_sec);
    }
  }
}

void get_order_two_monomials(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int Ly, int Lx, int ry, int rx, int starty, int startx, bool use_symm)
{

  std::vector<std::string> dirs = {"x", "y", "z"};
  // for (auto s : dirs)
  // {
  //   op_vec v0 = {spin_op(s, {0, 0}, {Ly, Lx}), spin_op(s, {1, 0}, {Ly, Lx})};

  //   auto [fac, vec] = get_normal_form(v0);

  //   //           if (use_symm)
  //   //           {
  //   add_state_with_symmetries(states, vec, map_sec, Ly, Lx);
  // }

  // for (auto s : dirs)
  // {
  //   op_vec v0 = {spin_op(s, {0, 0}, {Ly, Lx}), spin_op(s, {0, 1}, {Ly, Lx})};

  //   auto [fac, vec] = get_normal_form(v0);

  //   //           if (use_symm)
  //   //           {
  //   add_state_with_symmetries(states, vec, map_sec, Ly, Lx);
  // }
  for (int i = starty; i <= ry; i++)
  {

    for (int j = startx; j <= rx; j++)
    {
      for (auto s1 : dirs)
      {

        for (auto s2 : dirs)
        {

          if (i != 0 or j != 0)
          {
            int ind1 = (Ly + i) % Ly;
            int ind2 = (Lx + j) % Lx;
            op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {ind1, ind2}, {Ly, Lx})};

            auto [fac, vec] = get_normal_form(v0);

            if (use_symm)
            {
              add_state_with_symmetries(states, vec, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, vec, map_sec);
            }
          }
        }
      }
    }
  }
}

void get_order_three_monomials(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int Ly, int Lx, bool use_symm)
{
  // op_vec v0 = {spin_op("x", {0, 0}, {Ly, Lx}), spin_op("z", {0, 1}, {Ly, Lx}), spin_op("y", {1, 1}, {Ly, Lx})};
  // auto [fac, vec] = get_normal_form(v0);

  // add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
  std::vector<std::string> dirs = {"x", "y", "z"};
  for (auto s1 : dirs)
  {
    for (auto s2 : dirs)
    {
      for (auto s3 : dirs)
      {

        {
          op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {0, 1}, {Ly, Lx}), spin_op(s3, {1, 1}, {Ly, Lx})};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
        {
          op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {0, 1}, {Ly, Lx}), spin_op(s3, {Ly - 1, 1}, {Ly, Lx})};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
        {
          op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {1, 0}, {Ly, Lx}), spin_op(s3, {1, 1}, {Ly, Lx})};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
        {
          op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {Ly - 1, 0}, {Ly, Lx}), spin_op(s3, {Ly - 1, 1}, {Ly, Lx})};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
        {
          op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {1, 0}, {Ly, Lx}), spin_op(s3, {2, 0}, {Ly, Lx})};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
        {
          op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {0, 1}, {Ly, Lx}), spin_op(s3, {0, 2}, {Ly, Lx})};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
      }
    }
  }
}
void get_order_four_monomials(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int Ly, int Lx, bool use_symm)
{
  std::vector<std::string> dirs = {"x", "y", "z"};
  for (auto s1 : dirs)
  {
    for (auto s2 : dirs)
    {
      for (auto s3 : dirs)
      {

        {
          for (auto s4 : dirs)
          {

            {
              op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {1, 0}, {Ly, Lx}), spin_op(s3, {0, 1}, {Ly, Lx}), spin_op(s4, {1, 1}, {Ly, Lx})};
              auto [fac, vec] = get_normal_form(v0);

              if (use_symm)
              {
                add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
              }
              else
              {
                add_state(states, v0, map_sec);
              }
            }
            {
              op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {1, 0}, {Ly, Lx}), spin_op(s3, {2, 0}, {Ly, Lx}), spin_op(s4, {3, 0}, {Ly, Lx})};
              auto [fac, vec] = get_normal_form(v0);

              if (use_symm)
              {
                add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
              }
              else
              {
                add_state(states, v0, map_sec);
              }
            }
          }
        }
      }
    }
  }
}
// some 1d functions
void get_order_three_monomials_1d(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int Ly, int Lx, bool use_symm)
{

  std::vector<std::string> dirs = {"x", "y", "z"};
  for (auto s1 : dirs)
  {
    for (auto s2 : dirs)
    {
      for (auto s3 : dirs)
      {

        {
          op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {0, 1}, {Ly, Lx}), spin_op(s3, {0, 2}, {Ly, Lx})};
          auto [fac, vec] = get_normal_form(v0);

          {

            if (use_symm)
            {
              add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
            }
            else
            {
              add_state(states, v0, map_sec);
            }
          }
        }
      }
    }
  }
}
void get_order_four_monomials_1d(basis_structure &states, std::map<std::pair<int, int>, int> &map_sec, int Ly, int Lx, bool use_symm)
{
  std::vector<std::string> dirs = {"x", "y", "z"};
  for (auto s1 : dirs)
  {
    for (auto s2 : dirs)
    {
      for (auto s3 : dirs)
      {

        {
          for (auto s4 : dirs)
          {

            {
              op_vec v0 = {spin_op(s1, {0, 0}, {Ly, Lx}), spin_op(s2, {0, 1}, {Ly, Lx}), spin_op(s3, {0, 2}, {Ly, Lx}), spin_op(s4, {0, 3}, {Ly, Lx})};
              auto [fac, vec] = get_normal_form(v0);

              if (use_symm)
              {
                add_state_with_symmetries(states, v0, map_sec, Ly, Lx);
              }
              else
              {
                add_state(states, v0, map_sec);
              }
            }
          }
        }
      }
    }
  }
}
