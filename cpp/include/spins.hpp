#pragma once
#include "vector"
#include "string"
#include <iostream>
#include <complex>
#include <cassert>
#include <map>
#include <algorithm>
using cpx = std::complex<double>;

class spin_op_parent
{
  std::string dir_;
  std::vector<int> site_;

public:
  std::vector<int> offset_;
  std::vector<int> shifts_;
  std::string exp_;
  std::string symbol_;
  bool unit = false;
  spin_op_parent(std::string dir, std::vector<int> site, std::vector<int> offset, std::string symbol) : dir_(dir), site_(site), offset_(offset), symbol_(symbol)
  {
    if (site.size() < offset.size())
    {
      std::cout << "Need sufficient offsets" << std::endl;
    }
    shifts_.push_back(1);
    for (int i = 0; i < offset_.size() - 1; i++)
    {
      shifts_.push_back(shifts_.back() * offset_[i]);
    }
    compute_expression();
  };
  void compute_expression()
  {
    exp_ = symbol_;

    exp_ += "_[";

    exp_ += dir_ + ",(";
    for (int i = 0; i < site_.size() - 1; i++)
    {
      exp_ += std::to_string(site_[i]) + ",";
    }
    exp_ += std::to_string(site_[site_.size() - 1]);
    exp_ += ")]";
  };
  std::vector<int> get_site() const
  {
    return site_;
  }
  void set_site(std::vector<int> new_site)
  {
    site_ = new_site;
    compute_expression();
  }
  template <typename T>
  void set_dir(T &new_dir)
  {
    dir_ = new_dir;
    compute_expression();
  }
  std::string get_dir() const
  {
    return dir_;
  }
  spin_op_parent() { unit = true; };
  // can this function be improved?
  int const pos() const
  {

    int position = 0;
    for (int i = 0; i < offset_.size(); i++)
    {
      // std::cout<< "offset "<<offset_[i]<< " and factor "<<offset_.size()-1<<std::endl;
      position += site_[i] * shifts_[i];
    }

    return position;
  }
  spin_op_parent get_mirror_parent()
  {
    // assumes LxL lattice and order ....y,x
    auto new_sites = site_;

    new_sites[new_sites.size() - 1] = site_[site_.size() - 2];
    new_sites[new_sites.size() - 2] = site_[site_.size() - 1];

    return spin_op_parent(dir_, new_sites, offset_, symbol_);
  }
  spin_op_parent get_flipped_layer_parent()
  {
    // assumes LxL lattice and order ....y,x
    assert(site_.size() == 3);
    auto new_sites = site_;

    new_sites[0] = (site_[0] + 1) % 2;

    return spin_op_parent(dir_, new_sites, offset_, symbol_);
  }
  spin_op_parent get_translated_parent(int j, int L)
  {
    auto new_sites = site_;
    new_sites[new_sites.size() - 1] = (new_sites[new_sites.size() - 1] + j) % L;

    return spin_op_parent(dir_, new_sites, offset_, symbol_);
  }

  spin_op_parent get_translated_y_parent(int j, int L)
  {
    auto new_sites = site_;

    new_sites[new_sites.size() - 2] = (new_sites[new_sites.size() - 2] + j) % L;

    return spin_op_parent(dir_, new_sites, offset_, symbol_);
  }

  std::string expression() const
  {
    return exp_;
  }
  bool operator<(const spin_op_parent &obj) const
  {

    return this->expression() < obj.expression();
  }
  bool operator>(const spin_op_parent &obj) const
  {
    return this->expression() > obj.expression();
  }
  bool operator==(const spin_op_parent &obj)
  {

    return (expression() == obj.expression());
  }

  friend bool operator==(const spin_op_parent &c1, const spin_op_parent &c2);
  friend bool operator!=(const spin_op_parent &c1, const spin_op_parent &c2);
};
bool operator==(const spin_op_parent &c1, const spin_op_parent &c2)
{
  return (c1.expression() == c2.expression());
}

bool operator!=(const spin_op_parent &c1, const spin_op_parent &c2)
{
  return (c1.expression() != c2.expression());
}

std::ostream &operator<<(std::ostream &os, const spin_op_parent &op)
{
  os << op.expression();
  return os;
}
class spin_op : public spin_op_parent
{
public:
  using spin_op_parent::spin_op_parent;
  spin_op(std::string dir, std::vector<int> site, std::vector<int> offset) : spin_op_parent(dir, site, offset, "s") {}
  spin_op get_translated_y(int j, int L)
  {
    auto base = get_translated_y_parent(j, L);
    return spin_op(base.get_dir(), base.get_site(), base.offset_);
  }
  spin_op get_translated(int j, int L)
  {
    auto base = get_translated_parent(j, L);
    return spin_op(base.get_dir(), base.get_site(), base.offset_);
  }
  spin_op get_flipped_layer()
  {
    auto base = get_flipped_layer_parent();
    return spin_op(base.get_dir(), base.get_site(), base.offset_);
  }
  spin_op get_mirror()
  {
    auto base = get_mirror_parent();
    return spin_op(base.get_dir(), base.get_site(), base.offset_);
  }
};

/////////////////////////////////////////////////////////////////////////
using op_vec = std::vector<spin_op>;
using sector_structure = std::map<int, std::vector<op_vec>>;
using basis_structure = std::map<int, sector_structure>;

/////////////////////////////////////////////////////////////////
template <typename T>
std::string print_op(const std::vector<T> &oper)
{
  if (oper.empty())
    return "1";

  // Estimate total size to reserve
  size_t total_size = 0;
  for (const auto &O : oper)
    total_size += O.expression().size(); // sum of lengths

  std::string s;
  s.reserve(total_size); // avoid reallocations

  for (const auto &O : oper)
    s += O.expression();

  return s;
}

template <typename T>
std::vector<T> dagger_operator(std::vector<T> oper)
{
  double coeff = 1;
  std::vector<T> new_op;
  reverse(oper.begin(), oper.end());
  return oper;
}
template <typename T>
T apply_commutator(T &arr, int i)
{
  // if empty vector, then nothing happens (already c^dagc)
  // else change arr and return with deleted element

  T arr_copy;
  copy(arr.begin(), arr.begin() + i, back_inserter(arr_copy));

  copy(arr.begin() + i + 2, arr.end(), back_inserter(arr_copy));
  std::swap(arr[i], arr[i + 1]);
  return arr_copy;
}

std::pair<cpx, std::string> get_direction(std::string a, std::string b)
{
  std::pair<cpx, std::string> p1(cpx(0, 0), std::string("0"));
  if (a == "x" and b == "y")
  {
    return {cpx(0, 1), std::string("z")};
  }
  else if (a == "y" and b == "x")
  {
    return {cpx(0, -1), std::string("z")};
  }
  else if (a == "z" and b == "x")
  {
    return {cpx(0, 1), std::string("y")};
  }
  else if (a == "x" and b == "z")
  {
    return {cpx(0, -1), std::string("y")};
  }
  else if (a == "y" and b == "z")
  {
    return {cpx(0, 1), std::string("x")};
  }
  else if (a == "z" and b == "y")
  {
    return {cpx(0, -1), std::string("x")};
  }
  else
  {
    std::cout << "error commutators " << std::endl;
    return {cpx(0, 0), std::string("0")};
  }
}
std::pair<cpx, op_vec> run_loop(op_vec op)
{
  cpx pref = 1;
  std::map<int, spin_op> copied;
  for (int i = 0; i < op.size(); i++)
  {
    copied.insert({i, op[i]});
  }
  for (int i = 0; i < op.size() - 1; i++)
  {
    auto it = copied.find(i);

    if (it != copied.end())
    {
      auto o1 = copied[i];
      auto o2 = copied[i + 1];
      if (o1.pos() != o2.pos())
      {
        continue;
      }
      else if (o1.get_dir() == o2.get_dir())
      {
        copied.erase(i);
        copied.erase(i + 1);
      }
      else
      {
        auto [new_coeff, new_op] = get_direction(o1.get_dir(), o2.get_dir());
        copied.erase(i);
        copied.erase(i + 1);
        copied.insert({i, spin_op(new_op, o1.get_site(), o1.offset_)});
        pref *= new_coeff;
      }
    }
  }

  op_vec final_vec;
  for (auto it = copied.begin(); it != copied.end(); ++it)
  {
    final_vec.push_back(it->second);
  }

  return {pref, final_vec};
}

std::pair<cpx, op_vec> get_normal_form(op_vec op)
{
  if (op.size() < 1)
  {
    return std::pair<cpx, op_vec>(1.0, {});
    // std::cout << "error: unit sent to normal form" << std::endl;
  }

  cpx pref(1., 0);
  bool change = false;
  auto new_list = op;
  auto new_list_x = op;
  // deleted
  // sort(new_list_x.begin(), new_list_x.end(), [](spin_op &o1, spin_op &o2)
  //      { return o1.pos() < o2.pos(); });
  std::stable_sort(new_list.begin(), new_list.end(),
                   [](const spin_op &lhs, const spin_op &rhs)
                   {
                     return lhs.pos() < rhs.pos();
                   });
  // if (print_op(new_list) != print_op(new_list_x))
  // {
  //   std::cout << "sorting issue" << std::endl;
  // }
  // for (auto a : new_list)
  // {
  //   // std::cout<<a.pos()<<std::endl;
  // }

  while (not change)
  {

    auto [new_coeff, new_op] = run_loop(new_list);
    pref *= new_coeff;
    if (new_op.size() == 0 or new_op == new_list)
    {
      new_list = new_op;
      change = true;
    }
    else
    {
      new_list = new_op;
    }
  }
  return std::pair<cpx, op_vec>{pref, new_list};
}

std::vector<std::pair<cpx, op_vec>> generate_all_terms(op_vec op, bool &unit_found)
{
  std::vector<std::pair<cpx, op_vec>> terms;

  return terms;
}
std::pair<std::complex<double>, op_vec> sdp_get_form(op_vec op)
{

  auto [coeff_, nf] = get_normal_form(op);
  // std::cout << "fac " << coeff_ << std::endl;
  assert(std::abs(coeff_) - 1 < 1e-8);
  return std::pair<std::complex<double>, op_vec>(coeff_, nf);
}