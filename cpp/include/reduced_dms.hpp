#pragma once
#include <vector>
// Import Eigen
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
using mat_type = Eigen::MatrixXcd;
std::string find_state(op_vec nf, std::map<std::string, std::pair<std::string, std::complex<double>>> mp, int Lx)
{
  // pass in normal form

  auto all_t = generate_all_translations(nf, Lx);
  bool found = false;

  for (auto op_t : all_t)
  {

    // 	     // generate all y translations, inc=1
    auto all_ty = generate_all_translations_y(op_t, Lx, 1);

    for (auto op_ty : all_ty)
    {

      auto key = print_op(op_ty);
      auto it = mp.find(key);

      if (it != mp.end())
      {
        return key;
      }
      else
      {
      }
    }
  }
  std::cout << "RDM not found " << print_op(nf) << std::endl;
  std::string s = "";
  return s;
}
Matrix::t get_sparse_from_eigen(Eigen::MatrixXcd mat)
{
  std::vector<double> values;
  std::vector<int> rows;
  std::vector<int> cols;

  for (int i = 0; i < mat.rows(); i++)
  {

    for (int j = 0; j < mat.cols(); j++)
    {

      if (std::abs(mat.coeff(i, j)) > 1e-9)
      {
        // std::cout<< "her"<<std::endl;
        if (std::abs(mat.coeff(i, j).real()) > 1e-9)
        {
          values.push_back(mat.coeff(i, j).real());
          rows.push_back(i); // row index
          cols.push_back(j); // col index (here it is equal to k)
          values.push_back(mat.coeff(i, j).real());
          rows.push_back(i + mat.cols()); // row index
          cols.push_back(j + mat.cols()); // col index (here it is equal to k)
                                          // it.index(); // inner index, here it is equal to it.row()
        }

        if (std::abs(mat.coeff(i, j).imag()) > 1e-9)
        {
          values.push_back(mat.coeff(i, j).imag());
          rows.push_back(i + mat.cols()); // row index
          cols.push_back(j);              // col index (here it is equal to k) // correct with minus?
          values.push_back(-1. * mat.coeff(i, j).imag());
          rows.push_back(j + mat.cols()); // row index
          cols.push_back(i + mat.cols());
        }
      }
    }
  }

  Matrix::t Alpha = Matrix::sparse(2 * mat.rows(), 2 * mat.cols(), nint(rows), nint(cols), ndou(values));
  return Alpha;
}
// xt::xarray<std::complex<double>, xt::layout_type::dynamic>;
void generate_rmds_primal(std::vector<std::pair<int, int>> sites, std::map<std::string, int> refs, std::map<std::string, std::pair<std::string, std::complex<double>>> map, Variable::t var, int Lx, Model::t M)
{

  // mat_type sigz({2,2});
  std::cout << "generate RDM" << std::endl;
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

  // pauliX.makeCompressed();

  mat_type pauliY = mat_type::Zero(2, 2);
  pauliY(0, 1) = std::complex<double>(0, -1);
  pauliY(1, 0) = std::complex<double>(0, 1);

  // pauliX.makeCompressed();

  // std::cout<< sigy<<std::endl;
  int degree = sites.size();
  std::vector<std::string> terms;
  std::map<std::string, mat_type> sigma_map;
  sigma_map.insert({"1", pauliI});
  sigma_map.insert({"x", pauliX});
  sigma_map.insert({"y", pauliY});
  sigma_map.insert({"z", pauliZ});

  auto dirs = std::vector<std::string>{"1", "x", "y", "z"};
  std::vector<std::string> tots;

  for (auto d1 : dirs)
  {

    if (degree == 1)
    {
      tots.push_back(d1);
      continue;
    }
    for (auto d2 : dirs)
    {
      if (degree == 2)
      {
        tots.push_back(d1 + d2);
        continue;
      }
      for (auto d3 : dirs)
      {
        if (degree == 3)
        {
          tots.push_back(d1 + d2 + d3);
          continue;
        }
        for (auto d4 : dirs)
        {
          if (degree == 4)
          {
            tots.push_back(d1 + d2 + d3 + d4);
            continue;
          }
        }
      }
    }
  }
  int dim = std::pow(2, degree);

  std::vector<Expression::t> matrices;

  double prefac = 1;
  std::cout << "tots ;enght " << tots.size() << std::endl;
  for (auto t : tots)
  {
    mat_type mat;
    op_vec state;

    for (int i = 0; i < t.size(); i++)
    {

      std::string key = t.substr(i, 1);

      if (key != "1")
      {
        state.push_back(spin_op(key, {0, i}, {Lx, Lx}));
      }
      if (i == 0)
      {
        mat = sigma_map[key];
      }
      else
      {
        // std::cout<< "mat 2"<<sigma_map[key]<<std::endl;
        mat = Eigen::KroneckerProduct(mat, sigma_map[key]).eval();
      }
    }

    if (print_op(state) == "1")
    {
      // std::cout<< "mat 1"<<mat<<std::endl;
      mat = mat / std::
                      pow(2, degree);
      auto Alpha = get_sparse_from_eigen(mat);

      auto el = refs.at(print_op(state));
      matrices.push_back(Expr::mul(var->index(el), Alpha));
    }
    else
    {

      auto [fac, nf] = get_normal_form(state);
      std::string state_str = find_state(nf, map, Lx);
      // auto state_key=
      // find_state(nf,refs,map,Lx );

      auto [state_from_map, coeff] = map.at(state_str);
      if (state_from_map == "0")
      {
        continue;
      }
      if (std::abs((fac * coeff).imag()) > 1e-9)
      {
        std::cout << "value error" << std::endl;
      }
      auto el = refs.at(state_from_map);

      mat = mat * (fac * coeff).real() / std::pow(2, degree);
      auto Alpha = get_sparse_from_eigen(mat);
      //          //if(fac.imag())

      matrices.push_back(Expr::mul(Alpha, var->index(el)));
    }
  }
  M->constraint((Expr::add(new_array_ptr(matrices))), Domain::inPSDCone());
  return;
}
