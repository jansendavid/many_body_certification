#pragma once
#include "fusion.h"
#include"spins.hpp"
#include<unordered_map>
#include <memory>
using namespace mosek::fusion;
using namespace monty;
std::shared_ptr<ndarray<int,1>>    nint(const std::vector<int> &X)    { return new_array_ptr<int>(X); }
std::shared_ptr<ndarray<double,1>> ndou(const std::vector<double> &X) { return new_array_ptr<double>(X); }
template<class T>
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
    else { 
        // If the element is not 
        // present in the vector 
        return -1; 
    } 
} 
struct mom_ref
{
    mom_ref(Variable::t var, int i1, op_vec vec):var_(var),i1_(i1),vec_(vec){}
  Variable::t var_;
   int i1_;
  op_vec vec_;
  mom_ref()=default;
};
op_vec translation_y(op_vec op, int j, int L){
  op_vec vec;
  for(int i=0; i<op.size(); i++)
    {
      vec.push_back(op[i].get_translated_y(j,L));
    }
  return vec;
}
std::vector<op_vec> generate_all_translations_y(op_vec op, int L, int inc=1)
{
  std::vector<op_vec> all_T;
  all_T.push_back(op);
  for(int i=1; i<L; i+=inc)
    {
      auto new_op=translation_y(op, i,L);
      //      auto coeff=bubbleSort(new_op, new_op.size());
      auto [fac, vec] =get_normal_form(new_op);
      if(std::abs(fac.imag())>1e-9 or std::abs(fac.real()-1)>1e-9 ){std::cout<< "error in all trans"<<std::endl;}
      all_T.push_back(vec);
	
    }
  return all_T;
}

bool is_zero_signsym(op_vec op)
{

  std::vector<std::string> dirs={"x", "y", "z"};
    for(auto dir_: dirs)
	{
	  int fac=1;
    for(auto a: op)
    {
      
      if(a.dir_==dir_){fac*=-1;}
	}
      if(fac<0)
	{return true;}
    }
    return false;

}


op_vec translation(op_vec op, int j, int L){
  op_vec vec;
  for(int i=0; i<op.size(); i++)
    {
      vec.push_back(op[i].get_translated(j,L));
    }
  return vec;
}
std::vector<op_vec> generate_all_translations(op_vec op, int L)
{
 std::vector<op_vec> all_T;
  all_T.push_back(op);
  if(print_op(op)=="1"){return all_T;}
  

  for(int i=1; i<L; i++)
    {
      auto new_op=translation(op, i,L);
      //    do we need auto
      auto [fac, vec] =get_normal_form(new_op);
      if(std::abs(fac.imag())>1e-9 or std::abs(fac.real()-1)>1e-9 ){std::cout<< "error in all trans"<<std::endl;}
      all_T.push_back(vec);
	
    }
  return all_T;
}
  void display(char a[], int n) 
{ 
  for (int i = 0; i < n; i++) { 
    std::cout << a[i] << " "; 
  } 
  std::cout << std::endl; 
}
op_vec mirror(op_vec op){
  op_vec vec;
  for(int i=0; i<op.size(); i++)
    {
      vec.push_back(op[i].get_mirror());
    }
  return vec;
}
 void findPermutations(char a[], int n) 
{ 
 
  // Sort the given array 
  std::sort(a, a + n); 
 
  // Find all possible permutations 
  std::cout << "Possible permutations are:\n"; 
  do { 
    display(a, n); 
  } while (std::next_permutation(a, a + n)); 
} 
std::vector<op_vec> generate_all_permutations_xy(op_vec op)
{
 
 
 std::vector<op_vec> all_P;
 all_P.push_back(op);


  std::vector<std::map<std::string, std::string>> permutations(1);
 permutations[0].insert({"x","y"});
 permutations[0].insert({"y","x"});
 permutations[0].insert({"z","z"});

 for(auto& a: permutations)
   {
     auto new_op=op;
     std::for_each(new_op.begin(), new_op.end(), [a](spin_op &n){
       auto old_d=n.dir_;
       n.dir_=a.at(old_d); });
     all_P.push_back(new_op);
   }
 

  return all_P;
}
std::vector<op_vec> generate_all_permutations_xyz(op_vec op)
{
 
 
 std::vector<op_vec> all_P;
 all_P.push_back(op);
 std::vector<std::map<std::string, std::string>> permutations(5);
 permutations[0].insert({"x","x"});
 permutations[0].insert({"y","z"});
 permutations[0].insert({"z","y"});

 permutations[1].insert({"x","y"});
 permutations[1].insert({"y","x"});
 permutations[1].insert({"z","z"});
 
 permutations[2].insert({"x","y"});
 permutations[2].insert({"y","z"});
 permutations[2].insert({"z","x"});

 permutations[3].insert({"x","z"});
 permutations[3].insert({"y","x"});
 permutations[3].insert({"z","y"});

 permutations[4].insert({"x","z"});
 permutations[4].insert({"y","y"});
 permutations[4].insert({"z","x"});


 for(auto& a: permutations)
   {
     auto new_op=op;
     std::for_each(new_op.begin(), new_op.end(), [a](spin_op &n){
       auto old_d=n.dir_;
       n.dir_=a.at(old_d); });
     all_P.push_back(new_op);
   }
 

  return all_P;
}
