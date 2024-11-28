#pragma once
#include "fusion.h"
#include"spins.hpp"
#include<unordered_map>
#include <memory>
using namespace mosek::fusion;
using namespace monty;
std::shared_ptr<ndarray<int,1>>    nint(const std::vector<int> &X)    { return new_array_ptr<int>(X); }
std::shared_ptr<ndarray<double,1>> ndou(const std::vector<double> &X) { return new_array_ptr<double>(X); }
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
std::pair<int,int> get_sec(op_vec op)
{
  // computes the sector of a given vector of operators
  // NOTE: all are of the form vec{S_i}vec{S_j} in the Hamiltonian is needed. For e.g., TFI, small modifications must be made
  int sxy=1;
  int syz=1;
  for(auto a: op)
    {
      if(a.dir_=="x"){sxy*=-1;}
      if(a.dir_=="y"){syz*=-1;
	sxy*=-1;}
      if(a.dir_=="z"){syz*=-1;}
      
    }
  return std::pair<int,int>(sxy, syz);
}

void add_state(basis_structure& states, op_vec op, std::map<std::pair<int,int>, int> map_sec)
{
  // adds a state to a basis
      auto [fac, nf] =get_normal_form(op);
	     auto sign=get_sec( nf);
	     if(nf.size()>0)
	       {
	     states.at(map_sec.at(sign)).push_back(nf);
	       }
	     
	    

  return;}

void add_state_symm(basis_structure& states, op_vec op, std::map<std::pair<int,int>, int> map_sec, int L, std::string permuts)
{
  // adds a state to a basis
      auto [fac, nf] =get_normal_form(op);
      auto sign=get_sec( nf);
      if(nf.size()>0)
	{
	  auto all_t=generate_all_translations(nf, L);
	  bool found=false;
	  
	     	 for(auto op_t: all_t)
	   {
	     
	     
      // 	     // generate all y translations, inc=1
	      auto all_ty=generate_all_translations_y(op_t, L,1);
	     for(auto op_ty: all_ty)
	      {
      // 	     // generate all permutations
	   // 	 std::vector<op_vec> all_p;	 
	   // 	 if(permuts=="xyz")
	   // 	   {
	   //    all_p=generate_all_permutations_xyz(op_ty);
	   // 	   }
	   // 	 else if(permuts=="xy")
	   // 	   {
	   // 	     all_p=generate_all_permutations_xy(op_ty);
	   // 	   }

	   //   for(auto op_p: all_p)
	   // { 
	     
	     int cnt = count(states.at(map_sec.at(sign)).begin(), states.at(map_sec.at(sign)).end(), op_ty);
	     if(cnt>0){
	       found=true;
	       //}
	        }
	       }
       	   }
	     if(not found)
	       {   states.at(map_sec.at(sign)).push_back(nf);}
       	   }
	     
	    

  return;}
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
  basis_structure get_basis_2d(int L, int r)
  {
std::set <std::string> states_strings;
 basis_structure states;
      std::vector<op_vec> v_block_0;
  std::vector<op_vec> v_block_1;
  std::vector<op_vec> v_block_2;
  std::vector<op_vec> v_block_3;
  states.insert({0, v_block_0});
  states.insert({1, v_block_1});
  states.insert({2, v_block_2});
  states.insert({3, v_block_3});
  std::map<std::pair<int,int>, int> map_sec;
      map_sec.insert({std::pair<int,int>(1,1), 0});
  map_sec.insert({std::pair<int,int>(1,-1), 1});
  map_sec.insert({std::pair<int,int>(-1,1), 2});
    map_sec.insert({std::pair<int,int>(-1,-1), 3});
 
       std::vector<std::string> dirs={"x", "y", "z"};


	   for(auto s: dirs){
	    
	     op_vec v0={spin_op(s, {0,0}, L)};
       const bool is_in = states_strings.find(print_op(v0)) != states_strings.end();
       if(not is_in)
       {add_state(states, v0, map_sec);}


	    }
      int SS=0;
      for(int i=0; i<=r; i++)
      {
     
          for(int j=0; j<=r; j++)
      {
        for(auto s1: dirs){
	     
       for(auto s2: dirs){
       {
        
        if(i!=0 or j!=0)
        {
          int ind1=(L+i)%L;
          int ind2=(L+j)%L;
        op_vec v0={spin_op(s1, {0,0}, L),spin_op(s2, {ind1, ind2}, L)};

	
        auto [fac, vec] =get_normal_form(v0);
         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();

       if(not is_in)
       {
    
        add_state(states, vec, map_sec);
        
        }
      

	     
        }
		 }
		   
	    }
     }
      }
      }
 
       	   for(auto s1: dirs){
	      for(auto s2: dirs){
           for(auto s3: dirs){
	     
       {
        op_vec v0={spin_op(s1, {0,0}, L),spin_op(s2, {0,1}, L),spin_op(s3, {1,1}, L)};
        auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
       if(not is_in)
       {add_state(states, vec, map_sec);}

		 }
      {
        op_vec v0={spin_op(s1, {0,0}, L),spin_op(s2, {0,1}, L),spin_op(s3, {L-1,1}, L)};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
       if(not is_in)
       {add_state(states, vec, map_sec);}
		 }
     {
          op_vec v0={spin_op(s1, {0,0}, L),spin_op(s2, {1,0}, L),spin_op(s3, {1,1}, L)};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
       if(not is_in)
       {add_state(states, vec, map_sec);}
		 }
        {
          op_vec v0={spin_op(s1, {0,0}, L),spin_op(s2, {L-1,0}, L),spin_op(s3, {L-1,1}, L)};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
       if(not is_in)
       {add_state(states, vec, map_sec);}
		 }
      {
          op_vec v0={spin_op(s1, {0,0}, L),spin_op(s2, {1,0}, L),spin_op(s3, {2,0}, L)};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
       if(not is_in)
       {add_state(states, vec, map_sec);}
		 }
        {
          op_vec v0={spin_op(s1, {0,0}, L),spin_op(s2, {0,1}, L),spin_op(s3, {0,2}, L)};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
       if(not is_in)
       {add_state(states, vec, map_sec);}
		 }
		 
	    }
        }
     }

 for(auto s1: dirs){
	      for(auto s2: dirs){
           for(auto s3: dirs){
	     
       {
            for(auto s4: dirs){
	     
       {
        op_vec v0={spin_op(s1, {0,0}, L),spin_op(s2, {1,0}, L),spin_op(s3, {0,1}, L),spin_op(s4, {1,1}, L)};
        auto [fac, vec] =get_normal_form(v0);
	     add_state(states, vec, map_sec);
		 }
            } }}}}

       return states;
  }