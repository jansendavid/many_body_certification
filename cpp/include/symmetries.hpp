#pragma once
#include "fusion.h"
#include"spins.hpp"
#include<unordered_map>
#include <memory>
using namespace mosek::fusion;
using namespace monty;
std::shared_ptr<ndarray<int,1>>    nint(const std::vector<int> &X)    { return new_array_ptr<int>(X); }
std::shared_ptr<ndarray<double,1>> ndou(const std::vector<double> &X) { return new_array_ptr<double>(X); }
struct rdm_operator{
  
  std::vector<std::vector<int>> op_;
  rdm_operator(std::vector<std::vector<int>> op): op_(op){};
rdm_operator(){};
unsigned int size()
{
  return op_.size();
}
std::vector<int> at(int i)
{
  return op_[i];
}
 bool operator <(const rdm_operator& rhs) const
    {
        return op_ < rhs.op_;
    }

};
struct rdms_struct{
  // struct to manage the reduced density matrices
  // store a vector with operator and each opertor has vector of its sites
std::vector<rdm_operator> rdms;
rdms_struct(std::vector<rdm_operator> rdms): rdms(rdms){};
rdms_struct(){};
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
rdms_struct  get_rdms(int Lx, int dim)
{
  rdms_struct data;
     
        rdm_operator newstate({{0,0}, {0,1}});
data.add_operator(newstate);
if(dim>=4)
{
    rdm_operator newstate({{0,0}, {0,1}, {0,2}});
data.add_operator(newstate);
 rdm_operator newstate_1({{0,0}, {0,1}, {0,2},{0,3}});
data.add_operator(newstate_1);
 rdm_operator newstate_2({{0,0}, {0,1}, {0,2},{0,3},{1,1}});
data.add_operator(newstate_2);
//  rdm_operator newstate_3({{0,0}, {0,1}, {0,2},{0,3},{1,1},{1,2}});
// data.add_operator(newstate_3);
//  rdm_operator newstate_4({{0,0}, {0,1}, {0,2},{0,3},{1,1},{1,2},{1,3}});
// data.add_operator(newstate_4);
//  rdm_operator newstate_5({{0,0}, {0,1}, {0,2},{0,3},{1,1},{1,2},{1,3},{2,2}});
// data.add_operator(newstate_5);

//  rdm_operator newstate_55({{0,0}, {2,1}, {0,2},{2,3},{1,1},{1,2},{1,3},{2,2}});
//data.add_operator(newstate_55);
//  rdm_operator newstate_6({{0,0}, {0,1}, {0,2},{0,3},{1,1},{1,2},{1,3},{2,2},{3,3}});
// data.add_operator(newstate_6);
//  rdm_operator newstate_7({{0,0}, {0,1}, {0,2},{0,3},{1,1},{1,2},{1,3},{2,2},{3,3}});
// data.add_operator(newstate_7);
}
if(dim>=6)
{
    rdm_operator newstate({{0,0}, {0,1}, {0,2},{0,3},{0,4}});
data.add_operator(newstate);
 rdm_operator newstate_1({{0,0}, {0,1}, {0,2},{0,3},{0,4},{0,5}});
data.add_operator(newstate_1);
}
  
return data;

}
bool found_operator(std::set<std::string> list_of_op, op_vec vec, int Lx, int Ly)
{
  bool results=false;
  //   auto all_ty=generate_all_translations_y(op, Ly,1);
		  
	
	//  for(auto op_ty: all_ty)
	//    {
	//     auto all_t=generate_all_translations(op_ty, Lx);
	
	//      for(auto op_t: all_t)
	//        {
	// 		auto it=mat_terms.find(print_op(op_t));
	// 		if(it != mat_terms.end())
	// 		{
				
	//     TI_map_.insert({print_op(op), { it->first,1}});
	// 	return;
	// 		}
	// 		std::vector<op_vec> all_p;	
	// 		if(permuts_=="xyz" or permuts_=="yxz" or permuts_=="zxy" or permuts_=="zyx")
	// 	   {
	//       all_p=generate_all_permutations_xyz(op_t);
	// 	   }
	// 	 else if(permuts_=="xy")
	// 	   {
	// 	     all_p=generate_all_permutations_xy(op_t);
	// 	   }
	// 	        	 for(auto op_p: all_p)
	//    { 
	//      auto all_d8sym=generate_all_d8(op_p,  L_);
		
	// 	 for(auto d8s: all_d8sym)
	// 	 {	
	//        auto it=mat_terms.find(print_op(d8s));
	// 	  if(it != mat_terms.end())
	// 	 	{
				
	//     TI_map_.insert({print_op(op), { it->first,1}});
	// 	 return;
	// 	 	}
		
	// 		      auto op_mirror=mirror(d8s);
				
	//        it=mat_terms.find(print_op(op_mirror));
	//        if(it != mat_terms.end())
	//        {
	// 		TI_map_.insert({print_op(op), { it->first,1}});
	// 		 return;
	//    		}
	       
	// 	   }
			
		
	//    }
	//    }
  //    }

  return results;
}
rdms_struct  translation_invariant_rdms_2nd(int Lx, int Ly)
{
  rdms_struct data;
  std::string s="x";
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
  for(int i=0; i<Lx; i++)
  {
    for(int j=0; j<Lx; j++)
    {
      if(i!=0 and j!=0)
      {
      combs.push_back({i,j});
    }
    }
  }
    for(auto p: combs)
    {
         for(auto o: combs)
    {
        for(auto l: combs)
    {
      if(p!=o and p!=l and o!=l)
      {
        rdm_operator newstate({{0,0}, {o}, {p}, {l}});
    data.add_operator(newstate);
      }
    }
    }
    }


 return data;

}
rdms_struct  translation_invariant_rdms_4th(int Lx, int Ly)
{
  rdms_struct  data;
  std::vector<std::vector<int>> combs;
  for(int i=0; i<Lx; i++)
  {
    for(int j=0; j<Lx; j++)
    {
      if(i!=0 and j!=0)
      {
      combs.push_back({i,j});
    }
    }
  }
    for(auto p: combs)
    {
         for(auto o: combs)
    {
        for(auto l: combs)
    {
          for(auto m: combs)
    {
      std::set<std::vector<int>> set_={p,o,l,m};

      if(set_.size()==4)
      {
        rdm_operator newstate({{0,0}, {o}, {p}, {l}, {m}});
    data.add_operator(newstate);
      }
    }
    }
    }
    }

return data;

}
// std::vector<std::vector<std::vector<int>>>  translation_invariant_rdms_5th(int Lx, int Ly)
// {
//   std::vector<std::vector<std::vector<int>>> data;
//   std::vector<std::vector<int>> combs;
//   for(int i=0; i<Lx; i++)
//   {
//     for(int j=0; j<Lx; j++)
//     {
//       if(i!=0 and j!=0)
//       {
//       combs.push_back({i,j});
//     }
//     }
//   }
//     for(auto m: combs)
//     {
//          for(auto l: combs)
//     {
//         for(auto n: combs)
//     {
//           for(auto o: combs)
//     {
//               for(auto p: combs)
//     {
//       std::set<std::vector<int>> set_={m,l,n,o,p};

//       if(set_.size()==5)
//       {
//         std::vector<std::vector<int>> newstate={{0,0}, {m}, {l}, {n}, {o}, {p}};
//     data.push_back(newstate);
//       }
//     }
//     }
//     }
//     }
//     }

// return data;

// }
// std::vector<std::vector<std::vector<int>>>  translation_invariant_rdms_6th(int Lx, int Ly)
// {
//   std::vector<std::vector<std::vector<int>>> data;
//   std::vector<std::vector<int>> combs;
//   for(int i=0; i<Lx; i++)
//   {
//     for(int j=0; j<Lx; j++)
//     {
//       if(i!=0 and j!=0)
//       {
//       combs.push_back({i,j});
//     }
//     }
//   }
//     for(auto m: combs)
//     {
//          for(auto l: combs)
//     {
//         for(auto n: combs)
//     {
//           for(auto o: combs)
//     {
//               for(auto p: combs)
//     {
//                  for(auto q: combs)
//     {
//       std::set<std::vector<int>> set_={m,l,n,o,p, q};

//       if(set_.size()==6)
//       {
//         std::vector<std::vector<int>> newstate={{0,0}, {m}, {l}, {n}, {o}, {p},{q}};
//     data.push_back(newstate);
//       }
//     }
//     }
//     }
//     }
//     }
//     }

// return data;

// }
op_vec apply_group_trafo(op_vec op,Eigen::Matrix2i mat, int L )
{

   op_vec vec;
  for(int i=0; i<op.size(); i++)
    {
      
      int x_cor=op[i].site_[op[i].site_.size()-1];
      int y_cor=op[i].site_[op[i].site_.size()-2];
      if(x_cor>L/2){
        x_cor=-L+x_cor;
      }
      if(y_cor>L/2){
        y_cor=-L+y_cor;
      }

      Eigen::Vector2i location(x_cor, y_cor);

      auto new_location=mat*location;
      int x_cor_new=(L+new_location(0))%L;
      int y_cor_new=(L+new_location(1))%L;

      auto sites_new=op[i].site_;
      sites_new[op[i].site_.size()-1]=x_cor_new;
      sites_new[op[i].site_.size()-2]=y_cor_new;
      vec.push_back(spin_op(op[i].dir_,sites_new, op[i].offset_));
  
    }
  return vec;


}

std::vector<op_vec> generate_all_d8(op_vec op, int L)
{
std::vector<op_vec> all_d8;
Eigen::Matrix2i a(2,2);
a(0,0)=0;
a(0,1)=-1;
a(1,0)=1;
a(1,1)=0;

Eigen::Matrix2i x;
x(0,0)=1;
x(0,1)=0;
x(1,0)=0;
x(1,1)=-1;


   auto [fac, vec] =get_normal_form(op);
   // unit element
  all_d8.push_back(vec);
  // a (rotation pi/2)
all_d8.push_back(apply_group_trafo(vec,a, L ));
//a^2
all_d8.push_back(apply_group_trafo(vec,a*a, L ));

// a^3
all_d8.push_back(apply_group_trafo(vec,a*a*a, L ));


// // x
all_d8.push_back(apply_group_trafo(vec,x, L ));

// ax
all_d8.push_back(apply_group_trafo(vec,a*x, L ));

//a^2x

all_d8.push_back(apply_group_trafo(vec,a*a*x, L ));

//a^3x
all_d8.push_back(apply_group_trafo(vec,a*a*a*x, L ));



return all_d8;
}

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
    if(print_op(op)=="1"){
      all_T.push_back(op);
    return all_T;
    }
    
		

   auto [fac, vec] =get_normal_form(op);
  all_T.push_back(vec);
  for(int i=1; i<L; i+=inc)
    {
      auto new_op=translation_y(op, i,L);
     
		
      //      auto coeff=bubbleSort(new_op, new_op.size());
      auto [fac, vec] =get_normal_form(new_op);
   
    
      if(std::abs(fac.imag())>1e-9 or std::abs(fac.real()-1)>1e-9 ){
       
        std::cout<< "error in all trans"<<std::endl;
      
        }
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

  if(print_op(op)=="1"){
      all_T.push_back(op);
    return all_T;
    }

   auto [fac, vec] =get_normal_form(op);
  all_T.push_back(vec);


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
  void add_state_with_symmetries_1d(basis_structure& states, op_vec op, std::map<std::pair<int,int>, int> map_sec, int L)
{
  // adds a state to a basis

      auto [fac, nf] =get_normal_form(op);
	     auto sign=get_sec( nf);
       bool found=false;
	     if(nf.size()>0)
	       {
      auto all_t=generate_all_translations(nf, L);
	  bool found=false;
	  
	     	 for(auto op_t: all_t)
	   {
   
        auto it = std::find(states.at(map_sec.at(sign)).begin(), states.at(map_sec.at(sign)).end(), op_t);
    if (it!=states.at(map_sec.at(sign)).end())
        {found=true;
        break;
        }
       

        
     }
          if(!found)
     {states.at(map_sec.at(sign)).push_back(nf);}
         }


      

   
     
        

  return;}

  void add_state_with_symmetries(basis_structure& states, op_vec op, std::map<std::pair<int,int>, int> map_sec, int L)
{
  // adds a state to a basis

      auto [fac, nf] =get_normal_form(op);
      bool print=false;
      
	     auto sign=get_sec( nf);
       bool found=false;
	     if(nf.size()>0)
	       {
      auto all_t=generate_all_translations(nf, L);
	  bool found=false;
	  
	     	 for(auto op_t: all_t)
	   {
  
        auto all_ty=generate_all_translations_y(op_t, L,1);

        	     for(auto op_ty: all_ty)
	      {
          if(print)
          {
            std::cout<<print_op(op_ty)<<std::endl;
          }
        auto it = std::find(states.at(map_sec.at(sign)).begin(), states.at(map_sec.at(sign)).end(), op_ty);
    if (it!=states.at(map_sec.at(sign)).end())
        {
          //std::cout<<print_op(nf)<< " was "<< print_op(*it)<<std::endl;
          found=true;
        break;}
        }
       

        
	     

     }
   

    //     auto it = std::find(states.at(map_sec.at(sign)).begin(), states.at(map_sec.at(sign)).end(), op_ty);
    // if (it!=states.at(map_sec.at(sign)).end())
    //     {std::cout<<print_op(nf)<< " was "<< print_op(*it)<<std::endl;
    //       found=true;
    //     break;}
       

    //     }
	     

     
      if(!found)
     {states.at(map_sec.at(sign)).push_back(nf);}
         
         
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
	{
 
    return true;}
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
     auto [fac, nf] =get_normal_form(vec);
    
  
  assert(fac.imag()<1e-9);
  return nf;
}
op_vec flip_layer(op_vec op){
  op_vec vec;
  for(int i=0; i<op.size(); i++)
    {
      vec.push_back(op[i].get_flipped_layer());
    }
     auto [fac, nf] =get_normal_form(vec);
    
  
  assert(fac.imag()<1e-9);
  return nf;
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
  auto [fac, vec] =get_normal_form(op);
  assert(fac.imag()<1e-9);
 all_P.push_back(vec);
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
         auto [fac, vec] =get_normal_form(new_op);
    assert(fac.imag()<1e-9);
     all_P.push_back(vec);
   }
 

  return all_P;
}
//   basis_structure get_basis_1d(int L, int r, int start=0)
//   {
// std::set <std::string> states_strings;
//  basis_structure states;
//       std::vector<op_vec> v_block_0;
//   std::vector<op_vec> v_block_1;
//   std::vector<op_vec> v_block_2;
//   std::vector<op_vec> v_block_3;
//   states.insert({0, v_block_0});
//   states.insert({1, v_block_1});
//   states.insert({2, v_block_2});
//   states.insert({3, v_block_3});
//   std::map<std::pair<int,int>, int> map_sec;
//       map_sec.insert({std::pair<int,int>(1,1), 0});
//   map_sec.insert({std::pair<int,int>(1,-1), 1});
//   map_sec.insert({std::pair<int,int>(-1,1), 2});
//     map_sec.insert({std::pair<int,int>(-1,-1), 3});
 
//        std::vector<std::string> dirs={"x", "y", "z"};
       


// 	   for(auto s: dirs){
	    
// 	     op_vec v0={spin_op(s, {0})};
//        const bool is_in = states_strings.find(print_op(v0)) != states_strings.end();
//        if(not is_in)
//        {
//         //add_state_with_symmetries_1d(states, v0, map_sec, L);
//        add_state(states, v0, map_sec);
//        }


// 	    }
//       int SS=0;
//       for(int i=1; i<=r; i++)
//       {
     
  
//         for(auto s1: dirs){
	     
//        for(auto s2: dirs){
//        {
        
      
//         op_vec v0={spin_op(s1, {0}),spin_op(s2, {i})};

	
//         auto [fac, vec] =get_normal_form(v0);
//          const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();

//        if(not is_in)
//        {
    
//         //add_state_with_symmetries_1d(states, vec, map_sec, L);
//         add_state(states, vec, map_sec);
        
//         }
      
// 	    }
//      }
//       }
//       }
 
//        	   for(auto s1: dirs){
// 	      for(auto s2: dirs){
//            for(auto s3: dirs){
	     
//        {
//         op_vec v0={spin_op(s1, {0}),spin_op(s2, {1}),spin_op(s3, {2})};
//         auto [fac, vec] =get_normal_form(v0);
// 	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
//        if(not is_in)
//        {
//         //add_state_with_symmetries_1d(states, vec, map_sec, L);
//          add_state(states, vec, map_sec);
//        }

// 		 }

		 
// 	    }
//         }
//      }

//  for(auto s1: dirs){
// 	      for(auto s2: dirs){
//            for(auto s3: dirs){
	     
//        {
//             for(auto s4: dirs){
	     
//        {
//         op_vec v0={spin_op(s1, {0}, L),spin_op(s2, {1}, L),spin_op(s3, {2}, L),spin_op(s4, {3}, L)};
//         auto [fac, vec] =get_normal_form(v0);
// 	     //add_state_with_symmetries_1d(states, vec, map_sec, L);
//        add_state(states, vec, map_sec);
// 		 }
//             } }}}}

//        return states;
//   }



  basis_structure get_basis_2d(int L, int r, int start,bool use_symm)
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
	    
	     op_vec v0={spin_op(s, {0,0}, {L,L})};
       const bool is_in = states_strings.find(print_op(v0)) != states_strings.end();
       //if(not is_in)
       {
        if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
        
        
       }


	    }
   
      
 
       	   for(auto s1: dirs){
	      for(auto s2: dirs){
           for(auto s3: dirs){
	     
       {
        op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {0,1},{L,L}),spin_op(s3, {1,1}, {L,L})};
        auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
      // if(not is_in)
       {
      
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
		 }
        }
      {
        op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {0,1}, {L,L}),spin_op(s3, {L-1,1}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
      // if(not is_in)
       {
     
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
       }
		 }
     {
          op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {1,0}, {L,L}),spin_op(s3, {1,1}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
      // if(not is_in)
       {
      
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
       }
		 }
        {
          op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {L-1,0}, {L,L}),spin_op(s3, {L-1,1}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
       //if(not is_in)
       {
       
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
        }
		 }
      {
          op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {1,0}, {L,L}),spin_op(s3, {2,0}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
       //if(not is_in)
       {
      
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
        }
		 }
        {
          op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {0,1}, {L,L}),spin_op(s3, {0,2}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	         const bool is_in = states_strings.find(print_op(vec)) != states_strings.end();
       //if(not is_in)
       {
      
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
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
            for(auto s4: dirs){
	     
       {
        op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {1,0}, {L,L}),spin_op(s3, {0,1}, {L,L}),spin_op(s4, {1,1}, {L,L})};
        auto [fac, vec] =get_normal_form(v0);
	 
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
		 }
            } }}}}

       return states;
  }


//    basis_structure get_basis_bilayer_2d(int Lx, int Ly, int layers, int r, int start,bool use_symm)
//   {
// std::set <std::string> states_strings;
//  basis_structure states;
//       std::vector<op_vec> v_block_0;
//   std::vector<op_vec> v_block_1;
//   std::vector<op_vec> v_block_2;
//   std::vector<op_vec> v_block_3;
//   states.insert({0, v_block_0});
//   states.insert({1, v_block_1});
//   states.insert({2, v_block_2});
//   states.insert({3, v_block_3});
//   std::map<std::pair<int,int>, int> map_sec;
//       map_sec.insert({std::pair<int,int>(1,1), 0});
//   map_sec.insert({std::pair<int,int>(1,-1), 1});
//   map_sec.insert({std::pair<int,int>(-1,1), 2});
//     map_sec.insert({std::pair<int,int>(-1,-1), 3});
 
//        std::vector<std::string> dirs={"x","y","z"};
       
 
// for(int i=0; i<2;i++)
// {
// 	   for(auto s: dirs){
	    
// 	     op_vec v0={spin_op(s, {i,0,0}, {layers,Ly, Lx})};
//        const bool is_in = states_strings.find(print_op(v0)) != states_strings.end();

//        {
//         if(use_symm)
//         {
// add_state_with_symmetries(states, v0, map_sec, Ly);
//         }
//         else{
//           add_state(states, v0, map_sec);
//         }
        
        
//        }


    
  
		 

struct monomial_struct_2d
{
int Lx_;
int Ly_;
// vector where all Lx*Ly monomials are stored
std::vector<std::vector<op_vec>> monomials;
monomial_struct_2d(op_vec first_op, int Lx, int Ly): Lx_(Lx), Ly_(Ly)
{
	for(int i=0;i<Ly_; i++)
	{
monomials.push_back({});
	for(int j=0;j<Ly_; j++)
	{
		monomials[i].push_back({});
	}	
	}
	monomials[0][0]=first_op;
}
op_vec get(int i, int j)
	{return monomials[i][j];}

};
bool see_if_found(std::pair<int,int> pos, std::set<std::pair<int, int>> elements, int Lx, int Ly)
{
      for(int l=0; l<Ly; l++)
           {
            for(int m=0; m<Lx; m++)
            {
            pos.first=(pos.first+l)%Ly;
           // pos.second=(pos.second+m)%Lx;
     if(elements.find(pos)!=elements.end())
     {
      return true;
     }
            }
           }
           return false;
}


std::map<std::pair<int,int>, int> get_sector_map()
{
 std::map<std::pair<int,int>, int> map_sec;
      map_sec.insert({std::pair<int,int>(1,1), 0});
  map_sec.insert({std::pair<int,int>(1,-1), 1});
  map_sec.insert({std::pair<int,int>(-1,1), 2});
    map_sec.insert({std::pair<int,int>(-1,-1), 3});
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
void get_order_one_monomials(basis_structure& states, std::map<std::pair<int,int>, int>& map_sec, int L,bool use_symm)
{

   std::vector<std::string> dirs={"x","y","z"};
       
 

	   for(auto s: dirs){
	    
 	     op_vec v0={spin_op(s, {0,0}, {L,L})};

         auto sign=get_sec( v0);

	         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
 


 	    }
}

void get_order_two_monomials(basis_structure& states, std::map<std::pair<int,int>, int>& map_sec,int L, int r, int start,bool use_symm)
{
std::vector<std::string> dirs={"x","y","z"};
   //    dirs={"x", "y", "z"};
    //  std::vector<std::pair<int, int>> rvals={{1,0}, {0,1}, {1,1}, {1,3}, {2,0}, {0,2}, {2,1}, {1,2}, {2,2}, {2,3}};
    //  for(auto s1: dirs){
	     
    //    for(auto s2: dirs)
    //    {
    //    //if(s1!=s2)
    //    {
    //     for(auto b: rvals)
    //     {
    //               op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {b.first, b.second}, {L,L})};

	
    //     auto [fac, vec] =get_normal_form(v0);
    //     add_state(states, vec, map_sec);
    //     }
    //    }}
    //  }
//       int SS=0;

      for(int i=start; i<=r; i++)
      {
     
          for(int j=start; j<=r; j++)
      {
        for(auto s1: dirs){
	     
       for(auto s2: dirs)
       {
   
        if(i!=0 or j!=0)
        {
          int ind1=(L+i)%L;
          int ind2=(L+j)%L;
        op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {ind1, ind2}, {L,L})};

	
        auto [fac, vec] =get_normal_form(v0);
   
   
    
         if(use_symm)
        {
add_state_with_symmetries(states, vec, map_sec, L);
        }
        else{
          add_state(states, vec, map_sec);
        }
        
      

	     
        }
		 }
		   
	    }
     }
      }
}

void get_order_three_monomials(basis_structure& states, std::map<std::pair<int,int>, int>& map_sec, int L, bool use_symm)
{

std::vector<std::string> dirs={"x","y","z"};
       	   for(auto s1: dirs){
	      for(auto s2: dirs){
           for(auto s3: dirs){
	     
       {
        op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {0,1},{L,L}),spin_op(s3, {1,1}, {L,L})};
        auto [fac, vec] =get_normal_form(v0);
	       
       {
      
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
		 }
        }
      {
        op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {0,1}, {L,L}),spin_op(s3, {L-1,1}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	    
       {
     
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
       }
		 }
     {
          op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {1,0}, {L,L}),spin_op(s3, {1,1}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	     
       {
      
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
       }
		 }
        {
          op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {L-1,0}, {L,L}),spin_op(s3, {L-1,1}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	       
       {
       
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
        }
		 }
      {
          op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {1,0}, {L,L}),spin_op(s3, {2,0}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	       
       {
      
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
        }
		 }
        {
          op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {0,1}, {L,L}),spin_op(s3, {0,2}, {L,L})};
	     auto [fac, vec] =get_normal_form(v0);
	      
       {
      
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
        }
		 }
		     
     
	    }
        }
     }
}
  void get_order_four_monomials(basis_structure& states, std::map<std::pair<int,int>, int>& map_sec, int L, bool use_symm)
{
std::vector<std::string> dirs={"x","y","z"};
  for(auto s1: dirs){
	      for(auto s2: dirs){
           for(auto s3: dirs){
	     
       {
            for(auto s4: dirs){
	     
       {
        op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {1,0}, {L,L}),spin_op(s3, {0,1}, {L,L}),spin_op(s4, {1,1}, {L,L})};
        auto [fac, vec] =get_normal_form(v0);
	 
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
		 }
          {
        op_vec v0={spin_op(s1, {0,0}, {L,L}),spin_op(s2, {1,0}, {L,L}),spin_op(s3, {2,0}, {L,L}),spin_op(s4, {3,0}, {L,L})};
        auto [fac, vec] =get_normal_form(v0);
	 
         if(use_symm)
        {
add_state_with_symmetries(states, v0, map_sec, L);
        }
        else{
          add_state(states, v0, map_sec);
        }
		 }
            } }}}}
}