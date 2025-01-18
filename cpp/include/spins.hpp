#pragma once
#include"vector"
#include"string"
#include <iostream>
using cpx =std::complex<double>;

class spin_op{
public:
  std::string dir_;
  std::vector<int> site_;

  std::vector<int> offset_;
  bool unit=false;
  spin_op(std::string dir, std::vector<int> site, std::vector<int> offset):dir_(dir), site_(site), offset_(offset){
    if(site.size()<offset.size()){std::cout<< "Need sufficient offsets"<<std::endl;}};
  spin_op(){ unit=true;};
  int const pos()const {
    
    int position=0;
    for(int i=site_.size()-1; i>-1; i--)
      {
       // std::cout<< "offset "<<offset_[i]<< " and factor "<<offset_.size()-1<<std::endl;
	    position+=site_[i]*std::pow(offset_[i],i);
      }
     
      return position;

  }
    spin_op get_mirror()
  {
    // assumes LxL lattice and order ....y,x
    auto new_sites=site_;

    new_sites[new_sites.size()-1]=site_[site_.size()-2];
    new_sites[new_sites.size()-2]=site_[site_.size()-1];
    
    return spin_op(dir_,new_sites, offset_ );
  }
     spin_op get_flipped_layer()
  {
    // assumes LxL lattice and order ....y,x
    assert(site_.size()==3);
    auto new_sites=site_;

    new_sites[0]=(site_[0]+1)%2;
 
    return spin_op(dir_,new_sites, offset_ );
  }
 spin_op get_translated(int j, int L)
  {
    auto new_sites=site_;
    new_sites[new_sites.size()-1]=(new_sites[new_sites.size()-1]+j)%L;
    
    return spin_op(dir_,new_sites, offset_ );
  }

   spin_op get_translated_y(int j, int L)
  {
    auto new_sites=site_;
    
    new_sites[new_sites.size()-2]=(new_sites[new_sites.size()-2]+j)%L;
    
    return spin_op(dir_,new_sites, offset_ );
  }
  
  std::string expression() const
  {
    std::string exp="s";

    exp+="_[";
    
    exp+=dir_+",(";
    for(int i=0; i< site_.size()-1;i++)
    {
      exp+=std::to_string(site_[i])+",";
    }
      exp+=std::to_string(site_[site_.size()-1]);
  exp+=")]";
  return exp;
}
   bool operator<(const spin_op& obj) const
    { 

          
	    return this->expression()<obj.expression();
	  
    } 
   bool operator>(const spin_op& obj) const
    { 
	    return this->expression()>obj.expression();        

    } 
bool operator == (const spin_op& obj)
{
   
  return (expression() == obj.expression());
}
  
    friend bool operator== (const spin_op& c1, const spin_op& c2);
    friend bool operator!= (const spin_op& c1, const spin_op& c2);
};
bool operator== (const spin_op& c1, const spin_op& c2)
{
    return (c1.expression() == c2.expression());
}

bool operator!= (const spin_op& c1, const spin_op& c2)
{
    return (c1.expression() != c2.expression());
}

std::ostream& operator<<(std::ostream& os, const spin_op& op)
{
  os << op.expression();
    return os;
}


/////////////////////////////////////////////////////////////////////////
using op_vec = std::vector<spin_op>;
using basis_structure =std::map<int, std::vector<op_vec>>;
/////////////////////////////////////////////////////////////////
std::string print_op(const op_vec& oper)
{
  std::string s="";
  if(oper.size()<1)
    {s+="1";}
  else{
  for(auto& O : oper)
    {
      s+=O.expression();
    }}
  return s;
}

op_vec dagger_operator(op_vec oper)
{
 double coeff=1;
  op_vec new_op;
  reverse(oper.begin(),oper.end()); 
 return oper;

}
template<typename T>
T apply_commutator(T& arr, int i)
{
  // if empty vector, then nothing happens (already c^dagc)
  // else change arr and return with deleted element
  
  T arr_copy;
    copy(arr.begin(), arr.begin()+i, back_inserter(arr_copy));
 
  copy(arr.begin()+i+2, arr.end(), back_inserter(arr_copy));   
      std::swap(arr[i], arr[i + 1]);
      return arr_copy;
}


std::pair<cpx,std::string>get_dir(std::string a, std::string b)
{
  std::pair<cpx, std::string> p1(cpx(0,0), std::string("0"));
  if(a=="x" and b=="y")
    {return {cpx(0,1), std::string("z")};}
  else if(a=="y" and b=="x")
    {return {cpx(0,-1), std::string("z")};}
  else if(a=="z" and b=="x")
    {return {cpx(0,1), std::string("y")};}
  else if(a=="x" and b=="z")
    {return {cpx(0,-1), std::string("y")};}
  else if(a=="y" and b=="z")
    {return {cpx(0,1), std::string("x")};}
   else if(a=="z" and b=="y")
     {return {cpx(0,-1), std::string("x")};}
   else{
    std::cout<< "error commutators "<<std::endl;
    return {cpx(0,0), std::string("0")};
      }
}
std::pair<cpx,op_vec> run_loop(op_vec op)
{
  cpx pref=1;
  std::map<int, spin_op> copied;
    for(int i=0; i<op.size(); i++)
      {
	copied.insert({i, op[i]});
      }
  for(int i=0; i<op.size()-1; i++)
    {
      auto it = copied.find(i);

    if (it != copied.end()) {
      auto o1=copied[i];
      auto o2=copied[i+1];
      if(o1.pos()!=o2.pos())
	{
	  continue;
	}
      else if(o1.dir_==o2.dir_){copied.erase(i);copied.erase(i+1);}
      else{
        auto [new_coeff, new_op]=get_dir(o1.dir_,o2.dir_);
	copied.erase(i);
	copied.erase(i+1);
	copied.insert({i, spin_op(new_op,o1.site_, o1.offset_)});
	pref*=new_coeff;

      }
      }

    }
  
    op_vec final_vec;
    for(auto it=copied.begin(); it!=copied.end(); ++it)
      {
	final_vec.push_back(it->second);
      }
    
    return {pref, final_vec};
    
 }
  
std::pair<cpx,op_vec> get_normal_form(op_vec op)
{
  if(op.size()<1){std::cout<< "error: unit sent to normal form"<<std::endl;}

  cpx pref(1.,0);
  bool change=false;
  auto new_list=op;  
  sort(new_list.begin(), new_list.end(),[](spin_op& o1, spin_op& o2){return o1.pos()<o2.pos();} );
  for(auto a: new_list)
  {
    //std::cout<<a.pos()<<std::endl;
  }
  
  while(not change)
    {

      auto [new_coeff, new_op]=run_loop(new_list);
      pref*=new_coeff;
      if(new_op.size()==0 or new_op==new_list)
	{
	  new_list=new_op;
	   change=true;}
      else{
	new_list=new_op;
       }
      
    }
  return std::pair<cpx,op_vec>{pref, new_list};
 }


// double sort(op_vec& op, int n)
// {
//   sort(op.begin(), op.end(),[](spin_op o1, spin_op o2){return o1.pos()<o2.pos();} );
//  return 1;
// }
   std::vector<std::pair<cpx,op_vec>> generate_all_terms(op_vec op, bool& unit_found)
 {
      std::vector<std::pair<cpx,op_vec>> terms;
      // auto [fac, vec] get_normal_form(op)
     // terms.push_back({fac,vec});
     // unit_found=false;
     return terms;
 }
// std::pair<bool, double> check_constant(op_vec op, bool is_con)
// {
//   bool zero=false;
//   for(int i=0; i<op.size()-1; i++)
//     {
//       if(op[i]==op[i+1])
// 	{return true;}
//     }
    
//   return zero;


// }
//   std::vector<std::pair<double,op_vec>> terms;
  
//     terms.push_back({1.,op});

//   // // bool ended=false;
//   // // auto it = begin (terms);
//   //  int ind=0;
//   //  int size=terms.size();
  
//   //  while(ind<size)
//   //    {
//   //      auto n=check_one_term(terms[ind].second,terms[ind].first, terms, unit_found);
//   //  terms.insert(terms.end(), n.begin(), n.end());
  
//   //  ind+=1;

//   //  size=terms.size();
    
//   //    }
//   //  // for(int i=0; i<terms;)
//   //   // for(auto a: terms)
//   //   //    {
//   //   // 	 std::cout<< "coeff "<< a.first<< " term "<< print_op(a.second)<<std::endl;
//   //   //    }

//   return terms;

// }


