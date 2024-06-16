#pragma once
#include"vector"
#include"string"
#include <iostream>

class electron_op{
public:
  std::string spin_;
  std::vector<int> site_;
  bool dagger_;
  int offset_{0};
  bool unit=false;
  electron_op(std::string spin, std::vector<int> site,bool dagger ):spin_(spin), site_(site), dagger_(dagger){};
  electron_op(){ unit=true;};
  int const pos()const {
    int position=0;
    if(site_.size()>1)
      {
	position+=site_[0]*offset_;
      }
      position+=site_[site_.size()-1];
      return position;

  }
  electron_op get_translated(int j, int L)
  {
    auto new_sites=site_;
    new_sites[new_sites.size()-1]=(new_sites[new_sites.size()-1]+j)%L;
    
    return electron_op(spin_,new_sites,dagger_ );
  }
  
  std::string expression() const
  {
    std::string exp="c";
    if(dagger_)
      {
	exp+="^dag";}
    exp+="_[";
    
    exp+=spin_+",(";
    for(int i=0; i< site_.size()-1;i++)
    {
      exp+=std::to_string(site_[i])+",";
    }
      exp+=std::to_string(site_[site_.size()-1]);
  exp+=")]";
  return exp;
}
   bool operator<(const electron_op& obj) 
    { 

            if (this->spin_ == obj.spin_)
	  {
	    return this->pos()<obj.pos();
	  }
	            return this->spin_ < obj.spin_; 
    } 
   bool operator>(const electron_op& obj) 
    { 

      if (this->spin_ == obj.spin_)
	  {
	    return this->pos()>obj.pos();
	  }
        
        
        return this->spin_ > obj.spin_; 
    } 
bool operator == (const electron_op& obj)
{
   if (spin_ == obj.spin_ && dagger_ == obj.dagger_ and (site_ == obj.site_))
      return true;
  return false;
}
};

std::ostream& operator<<(std::ostream& os, const electron_op& op)
{
  os << op.expression();
    return os;
}
class unit_op: public electron_op{
public:
  unit_op(){};
  std::string expression() const
  {
    std::string s="1";
    return s;
  }
};

template<typename T>
int bubbleSort(T& arr, int n)
{
    int i, j;
    bool swapped;
    int coeff=1;
    
    for (i = 0; i < n - 1; i++) {
        swapped = false;
        for (j = 0; j < n - i - 1; j++) {
	  
            if (arr[j] > arr[j + 1]) {

	      std::swap(arr[j], arr[j + 1]);
                swapped = true;
		coeff*=-1;
            }
        }

        // If no two elements were swapped
        // by inner loop, then break
        if (swapped == false)
            break;
    }

    return coeff;
}


using op_vec = std::vector<electron_op>;
std::string print_op(const op_vec& oper)
{
  std::string s="";
  for(auto& O : oper)
    {
      s+=O.expression();
    }
  return s;
}

op_vec dagger_operator(op_vec oper)
{
 double coeff=1;
  op_vec new_op;
  reverse(oper.begin(),oper.end());
  std::for_each(oper.begin(), oper.end(), [](electron_op &O){ O.dagger_=!O.dagger_; });
 
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
bool check_zero(op_vec op)
{
  bool zero=false;
  for(int i=0; i<op.size()-1; i++)
    {
      if(op[i]==op[i+1])
	{return true;}
    }
    
  return zero;


}

std::vector<std::pair<double,op_vec>> check_one_term(op_vec& op, double& coeff, std::vector<std::pair<double,op_vec>>& terms, bool& unit_found)
{
  auto it = op.begin();
  int index=0;
  auto end=op.size();
  int i=0;
  
  std::vector<std::pair<double,op_vec>> new_terms;
  while(i<end-1)
    {
 
       
      if(op[i].pos()==op[i+1].pos())
	{
      if(op[i].dagger_==false and op[i+1].dagger_==true)
	{
      auto m=apply_commutator(op, i);
      coeff*=-1;
      if(m.size()>0)
	{
	  bool zero=check_zero(m);
	  
	  new_terms.push_back({-1*coeff,m});
	      
	}
      else{
	unit_found=true;
      }
	}
	}
      i++;
    }
  
  return new_terms;
}


std::vector<std::pair<double,op_vec>> generate_all_terms(op_vec op, bool& unit_found)
{
  std::vector<std::pair<double,op_vec>> terms;
  
    terms.push_back({1.,op});

  // bool ended=false;
  // auto it = begin (terms);
   int ind=0;
   int size=terms.size();
  
   while(ind<size)
     {
       auto n=check_one_term(terms[ind].second,terms[ind].first, terms, unit_found);
   terms.insert(terms.end(), n.begin(), n.end());
  
   ind+=1;

   size=terms.size();
    
     }
   // for(int i=0; i<terms;)
    // for(auto a: terms)
    //    {
    // 	 std::cout<< "coeff "<< a.first<< " term "<< print_op(a.second)<<std::endl;
    //    }

  return terms;

}


