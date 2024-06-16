#pragma once
void test_sorting()
{
  electron_op c0("up", {0}, true);
  // std::cout<< c0<<std::endl;

  //  electron_op c1("dn", {0,3}, true);
  //std::cout<< c1<<std::endl;
  electron_op c1("up", {2}, false);
 
  electron_op c2("up", {1}, false);
 
  electron_op c3("dn", {1}, true);

    electron_op c4("up", {0}, false);

      electron_op c5("up", {2}, true);
      std::vector<electron_op> arr={c5, c1,c0,  c3, c2, c4};

      ///      std::vector<electron_op> arr_sorted={c3,};

      
  for(auto& a : arr)
    {
      std::cout<<a<<std::endl;

    }
    int coeff=bubbleSort( arr,arr.size());
    std::cout<< "sorted"<<std::endl;
    for(auto& a : arr)
    {
            std::cout<<a<<std::endl;

    }
    
    return;
}

void test_commutator()
{

//    electron_op c0("up", {0}, true);

//     electron_op c1("up", {1}, false);
//     electron_op c2("up", {1}, true);

//       electron_op c3("dn", {2}, true);
//       std::vector<electron_op> arr={c0, c1, c2, c3};
//     int coeff=bubbleSort( arr,arr.size());
// for(auto& a : arr)
//     {
//       std::cout<<a<<std::endl;

//     }
//   std::cout<< "end"<<std::endl;

//    auto new_arr= apply_commutator(arr, 2);
//    for(auto& a : arr)
//     {
//       std::cout<<a<<std::endl;

//     }
//   std::cout<< "end"<<std::endl;

//   for(auto& a : new_arr)
//     {
//       std::cout<<a<<std::endl;

//     }
//   std::cout<< "end"<<std::endl;


//    electron_op c0("up", {1}, true);

//     electron_op c1("up", {0}, false);
//     electron_op c2("up", {0}, true);

//       electron_op c3("dn", {2}, true);
//       std::vector<electron_op> arr={c0, c1, c2, c3};
//     int coeff=bubbleSort( arr,arr.size());
// for(auto& a : arr)
//     {
//       std::cout<<a<<std::endl;

//     }
//   std::cout<< "end"<<std::endl;

//    auto new_arr= apply_commutator(arr, 1);
//    for(auto& a : arr)
//     {
//       std::cout<<a<<std::endl;

//     }
//   std::cout<< "end"<<std::endl;

//   for(auto& a : new_arr)
//     {
//       std::cout<<a<<std::endl;

//     }
//   std::cout<< "end"<<std::endl;
   electron_op c0("up", {0}, true);

    electron_op c1("up", {1}, false);
    electron_op c2("up", {1}, true);

      electron_op c3("dn", {2}, true);
      std::vector<electron_op> arr={c0, c2, c1, c3};
    int coeff=bubbleSort( arr,arr.size());
for(auto& a : arr)
    {
      std::cout<<a<<std::endl;

    }
  std::cout<< "end"<<std::endl;

   auto new_arr= apply_commutator(arr, 2);
   for(auto& a : arr)
    {
      std::cout<<a<<std::endl;

    }
  std::cout<< "end"<<std::endl;

  for(auto& a : new_arr)
    {
      std::cout<<a<<std::endl;

    }
  std::cout<< "end"<<std::endl;

  return;
}
