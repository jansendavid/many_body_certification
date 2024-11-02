#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "fusion.h"
#include <bits/stdc++.h>
#include"spins.hpp"
#include<unordered_map>
#include <Eigen/Dense>

#include"test_momentum_spins_eff.hpp"
using namespace mosek::fusion;
using namespace monty;

int main()
{
  //     test_translation_2d();
  //     test_permutations();
   //  test_single_block();
   //   test_multiple_blocks();
   // test_multiple_blocks_higher_order();
     // test_multiple_blocks_2d();
  //test_multiple_blocks_higher_order_2d();
  test_multiple_blocks_higher_order_2d_rdm();
  
       // test_J1J2_2d();
  //  test_load_basis();
         return 0;
}
