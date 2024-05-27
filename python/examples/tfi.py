#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:22:42 2024

@author: david
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:03:35 2024

@author: david
"""

import numpy as np
import sys
sys.path.insert(0,"../")
from copy import deepcopy
import spin_class as spin_class
import npa_funcs as npa_funcs
from npa_funcs import sym
import picos
import npa_hams as npa_hams
import sys
from marginal_problem import *
import marginal_problem_hamiltonians as mag_ham
from generate_monomials import prepare_nfs_tfisign_symmetrie
PB=True
J=4
h=4
L=12
total_vecs,total_nfs=prepare_nfs_tfisign_symmetrie(L)

P = picos.Problem()
Ms=npa_funcs.get_mom_matrix_M(P=P, vector=total_vecs,normal_forms=total_nfs,)
shapes={}
for k in Ms.keys():
    shapes[k]=Ms[k].shape

Is=npa_hams.get_tfi_npa(normal_forms=total_nfs, J=J, h=h,shapes=shapes,L=L, PB=PB)
Is_picos={}
for k in Is.keys():
    Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
# use marginal problem#
m=4
rho_list=solve_marginal_problem(P,  m)
rho_1=get_rho_1(i=0,Ms=Ms,normal_forms=total_nfs)
rho_2=get_rho_2(0,1, Ms=Ms,normal_forms=total_nfs)
rho_list[1]=picos.HermitianVariable("rho{0}".format(1), (2,2))
P.add_constraint(rho_list[1]>>0)
P.add_constraint(rho_2==rho_list[2])
P.add_constraint(rho_1==rho_list[1])
P.add_constraint(picos.partial_trace(rho_list[2], 0,2)==rho_list[1])
P.add_constraint(picos.partial_trace(rho_list[2], 1,2)==rho_list[1])
H=mag_ham.TFI(J=J,h=h).reshape(4,4)
#print(H)
H_picos=picos.Constant("H", H)
P.add_constraint(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L==(H_picos| rho_list[2]))


#P.set_objective("min",(H_picos| rho_list[2]))
P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
P.solve(solver="mosek")
    #print([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()])
print(sum([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()])) 
print((H_picos| rho_list[2])) 