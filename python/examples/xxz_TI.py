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
from generate_monomials import prepare_nfs_xxz_T_symmetrie
import momentum_funcs as mfuncs
def get_xxz(obj, J):
    Is={}
    Is[0]=np.zeros(obj.blocks[0].shape)
    sym=["x", "y", "z"]
    #sym=["z"]
    for s in sym:          
        o_1=spin_class.spin_op(s,min(0,1))
        o_2=spin_class.spin_op(s,max(0,1))
        otot=spin_class.op_list([o_1,o_2])
        coeff,nf=spin_class.normal_form(otot, class_name=spin_class.spin_op, **{})
        if nf in obj.d_vector_map.keys():
            #print(nf.sym)
            co=obj.d_vector_map[nf]
            #print(nf.sym, co.get())
            Is[0][0,co.ind]=J/4 # plus one because it is in block [1,d,d
            Is[0][co.ind, 0]=J/4
    return Is[0]/2
L=10
PB=True
J=1
Delta=1
h=0
operators=prepare_nfs_xxz_T_symmetrie(L)
P = picos.Problem()
# for o in ops0:
#     print(o.sym)
# for o in ops1:
#     print(o.sym)

basis=mfuncs.momentum_basis( operators, L, P, spin_class.spin_op,**{})
Is=get_xxz( basis.sectors[(1,1)], J=1)

#print(Is)

Is_picos={}
#print(P)
Is_picos=picos.Constant("H{0}".format(0), Is)

P.set_objective("min",(basis.sectors[(1,1)].blocks[0] | Is_picos))
print("start solving")
P.solve(solver="mosek")
print((basis.sectors[(1,1)].blocks[0] | Is_picos).np/np.sqrt(L))