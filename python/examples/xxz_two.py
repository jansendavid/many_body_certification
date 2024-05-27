#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 13:14:29 2024

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
def prepare_nfs(L):
    ops0=[]
    ops1=[]
    ops2=[]
    ops3=[]
    # ops contains spin list and list of symmetry sector
    for i in range(L):
 
        s="z"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        #print(op.sym)
        ops0+=[(op, (1,1))]
        s="x"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops0+=[(op, (1,1))]
        s="y"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops0+=[(op, (1,1))]
        
    # for i in range(L):
    #     mn=min(i, (i+1)%L)
    #     mx=max(i, (i+1)%L)
    #     s="z"
    #     op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
    #     #print(op.sym)
    #     ops0+=[(op, (1,1))]
    #     s="x"
    #     op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
    #     ops0+=[(op, (1,1))]
    #     s="y"
    #     op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
    #     ops0+=[(op, (1,1))]

        
    vecs0=npa_funcs.get_normal_form_vector(ops0, (1,1), index=1)

    
    total_vecs={(1,1): vecs0}
    total_nfs=npa_funcs.get_NF_matrix_elements(total_vecs)
    return total_vecs,total_nfs

L=10
PB=True
J=1
Delta=1
h=0
total_vecs,total_nfs=prepare_nfs(L)
P = picos.Problem()
Ms=npa_funcs.get_mom_matrix_M(P=P, vector=total_vecs,normal_forms=total_nfs,)
shapes={}
for k in Ms.keys():
    shapes[k]=Ms[k].shape

Is=npa_hams.get_xxz_npa(normal_forms=total_nfs, J=J, Delta=Delta,shapes=shapes,L=L, PB=PB)
Is_picos={}
for k in Is.keys():
    Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
# use marginal problem
m=3
rho_list=solve_marginal_problem(P,  m)
rho_1=get_rho_1(i=0,Ms=Ms,normal_forms=total_nfs)
rho_2=get_rho_2(0,1, Ms=Ms,normal_forms=total_nfs)
rho_list[1]=picos.SymmetricVariable("rho{0}".format(1), (2,2))
P.add_constraint(rho_list[1]>>0)
P.add_constraint(rho_2==rho_list[2])
P.add_constraint(rho_1==rho_list[1])
P.add_constraint(picos.partial_trace(rho_list[2], 0,2)==rho_list[1])
P.add_constraint(picos.partial_trace(rho_list[2], 1,2)==rho_list[1])
H=mag_ham.XXZ_chain(J=J,Delta=Delta,h=h).reshape(4,4)
#print(H)
H_picos=picos.Constant("H", H)
#P.add_constraint(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L==(H_picos| rho_list[2]))


#P.set_objective("min",(H_picos| rho_list[2]))
print(rho_list[2])
P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
P.solve(solver="mosek")
print("trace", np.trace(( rho_list[2])))
print("is herm",np.allclose(rho_list[2].np ,np.transpose(np.conj(rho_list[2].np) )))
#print("trace", ( rho_list[2]).np)
    #print([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()])
print(sum([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()])) 
print((H_picos| rho_list[2])) 