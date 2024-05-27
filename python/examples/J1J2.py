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
from generate_monomials import prepare_nfs_xxz_T_symmetrie,prepare_nfs_xxz_T_symmetrie_2
import momentum_funcs as mfuncs
def get_J1J2(obj, J1,J2):
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
            Is[0][0,co.ind]=J1/4 # plus one because it is in block [1,d,d
            Is[0][co.ind, 0]=J1/4
    for s in sym:          
        o_1=spin_class.spin_op(s,min(0,2))
        o_2=spin_class.spin_op(s,max(0,2))
        otot=spin_class.op_list([o_1,o_2])
        coeff,nf=spin_class.normal_form(otot, class_name=spin_class.spin_op, **{})
        if nf in obj.d_vector_map.keys():
            #print(nf.sym)
            co=obj.d_vector_map[nf]
            #print(nf.sym, co.get())
            Is[0][0,co.ind]=J2/4 # plus one because it is in block [1,d,d
            Is[0][co.ind, 0]=J2/4
    return Is[0]/2
L=10
PB=True
J=1
Delta=1
h=0
ops=prepare_nfs_xxz_T_symmetrie(L)
P = picos.Problem()
# for o in ops0:
#     print(o.sym)
# for o in ops1:
#     print(o.sym)

# for o in ops0:
#     print(o.sym)
# for o in ops1:
#     print(o.sym)

kwargs={}
basis=mfuncs.momentum_basis(ops, L, P, spin_class.spin_op,**{})
Is=get_J1J2( basis.sectors[(1,1)], J1=1, J2=1)
#print(Is)

Is_picos={}
#print(P)

Is_picos=picos.Constant("H{0}".format(0), Is)
m=4
# rho_1=general_rho_TI([0],basis)
# rho_2=general_rho_TI([0,1],basis)
# rho_3=general_rho_TI([0,1,2],basis)
# rho_list=solve_marginal_problem(P,  m)
# P.add_constraint(rho_list[1]>>0)
# #P.add_constraint(rho_list[2]>>0)
# #P.add_constraint(rho_list[3]>>0)
# P.add_constraint(rho_2==rho_list[2])

# P.add_constraint(rho_3==rho_list[3])

# P.add_constraint(rho_1==rho_list[1])
#general_rho(i=0,Ms=Ms,normal_forms=basis.sectors)
#rho_1=get_rho_2(0,1, Ms=Ms,normal_forms=total_nfs)

P.set_objective("min",(basis.sectors[(1,1)].blocks[0] | Is_picos))

P.solve(solver="mosek")
print((basis.sectors[(1,1)].blocks[0] | Is_picos).np/np.sqrt(L))
# # use marginal problem



# rho_list[1]=picos.SymmetricVariable("rho{0}".format(1), (2,2))

# P.add_constraint(picos.partial_trace(rho_list[2], 0,2)==rho_list[1])
# P.add_constraint(picos.partial_trace(rho_list[2], 1,2)==rho_list[1])
# H=mag_ham.XXZ_chain(J=J,Delta=Delta,h=h).reshape(4,4)
# #print(H)
# H_picos=picos.Constant("H", H)
# P.add_constraint(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L==(H_picos| rho_list[2]))
# Oop=np.kron(np.eye(2), sx)+np.kron( sx,np.eye(2))
# P.add_constraint(Oop*rho_list[2]-rho_list[2]*Oop==np.zeros(shape=(4,4)))
# Oop=np.kron(np.eye(2), sz)+np.kron( sz,np.eye(2))
# P.add_constraint(Oop*rho_list[2]-rho_list[2]*Oop==np.zeros(shape=(4,4)))
# Oop=np.kron(np.eye(2), sy)+np.kron( sy,np.eye(2))
# P.add_constraint(Oop*rho_list[2]-rho_list[2]*Oop==np.zeros(shape=(4,4)))

# #P.set_objective("min",(H_picos| rho_list[2]))

# #P.set_objective("min",(H_picos| rho_list[2]))
# P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L)
# P.solve(solver="mosek")
#     #print([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()])
# print(sum([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()])) 
# print((H_picos| rho_list[2])) 