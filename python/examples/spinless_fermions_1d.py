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
from generate_monomials import prepare_nfs_spinlessfermions_T_symmetrie
def get_sign_ferm(op):
    coeff1=1
    for o in op.ops:
        if o.cor=="x":
            coeff1*=-1
        elif o.cor=="y":
            coeff1*=-1
    return (coeff1)

        



        
    return Is
def get_spinless(obj, J):
    Is={}
    Is[0]=np.zeros(obj.blocks[0].shape)
    for s in ["x", "y"]:
        l=[spin_class.spin_op(s,0)]
        l+=[spin_class.spin_op(s,1)]
        op=spin_class.op_list(l)

        coeff,nf=spin_class.normal_form(op, class_name=spin_class.spin_op, **{})
        print("sarch",nf.sym)
        if nf in obj.d_vector_map.keys():
            #print(nf.sym)
            print("found", nf)
            co=obj.d_vector_map[nf]
            #print(nf.sym, co.get())
            Is[0][0,co.ind]=J*(L-1) # plus one because it is in block [1,d,d
            Is[0][co.ind, 0]=J*(L-1)
    # print("done, x")           
    for s in ["x","y"]:
        l=[spin_class.spin_op(s,0)]
        for j in range(1,L-1):
            l+=[spin_class.spin_op("z",j)]
        l+=[spin_class.spin_op(s,L-1)]
        op=spin_class.op_list(l)
        coeff,nf=spin_class.normal_form(op, class_name=spin_class.spin_op, **{})
        sign=get_sign_ferm(nf)
        print("search", nf.sym)
        if nf in obj.d_vector_map.keys():
            print("found", coeff, nf.sym)
            #print(nf.sym)
            Is[0][0,co.ind]=J # plus one because it is in block [1,d,d
            Is[0][co.ind, 0]=J


    print("done")    
    return Is[0]/2

def get_oarticle_nr(obj, J):
    Is={}
    Is[0]=np.zeros(obj.blocks[0].shape)
    print("particle")
    for s in ["z"]:
        l=[spin_class.spin_op(s,0)]
        op=spin_class.op_list(l)

        coeff,nf=spin_class.normal_form(op, class_name=spin_class.spin_op, **{})
        if nf in obj.d_vector_map.keys():
            #print(nf.sym)
            print("found", nf.sym)
            co=obj.d_vector_map[nf]
            #print(nf.sym, co.get())
            Is[0][0,co.ind]=1. # plus one because it is in block [1,d,d
            Is[0][co.ind, 0]=1.


    
    return Is[0]/2


L=4
PB=True
J=1
Delta=1
h=0
operators=prepare_nfs_spinlessfermions_T_symmetrie(L)

P = picos.Problem()
basis=mfuncs.momentum_basis( operators, L, P, spin_class.spin_op,**{})
Is=get_spinless(basis.sectors[(1)], 1)
n=get_oarticle_nr(basis.sectors[(1)], 1)
Is_picos={}
#print(P)
Is_picos=picos.Constant("H{0}".format(0), Is)

n_picos={}
#print(P)
n_picos=picos.Constant("n{0}".format(0), n)

#P.add_constraint(L*(0.5*((basis.sectors[(1)].blocks[0] | n_picos)/np.sqrt(L)-1.))==2)

print( Is_picos)
P.set_objective("min",(basis.sectors[(1)].blocks[0] | Is_picos))


P.solve(solver="mosek")
print((basis.sectors[(1)].blocks[0] | Is_picos).np/np.sqrt(L))

print(basis.sectors[(1)].blocks[0].np[0,:]/np.sqrt(L))
print(0.5*(basis.sectors[(1)].blocks[0].np[0,1]/np.sqrt(L)-1.))


# print(basis.sectors[(1)].blocks[0].np[0,0])


