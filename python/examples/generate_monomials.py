#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 14:43:37 2024

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
# def prepare_nfs_tfisign_symmetrie(L):
#     ops1=[]
#     ops2=[]
#     # ops contains spin list and list of symmetry sector
#     for i in range(L):
#         for s in sym:
#             op=spin_class.op_list([spin_class.spin_op(s,i)])
#             if s=="x":
#                 ops1+=[(op, 1)]
#             else:
#                 ops2+=[(op,-1)]


#     vecs1=npa_funcs.get_normal_form_vector(ops1, (1), index=1)
#     vecs2=npa_funcs.get_normal_form_vector(ops2, (-1), index=0)
#     #for n in vecs2.keys():
#     #    print(n.sym, vecs2[n].get())
#     total_vecs={ (+1): vecs1,(-1):vecs2}#

#     total_nfs=npa_funcs.get_NF_matrix_elements(total_vecs)

#     return total_vecs,total_nfs
# def prepare_nfs_xxzsign_symmetrie(L):
#     ops0=[]
#     ops1=[]
#     ops2=[]
#     ops3=[]
#     # ops contains spin list and list of symmetry sector
#     ops0+=[(spin_class.op_list([]),(1,1))] # empty list gives unity
#     for i in range(L):
#         mn=min(i, (i+1)%L)
#         mx=max(i, (i+1)%L)
#         s="z"
#         op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
#         #print(op.sym)
#         ops0+=[(op, (1,1))]
#         s="x"
#         op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
#         ops0+=[(op, (1,1))]
#         s="y"
#         op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
#         ops0+=[(op, (1,1))]
#     for i in range(L):
#         s="z"
#         op=spin_class.op_list([spin_class.spin_op(s,i)])
#         ops1+=[(op, (1,-1))]
#         s="x"
#         op=spin_class.op_list([spin_class.spin_op(s,i)])
#         ops2+=[(op, (-1,1))]
#         s="y"
#         op=spin_class.op_list([spin_class.spin_op(s,i)])
#         ops3+=[(op, (-1,-1))]
        
#     for i in range(L):
#         mn=min(i,(i+1)%L)
#         mx=max(i,(i+1)%L)
#         s1="x"
#         s2="y"
#         op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
#         ops1+=[(op, (1,-1))]
        
#         op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
#         ops1+=[(op, (1,-1))]
#         s1="y"
#         s2="z"
#         op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
#         ops2+=[(op, (-1,1))]
#         op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
#         ops2+=[(op, (-1,1))]
#         s1="x"
#         s2="z"
#         op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
#         ops3+=[(op, (-1,-1))]
#         op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
#         ops3+=[(op, (-1,-1))]
        
#     vecs0=npa_funcs.get_normal_form_vector(ops0, (1,1))
#     vecs1=npa_funcs.get_normal_form_vector(ops1, (1,-1))
#     vecs2=npa_funcs.get_normal_form_vector(ops2, (-1,1))
#     vecs3=npa_funcs.get_normal_form_vector(ops3, (-1,-1))
    
#     total_vecs={(1,1): vecs0,(1,-1): vecs1, (-1,1):vecs2, (-1,-1): vecs3}
#     total_nfs=npa_funcs.get_NF_matrix_elements(total_vecs)
#     return total_vecs,total_nfs
def get_sign(op):
    coeff1=1
    coeff2=1
    for o in op.ops:
        if o.cor=="x":
            coeff1*=-1
        if o.cor=="y":
            coeff1*=-1 
            coeff2*=-1
        if o.cor=="z":
            coeff2*=-1
    return (coeff1, coeff2)
def prepare_nfs_xxz_T_symmetrie(L):
    ops0=[] # (1,1)
    ops1=[] #(1,-1)
    ops2=[] #(-1,1)
    ops3=[] #(-1,-1)
    ops={(1,1):[],(1,-1):[] , (-1,1):[], (-1,-1):[]}
    # ops contains spin list and list of symmetry sector
    r=int(L/2)

    for i in range(1,r+1):
        for s1 in ["z" , "x", "y"]:
            for s2 in ["z" , "x", "y"]:
                mn=min(0, (i)%L)
                mx=max(0, (i)%L)
    #             s="z"
                op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
                sign=get_sign(op)
                #print(op.sym, sign)
                ops[sign]+=[op]
 
                        
    for i in range(1):
        for s in ["z" , "x", "y"]:
            op=spin_class.op_list([spin_class.spin_op(s,i)])
            
            sign=get_sign(op)
            print(op.sym, sign)
            ops[sign]+=[op]

    return ops


def prepare_nfs_xxz_T_symmetrie_2(L):
    ops0=[] # (1,1)
    ops1=[] #(1,-1)
    ops2=[] #(-1,1)
    ops3=[] #(-1,-1)
    ops={(1,1):[],(1,-1):[] , (-1,1):[], (-1,-1):[]}
    # ops contains spin list and list of symmetry sector
    r=3#int(L/2)

    for i in range(1,r+1):
        for s1 in ["z" , "x", "y"]:
            for s2 in ["z" , "x", "y"]:
                mn=min(0, (i)%L)
                mx=max(0, (i)%L)
    #             s="z"
                op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
                sign=get_sign(op)
                #print(op.sym, sign)
                coeff,nf=spin_class.normal_form(op, spin_class.spin_op,**{})
                ops[sign]+=[nf]
    for i in range(2,r+1):
        for s1 in ["z" , "x", "y"]:
            for s2 in ["z" , "x", "y"]:
                mn=min(0, (i)%L)
                mx=max(0, (i)%L)
    #             s="z"
                op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s1,mn+1),spin_class.spin_op(s2,mx+1),spin_class.spin_op(s2,mx)])
                sign=get_sign(op)
                #print(op.sym, sign)
                coeff,nf=spin_class.normal_form(op, spin_class.spin_op,**{})
                ops[sign]+=[nf]
 
                        
    for i in range(1):
        for s in ["z" , "x", "y"]:
            op=spin_class.op_list([spin_class.spin_op(s,i)])
            
            sign=get_sign(op)
            print(op.sym, sign)
            ops[sign]+=[op]

    return ops




def get_sign_ferm(op):
    coeff1=1
    for o in op.ops:
        if o.cor=="x":
            coeff1*=-1
        elif o.cor=="y":
            coeff1*=-1
    return (coeff1)
def prepare_nfs_spinlessfermions_T_symmetrie(L):
    ops0=[] # (1,1)
    ops1=[] #(1,-1)
    ops2=[] #(-1,1)
    ops3=[] #(-1,-1)
    ops={1:[],-1:[]}
    for i in range(1):
        for s in ["z" , "x", "y"]:
            op=spin_class.op_list([spin_class.spin_op(s,i)])
            
            sign=get_sign_ferm(op)
            print(op.sym, sign)
            ops[sign]+=[op]
    for i in range(1):
        mn=min(i,(i+1)%L)
        mx=max(i,(i+1)%L)
        s1="x"
        s2="y"
        op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
        sign=get_sign_ferm(op)
        ops[sign]+=[op]
        
        op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
        sign=get_sign_ferm(op)
        ops[sign]+=[op]
        s1="y"
        s2="z"
        op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
        sign=get_sign_ferm(op)
        ops[sign]+=[op]
        op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
        sign=get_sign_ferm(op)
        ops[sign]+=[op]
        s1="x"
        s2="z"
        op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
        sign=get_sign_ferm(op)
        ops[sign]+=[op]
        op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
        sign=get_sign_ferm(op)
        ops[sign]+=[op]
    # terms to generat electron hopping
    # print("ere")
    # for i in range(1):
    #     for s1 in ["x", "y", "z"]:
    #         for s2 in ["x", "y", "z"]:
    #             op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i+1)])            
    #             sign=get_sign_ferm(op)
    #             print(op.sym, sign)
    #             ops[sign]+=[op]
    # for i in range(1):
    #     for s1 in ["x", "y"]:
    #         l=[spin_class.spin_op(s1,0)]
    #         for j in range(1,L-2):
    #             l+=[spin_class.spin_op("z",j)]
    #         l+=[spin_class.spin_op(s1,L-1)]
    #         op=spin_class.op_list(l)
    #         sign=get_sign_ferm(op)
    #         print(op.sym, sign)
    #         ops[sign]+=[op]
 
            
    for i in range(1):
        for s1 in ["x","y", "z"]:
            l=[spin_class.spin_op(s1,0)]
            l+=[spin_class.spin_op(s1,i+1)]
            op=spin_class.op_list(l)
            sign=get_sign_ferm(op)
            print(op.sym, sign)
            ops[sign]+=[op]               
    # get_sign_ferm(op)
    return ops
    
    
    