#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:28:41 2024

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
def get_sign_ferm(op):
    coeff1=1
    for o in op.ops:
        if o.cor=="x":
            coeff1*=-1
        elif o.cor=="y":
            coeff1*=-1
    return (coeff1)

def get_nr_particles_npa(normal_forms,L, PB=False):
    print("part number")
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k])
    end=L-1
    if PB:
        end=L
    for i in np.arange(0,L,1):
        s="z"
        o_1=spin_class.spin_op(s,i)
        otot=spin_class.op_list([o_1])
        for symm in normal_forms.keys():
            if otot in normal_forms[symm].keys():
                matel=normal_forms[symm][otot]
                Is[symm][matel.ind1,matel.ind2]+=1.
                Is[symm][matel.ind2,matel.ind1]+=1.
                break
    
    for k in Is.keys():
        Is[k]=0.5*Is[k]
    return Is    
    
def get_spinless_npa(normal_forms, J, U,alpha,shapes,L, PB=False):
    """
    Transverse field ising model Hamiltonian H=\sum_i J(S^{x}_i_S^{x}_{i+1}+S^{y}_i_S^{y}_{i+1})+ Delta S^{z}_i_S^{z}_{i+1}

    J=t_0/2, U=V/4    
    """
    
    print("Ham")
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k])
    end=L-1
    if PB:
        end=L
    for i in np.arange(0,L,1):
        for s in ["x","y"]:
            o_1=spin_class.spin_op(s,min(i, (i+1)%L))
            o_2=spin_class.spin_op(s,max(i, (i+1)%L))
            otot=spin_class.op_list([o_1,o_2])
            for symm in normal_forms.keys():
                if otot in normal_forms[symm].keys():
                    matel=normal_forms[symm][otot]
                    #print("SSS",symm,otot.sym, (matel.ind1,matel.ind2))
                    Is[symm][matel.ind1,matel.ind2]+=-1.*J*(L-1)/L
                    Is[symm][matel.ind2,matel.ind1]+=-1.*J*(L-1)/L
                    break

            
    #if L>2:
        
    for s in ["x","y"]:
        for i in np.arange(0,L,1):
        #print("a",s)
            l=[spin_class.spin_op(s,i)]
            for j in range(1,L-1):
                l+=[spin_class.spin_op("z",(i+j)%L)]
            l+=[spin_class.spin_op(s,(i+L-1)%L)]
            op=spin_class.op_list(l)
            coeff,nf=spin_class.normal_form(op, class_name=spin_class.spin_op, **{})
            for symm in normal_forms.keys():
                if nf in normal_forms[symm].keys():       
                    matel=normal_forms[symm][nf]
                #print("xx", symm, "op",nf.sym, (matel.ind1,matel.ind2))
                    Is[symm][matel.ind1,matel.ind2]+=-1.*J/L
                    Is[symm][matel.ind2,matel.ind1]+=-1.*J/L
                    break
            
    for i in np.arange(0,L,1):
        for s in ["z"]:
            o_1=spin_class.spin_op(s,min(i, (i+1)%L))
            o_2=spin_class.spin_op(s,max(i, (i+1)%L))
            otot=spin_class.op_list([o_1,o_2])
            for symm in normal_forms.keys():
                if otot in normal_forms[symm].keys():
                    matel=normal_forms[symm][otot]
                    #print("SSS",symm,otot.sym, (matel.ind1,matel.ind2))
                    Is[symm][matel.ind1,matel.ind2]+=1.*U
                    Is[symm][matel.ind2,matel.ind1]+=1.*U
                    break
                                  
    for i in np.arange(0,L,1):
        for s in ["z"]:
            o_1=spin_class.spin_op(s,i)
            otot=spin_class.op_list([o_1])
            for symm in normal_forms.keys():
                if otot in normal_forms[symm].keys():
                    matel=normal_forms[symm][otot]
                    #print("SSS",symm,otot.sym, (matel.ind1,matel.ind2))
                    Is[symm][matel.ind1,matel.ind2]+=-1.*alpha
                    Is[symm][matel.ind2,matel.ind1]+=-1.*alpha
                    break               
    for k in Is.keys():
        Is[k]=0.5*Is[k]
    return Is
# ops contains spin list and list of symmetry sector
L=4
Npart=2.
PB=True
J=1./2
U=1./4
C=L*U
alpha=2*U
ops={1: [], -1: []}
ops[(1)]+=[spin_class.op_list([])] # empty list gives unity
print(sym)
for i in range(L):
    for s in sym:
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        sign=get_sign_ferm(op)
        #print(op.sym)
        ops[sign]+=[op]
      #   print(s)
        for s2 in ["z","y","x"]:
            mn=min(i,(i+1)%L)
            mx=max(i,(i+1)%L)
            op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s2,mx)])
            sign=get_sign_ferm(op)
            #print("NN",op.sym, sign)
            if op not in ops[sign]:
    #         print(op.sym)
                ops[sign]+=[op]

if L>2:
   for s in ["x","y"]:
       for i in np.arange(0,L,1):
       #print("a",s)
           l=[spin_class.spin_op(s,i)]
           for j in range(1,L-1):
               l+=[spin_class.spin_op("z",(i+j)%L)]
           l+=[spin_class.spin_op(s,(i+L-1)%L)]
           op=spin_class.op_list(l)
           coeff,nf=spin_class.normal_form(op, class_name=spin_class.spin_op, **{})
           sign=get_sign_ferm(nf)
           if nf not in ops[sign]:
               #print("here",sign, op.sym )
               ops[sign]+=[nf]

            
   for s in ["x","y", "z"]:
        l=[spin_class.spin_op(s,0)]
        for j in range(1,L):
            l+=[spin_class.spin_op("z",j)]
        op=spin_class.op_list(l)
        sign=get_sign_ferm(op)
        if op not in ops[sign]:
            #print("here",sign, op.sym )
            ops[sign]+=[op]
    

    
# print("start")
# for s in ops.keys():
#     print("sy,",s)
#     for d in ops[s]:
#         print(d.sym)
# print("done")    
    
d={}
basis=npa_funcs.basis(ops=ops, class_name=spin_class.spin_op, **{})
P = picos.Problem()
Ms=basis.get_mom_matrix_M(P)
shapes={}
for k in Ms.keys():
    shapes[k]=Ms[k].shape




Is=get_spinless_npa(normal_forms=basis.total_nfs, J=J, U=U,alpha=alpha,shapes=shapes,L=L, PB=PB)
Is_picos={}
for k in Is.keys():
    #print("sector ", k)
    #print( Is[k])
    Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
    
Ns=get_nr_particles_npa(normal_forms=basis.total_nfs,L=L, PB=False)
Ns_picos={}
for k in Ns.keys():
    Ns_picos[k]=picos.Constant("N{0}".format(k), Ns[k])
    #print(Ns[k])

# enforce translation invariance

for i in range(L):
        s="z"
        o_1=spin_class.spin_op(s,i)
        otot=spin_class.op_list([o_1])
        sign=get_sign_ferm(otot)
        matel=basis.total_nfs[sign][otot]
        P.add_constraint(Ms[sign][matel.ind1,matel.ind2]==1.-2*(Npart/L))
        
    

P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
#print(P)
P.solve(solver="mosek")
#for k in Is.keys():
#   print(Is[k])
print(sum([(Ms[k] | Is_picos[k]).np for k in Ms.keys()])+C)
# # print("particle number",-sum([(Ms[k] | Ns_picos[k]) for k in Ms.keys()])*0.5+(L/2))
# # #print(L*Ms[(1)].np[:,0])
# # #print("particle number",L*(Ms[(1)])()
# # #print((Ms[1] | Ns_picos[1]).np)
# # #print(Ms[1] | Ns_picos[1])
# # #print(Ms[1] * Ns_picos[1])

