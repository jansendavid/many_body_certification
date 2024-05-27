#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 16:03:39 2024

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
def get_mag_npa(normal_forms,shapes):

    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k], dtype=complex)
    #end=L-1
    #if PB:
    #    end=L
    sym=["z"]
    for i in np.arange(0,1,1):
        
        for s in sym:          
            o_1=spin_class.spin_op(s,i)
            #o_2=spin_class.spin_op(s,max(i, (i+1)%L))
            otot=spin_class.op_list([o_1])
            add_term=True
            for symmetrie in normal_forms.keys():
                if otot in normal_forms[symmetrie].keys() and add_term:
                    matel=normal_forms[symmetrie][otot]
                    if s=="z":
                        Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.
                        Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.
                    break
    for k in Is.keys():
        #print("np", k, Is[k])
        Is[k]=0.5*Is[k]
    return Is

def get_Lb_npa(normal_forms,shapes, gx, gy):

    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k], dtype=complex)
    

    otot=spin_class.op_list([spin_class.spin_op("x",0), spin_class.spin_op("z",0), spin_class.spin_op("x",0)])
    coeff_tot,nf_tot=spin_class.normal_form(otot)

    add_term=True
    for symmetrie in normal_forms.keys():
        #if otot in normal_forms[symmetrie].keys() and add_term:
        #for n in normal_forms[symmetrie].keys():
        #    print(n.sym)
        matel=normal_forms[symmetrie][nf_tot]
        
    #         if s=="z":
        Is[symmetrie][matel.ind1,matel.ind2]+=gx*np.conj(matel.prefac)*coeff_tot*1.
        Is[symmetrie][matel.ind2,matel.ind1]+=gx*matel.prefac*coeff_tot*1.
        break

    otot=spin_class.op_list([spin_class.spin_op("y",0), spin_class.spin_op("z",0), spin_class.spin_op("y",0)])
    coeff_tot,nf_tot=spin_class.normal_form(otot)

    add_term=True
    for symmetrie in normal_forms.keys():
        #if otot in normal_forms[symmetrie].keys() and add_term:
        #for n in normal_forms[symmetrie].keys():
        #    print(n.sym)
        matel=normal_forms[symmetrie][nf_tot]
        
    #         if s=="z":
        Is[symmetrie][matel.ind1,matel.ind2]+=gy*np.conj(matel.prefac)*coeff_tot*1.
        Is[symmetrie][matel.ind2,matel.ind1]+=gy*matel.prefac*coeff_tot*1.
        break

    otot=spin_class.op_list([spin_class.spin_op("z",0)])
    coeff_tot,nf_tot=spin_class.normal_form(otot)
 
    add_term=True
    for symmetrie in normal_forms.keys():
        #if otot in normal_forms[symmetrie].keys() and add_term:
        #for n in normal_forms[symmetrie].keys():
        #    print(n.sym)
        matel=normal_forms[symmetrie][nf_tot]
        
    #         if s=="z":
        Is[symmetrie][matel.ind1,matel.ind2]+=-(gx+gy)*np.conj(matel.prefac)*coeff_tot*1.
        Is[symmetrie][matel.ind2,matel.ind1]+=-(gx+gy)*(matel.prefac)*coeff_tot*1.
        break


    for k in Is.keys():
        #print("np", k, Is[k])
        Is[k]=0.5*Is[k]
    return Is

def generate_mons():
    ops1=[]
    #for i in range(1):
        #mn=min(i, (i+1)%L)
        # mx=max(i, (i+1)%L)
        # s="z"
        # op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        # #print(op.sym)
        # ops0+=[(op, (1,1))]
        # s="x"
        # op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        # ops0+=[(op, (1,1))]
        # s="y"
        # op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        # ops0+=[(op, (1,1))]
    ops1+=[(spin_class.op_list([]),(1,1))]
    for i in range(1):
        s="z"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[(op, (1,1))]
        s="x"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[(op, (1,1))]
        s="y"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[(op, (1,1))]
        
    for i in range(1):
        s1="z"
        s2="y"
        op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i)])
        ops1+=[(op, (1,1))]
        s1="z"
        s2="x"
        op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i)])
        #op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[(op, (1,1))]
        s1="x"
        s2="y"
        op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i)])
        #op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[(op, (1,1))]
        
    for i in range(1):
        s1="z"
        s2="y"
        s3="x"
        op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i),spin_class.spin_op(s3,i)])
        ops1+=[(op, (1,1))]
        s1="z"
        s2="x"
        s3="z"
        op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i),spin_class.spin_op(s3,i)])
        #op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[(op, (1,1))]
        s1="x"
        s2="z"
        s3="x"
        op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i),spin_class.spin_op(s3,i)])
        #op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[(op, (1,1))]        
        
        
    for o in ops1:
        print(o[0].sym)
    print("len",len(ops1))
    vecs0=npa_funcs.get_normal_form_vector(ops1, (1,1))
    total_vecs={(1,1): vecs0}
    total_nfs=npa_funcs.get_NF_matrix_elements(total_vecs)
    #print("end")

    return total_vecs,total_nfs
total_vecs,total_nfs=generate_mons()
#for v in total_nfs[(1,1)].keys():
#    print(v.sym)
# for v in total_nfs[(1,1)].keys():
#     print(v.sym, v)
print(len(total_vecs[(1,1)]))
P = picos.Problem()
Ms=npa_funcs.get_mom_matrix_M(P=P, vector=total_vecs,normal_forms=total_nfs,)
print(Ms[(1,1)].shape)
shapes={}
for k in Ms.keys():
    shapes[k]=Ms[k].shape

mag=get_mag_npa(total_nfs,shapes=shapes)
Is_picos={}
for k in mag.keys():
    Is_picos[k]=picos.Constant("H{0}".format(k), mag[k])
Lb=get_Lb_npa(total_nfs,shapes=shapes, gx=1, gy=1)
Lb_picos={}
for k in mag.keys():
    print(Lb[k])
    Lb_picos[k]=picos.Constant("H{0}".format(k), Lb[k])
# print(P)
P.add_constraint(sum([(Ms[k] | Lb_picos[k]) for k in Ms.keys()])==0)
P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
P.solve(solver="mosek")
print(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]).np)
# print(Ms[(1,1)].np)