#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:26:42 2024

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
from generate_monomials import  prepare_nfs_tfisign_symmetrie
np.random.seed(1)
def get_obs(normal_forms, i,j):
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k], dtype=complex)
    s1="z"
    s2="z"        
    o_1=spin_class.spin_op(s1,i)
    o_2=spin_class.spin_op(s2,j)
    otot=spin_class.op_list([o_1,o_2])
    add_term=True
    for symmetrie in normal_forms.keys():
        if otot in normal_forms[symmetrie].keys() and add_term:
            matel=normal_forms[symmetrie][otot]
            Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1./4
            Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1./4
            break
    for k in Is.keys():
        Is[k]=0.5*Is[k]
    return Is
    
def test_ppt(rho):
    dim=2
    rho_red=rho.reshape((dim,dim,dim,dim))
    #rho_red_ppt=ncon([rho_red], [[-3,-2,-1,-4]])
    rho_red_ppt=np.transpose(rho_red, (2,1,0,3))
    v_red,w_red=np.linalg.eigh(rho_red_ppt.reshape(dim*dim,dim*dim))
    return v_red, w_red, rho_red_ppt.reshape(dim*dim,dim*dim)
mats=[]
for i in np.arange(0,4,1):
    m=np.zeros((2,2))
    m[0,0]=1
    m[1,1]=-1
    mat_r=np.random.rand(2,2)
    mat_r=mat_r+np.transpose(mat_r)
    v,w=np.linalg.eigh(mat_r)
    M=np.transpose(w)@m@w
    mats+=[M]
W=0.5
J=1
Delta=1
PB=True
L=20

hx=np.random.uniform(low=-W/2, high=W/2,size=(L) )
hz=np.random.uniform(low=-W/2, high=W/2,size=(L) )
#print(hx)

total_vecs,total_nfs=prepare_nfs_tfisign_symmetrie(L)
P = picos.Problem()
Ms=npa_funcs.get_mom_matrix_M(P=P, vector=total_vecs,normal_forms=total_nfs,)
shapes={}
for k in Ms.keys():
    shapes[k]=Ms[k].shape

Is=npa_hams.get_mbl_npa(normal_forms=total_nfs, J=J, hx=hx, hz=hz,shapes=shapes,L=L, PB=PB)

Is_squared=npa_hams.get_mbl_squared_npa(normal_forms=total_nfs, J=J, hx=np.zeros(L), hz=np.zeros(L),shapes=shapes,L=L, PB=PB)
Is_picos={}
for k in Is.keys():
    Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
Is_picos_squared={}
for k in Is.keys():
    Is_picos_squared[k]=picos.Constant("HH{0}".format(k), Is_squared[k])
m=4
rho_list=solve_marginal_problem_noTI(P,  m)

rho_list_2=solve_marginal_problem_noTI(P,  m,name="_2")
rho_1=general_rho([0], Ms=Ms,normal_forms=total_nfs)

rho_2=general_rho([0,1], Ms=Ms,normal_forms=total_nfs)

#rho_2_2=general_rho([1,2], Ms=Ms,normal_forms=total_nfs)
#rho_3_2=general_rho([1,2,3], Ms=Ms,normal_forms=total_nfs)
#rho_4_2=general_rho([1,2,3,4], Ms=Ms,normal_forms=total_nfs)


rho_3=general_rho([0,1,2], Ms=Ms,normal_forms=total_nfs)

rho_4=general_rho([0,1,2,3], Ms=Ms,normal_forms=total_nfs)

rho_list[1]=picos.HermitianVariable("rho{0}".format(1), (2,2))
P.add_constraint(rho_list[1]>>0)
P.add_constraint(rho_4==rho_list[4])
P.add_constraint(rho_3==rho_list[3])
P.add_constraint(rho_2==rho_list[2])
P.add_constraint(rho_1==rho_list[1])

#P.add_constraint(rho_2_2==rho_list_2[2])
#P.add_constraint(rho_3_2==rho_list_2[3])
# P.add_constraint(rho_4_2==rho_list_2[4])
# # obs
Obs_list=get_obs(normal_forms=total_nfs, i=0,j=1)
# #print(H)
# for i in range(0,L-1,1):
#     mn=min(i,(i+1)%L)
#     mx=min(i,(i+1)%L)
#     rho_X=general_rho([mn,mx], Ms=Ms,normal_forms=total_nfs)
#     H_loc=mag_ham.MBL_chain(J, Delta,[hx[mn],hx[mx]] ).reshape(4,4)
#     Is_local=npa_hams.get_mbl_npa_local(normal_forms=total_nfs,i=mn,j=mx, J=J, hx=hx, hz=hz,Delta=Delta,shapes=shapes,L=L, PB=PB)
#     P.add_constraint(sum([(Ms[k] | Is_local[k]) for k in Ms.keys()])==(rho_X | H_loc ))

E=0.
P.add_constraint(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L<=E)
P.add_constraint(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L>=-E)
P.add_constraint(sum([(Ms[k] | Is_picos_squared[k]) for k in Ms.keys()])/L**2<=E**2)
P.add_constraint(sum([(Ms[k] | Is_picos_squared[k]) for k in Ms.keys()])/L**2>=0)
Exp=((rho_list[2] |  np.kron(mats[0],mats[2]) ))+(rho_list[2] | np.kron(mats[1],mats[2]))+(rho_list[2] | np.kron(mats[1],mats[3]))-(rho_list[2] | np.kron(mats[0],mats[3]))
s="x"
o_1=spin_class.spin_op(s,0)
s="x"
o_2=spin_class.spin_op(s,1)
OP_v=spin_class.op_list([o_1,o_2])
Iv_left=npa_hams.get_mbl_commleft_npa(normal_forms=total_nfs, J=J, hx=hx, hz=hz,shapes=shapes,L=L, OP=OP_v,PB=PB)
Iv_right=npa_hams.get_mbl_commright_npa(normal_forms=total_nfs, J=J, hx=hx, hz=hz,shapes=shapes,L=L, OP=OP_v,PB=PB)
P.add_constraint(sum([(Ms[k] | Iv_right[k]) for k in Ms.keys()])==sum([(Ms[k] | Iv_left[k]) for k in Ms.keys()]))

s="y"
o_1=spin_class.spin_op(s,0)
s="x"
o_2=spin_class.spin_op(s,1)
OP_v=spin_class.op_list([o_1,o_2])
Iv_left=npa_hams.get_mbl_commleft_npa(normal_forms=total_nfs, J=J, hx=hx, hz=hz,shapes=shapes,L=L, OP=OP_v,PB=PB)
Iv_right=npa_hams.get_mbl_commright_npa(normal_forms=total_nfs, J=J, hx=hx, hz=hz,shapes=shapes,L=L, OP=OP_v,PB=PB)
P.add_constraint(sum([(Ms[k] | Iv_right[k]) for k in Ms.keys()])==sum([(Ms[k] | Iv_left[k]) for k in Ms.keys()]))


s="z"
o_1=spin_class.spin_op(s,0)
s="y"
o_2=spin_class.spin_op(s,1)
OP_v=spin_class.op_list([o_1,o_2])
Iv_left=npa_hams.get_mbl_commleft_npa(normal_forms=total_nfs, J=J, hx=hx, hz=hz,shapes=shapes,L=L, OP=OP_v,PB=PB)
Iv_right=npa_hams.get_mbl_commright_npa(normal_forms=total_nfs, J=J, hx=hx, hz=hz,shapes=shapes,L=L, OP=OP_v,PB=PB)
P.add_constraint(sum([(Ms[k] | Iv_right[k]) for k in Ms.keys()])==sum([(Ms[k] | Iv_left[k]) for k in Ms.keys()]))








#P.set_objective("max",sum([(Ms[k] | Obs_list[k]) for k in Ms.keys()]))
P.set_objective("max",Exp)
#P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))

P.solve(solver="mosek")
print("obs",sum([(Ms[k] | Obs_list[k]) for k in Ms.keys()]).np)
print("chsh",Exp.np/(2))
# v, w, rho=test_ppt(rho_2.np)
# print(v)
v_red,w_red=np.linalg.eigh(rho_2.np)
print(v_red)
print("en",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]).np/L)
# #print("YY:", sum([(Ms[k] | Is_picos_XXX[k]) for k in Ms.keys()]).np/L)
# print(rho_list[2].np)


