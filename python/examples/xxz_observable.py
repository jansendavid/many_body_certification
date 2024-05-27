#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 14:43:07 2024

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
from generate_monomials import prepare_nfs_xxzsign_symmetrie
def get_obs(normal_forms, i,j):
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k], dtype=complex)
    s1="x"
    s2="x"        
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
    
L=10
PB=True
J=1
Delta=1
h=0
E_upper=-0.4431471805599453
E_lower=-0.45539828492925705
total_vecs,total_nfs=prepare_nfs_xxzsign_symmetrie(L)
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
m=6
rho_list=solve_marginal_problem(P,  m)
#rho_1=get_rho_1(i=0,Ms=Ms,normal_forms=total_nfs)
#rho_2=get_rho_2(0,1, Ms=Ms,normal_forms=total_nfs)
rho_1=general_rho([0], Ms=Ms,normal_forms=total_nfs)

rho_2=general_rho([0,1], Ms=Ms,normal_forms=total_nfs)

rho_3=general_rho([0,1,2], Ms=Ms,normal_forms=total_nfs)

rho_4=general_rho([0,1,2,3], Ms=Ms,normal_forms=total_nfs)

rho_list[1]=picos.SymmetricVariable("rho{0}".format(1), (2,2))
P.add_constraint(rho_list[1]>>0)
P.add_constraint(rho_4==rho_list[4])
P.add_constraint(rho_3==rho_list[3])
P.add_constraint(rho_2==rho_list[2])
P.add_constraint(rho_1==rho_list[1])
P.add_constraint(picos.partial_trace(rho_list[2], 0,2)==rho_list[1])
P.add_constraint(picos.partial_trace(rho_list[2], 1,2)==rho_list[1])
H=mag_ham.XXZ_chain(J=J,Delta=Delta,h=h).reshape(4,4)
#print(H)
H_picos=picos.Constant("H", H)
Obs_list=get_obs(normal_forms=total_nfs, i=0,j=1)
P.add_constraint(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L==(H_picos| rho_list[2]))
#P.add_constraint(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L==(H_picos@H_picos| rho_list[4])/2)

P.add_constraint(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L<=E_upper)
P.add_constraint(sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()])/L>=E_lower)


Oop=np.kron(np.eye(2), sx)+np.kron( sx,np.eye(2))
#P.add_constraint(Oop*rho_list[2]-rho_list[2]*Oop==np.zeros(shape=(4,4)))
Oop=np.kron(np.eye(2), sz)+np.kron( sz,np.eye(2))
#P.add_constraint(Oop*rho_list[2]-rho_list[2]*Oop==np.zeros(shape=(4,4)))
#P.add_constraint(picos.trace(Oop*rho_list[2])-picos.trace(rho_list[2]*Oop)==0)
Oop=np.kron(np.eye(2), sy)+np.kron( sy,np.eye(2))
#P.add_constraint(Oop*rho_list[2]-rho_list[2]*Oop==np.zeros(shape=(4,4)))

Oop=np.kron(np.eye(2), np.kron(np.eye(2),sx))+np.kron( np.eye(2),np.kron(sx,np.eye(2)))+np.kron( sx,np.kron(np.eye(2),np.eye(2)))
#P.add_constraint(Oop*rho_list[3]-rho_list[3]*Oop==np.zeros(shape=(8,8)))
Oop=np.kron(np.eye(2), np.kron(np.eye(2),sz))+np.kron( np.eye(2),np.kron(sz,np.eye(2)))+np.kron( sz,np.kron(np.eye(2),np.eye(2)))
#P.add_constraint(Oop*rho_list[3]-rho_list[3]*Oop==np.zeros(shape=(8,8)))
#Oop=np.kron(np.eye(2), np.kron(np.eye(2),sy))+np.kron( np.eye(2),np.kron(sy,np.eye(2)))+np.kron( sy,np.kron(np.eye(2),np.eye(2)))
#P.add_constraint(Oop*rho_list[3]-rho_list[3]*Oop==np.zeros(shape=(8,8)))




Orand=np.random.rand(2,2)
Orand=Orand+np.transpose(Orand)
Orand=np.kron(np.eye(2), Orand)+np.kron(Orand,np.eye(2))
P.add_constraint(picos.trace(rho_list[2]*Orand)==picos.trace(Orand*rho_list[2]))


#P.set_objective("min",(H_picos| rho_list[2]))
P.set_objective("min",sum([(Ms[k] | Obs_list[k]) for k in Ms.keys()]))

P.solve(solver="mosek")
    #print([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()])
print(sum([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()])) 
print((H_picos| rho_list[2])) 
print("obs",sum([(Ms[k] | Obs_list[k]) for k in Ms.keys()]))
#print(np.trace(np.kron(H, np.eye(4))+np.kron(np.eye(4),H ), rho_list[4].np)/2)