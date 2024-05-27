#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 14:11:16 2024

@author: david
"""
import numpy as np
import sys
sys.path.insert(0,"../")
from copy import deepcopy
import spin_class as spin_class
import spin_class_2d as spin_class_2d
import npa_funcs as npa_funcs
from npa_funcs import sym
import picos
import npa_hams as npa_hams
import sys
import networkx as nx
import momentum_funcs as mfuncs
import os 
import npa_hams_2d as npa_hams_2d
def get_xxz2d(obj, J, Delta, Ly,Lx):
    Is={}
    Is[0]=np.zeros(obj.blocks[0].shape)
    sym=["x", "y", "z"]
    #sym=["z"]
    for j in range(Ly):
        
        for s in sym:  
            if s=="z":
                
                xxz_c=Delta
            else:
                xxz_c=J
            o_1=spin_class_2d.spin_op_2d(s,[min(j,(j+1)%Ly),0], Lx=Lx)
            o_2=spin_class_2d.spin_op_2d(s,[max(j,(j+1)%Ly),0], Lx=Lx)
            otot=spin_class.op_list([o_1,o_2])
            coeff,nf=spin_class.normal_form(otot, class_name=spin_class_2d.spin_op_2d, **{"Lx": Lx})
            if nf in obj.d_vector_map.keys():
                #print(nf.sym)
                co=obj.d_vector_map[nf]
                #print(nf.sym, co.get())
                Is[0][0,co.ind]=xxz_c/4 # plus one because it is in block [1,d,d
                Is[0][co.ind, 0]=xxz_c/4
    
        for s in sym:  
            if s=="z":
                xxz_c=Delta
            else:
                xxz_c=J
            o_1=spin_class_2d.spin_op_2d(s,[j,0], Lx=Lx)
            o_2=spin_class_2d.spin_op_2d(s,[j,1], Lx=Lx)
            otot=spin_class.op_list([o_1,o_2])
            coeff,nf=spin_class.normal_form(otot, class_name=spin_class_2d.spin_op_2d, **{"Lx": Lx})
            if nf in obj.d_vector_map.keys():
                #print(nf.sym)
                co=obj.d_vector_map[nf]
                #print(nf.sym, co.get())
                Is[0][0,co.ind]=xxz_c/4 # plus one because it is in block [1,d,d
                Is[0][co.ind, 0]=xxz_c/4
    #print(Is[0])       
    return Is[0]/2
def get_moments_1(Lx,Ly):
    ops0=[] # (1,1)
    ops1=[] #(1,-1)
    ops2=[] #(-1,1)
    ops3=[] #(-1,-1)
    for i in range(Ly):
        mn=min(i, (i+1)%Ly)
        mx=max(i, (i+1)%Ly)
        s="z"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops0+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[mx,0], Lx=Lx)])
        ops0+=[op]
        s="x"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops0+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[mx,0], Lx=Lx)])
        ops0+=[op]
        s="y"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops0+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[mx,0], Lx=Lx)])
        ops0+=[op]
    for i in range(Ly):
        s="z"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        ops1+=[op]
        s="x"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        ops2+=[op]
        s="y"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        ops3+=[op]
        
    for i in range(Ly):
        mn=min(i,(i+1)%Ly)
        mx=max(i,(i+1)%Ly)
        s1="x"
        s2="y"
      
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
        ops1+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
        ops1+=[op]
        
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
        ops1+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
        ops1+=[op]
        s1="y"
        s2="z"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
        ops2+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
        ops2+=[op]
        
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
        ops2+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
        ops2+=[op]
        s1="x"
        s2="z"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
        ops3+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
        ops3+=[op]
        
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
        ops3+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
        ops3+=[op]
    operators={(1,1): ops0,(1,-1): ops1,(-1,1): ops2,(-1,-1): ops3}
    return operators
def get_moments_2(Lx,Ly):
    ops0=[] # (1,1)
    ops1=[] #(1,-1)
    ops2=[] #(-1,1)
    ops3=[] #(-1,-1)
    for i in range(Ly):
        mn=min(i, (i+1)%Ly)
        mx=max(i, (i+1)%Ly)
        s="z"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops0+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[mx,0], Lx=Lx)])
        ops0+=[op]
        s="x"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops0+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[mx,0], Lx=Lx)])
        ops0+=[op]
        s="y"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops0+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[mx,0], Lx=Lx)])
        ops0+=[op]
    for i in range(Ly):
        s="z"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        ops1+=[op]
        s="x"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        ops2+=[op]
        s="y"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        ops3+=[op]
        
    # for i in range(Ly):
    #     mn=min(i,(i+1)%Ly)
    #     mx=max(i,(i+1)%Ly)
    #     s1="x"
    #     s2="y"
      
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
    #     ops1+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
    #     ops1+=[op]
        
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
    #     ops1+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
    #     ops1+=[op]
    #     s1="y"
    #     s2="z"
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
    #     ops2+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
    #     ops2+=[op]
        
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
    #     ops2+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
    #     ops2+=[op]
    #     s1="x"
    #     s2="z"
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
    #     ops3+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
    #     ops3+=[op]
        
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
    #     ops3+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
    #     ops3+=[op]
    operators={(1,1): ops0,(1,-1): ops1,(-1,1): ops2,(-1,-1): ops3}
    return operators
def get_moments_3(Lx,Ly):
    ops0=[] # (1,1)
    ops1=[] #(1,-1)
    ops2=[] #(-1,1)
    ops3=[] #(-1,-1)
    for i in range(Ly):
        mn=min(i, (i+1)%Ly)
        mx=max(i, (i+1)%Ly)
        s="z"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops0+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[mx,0], Lx=Lx)])
        ops0+=[op]
        s="x"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops0+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[mx,0], Lx=Lx)])
        ops0+=[op]
        s="y"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops0+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[mx,0], Lx=Lx)])
        ops0+=[op]
   # for i in range(Ly):
    #    s="z"
     #   op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
      #  ops1+=[op]
       # s="x"
       # op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        #ops2+=[op]
        #s="y"
        #op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        #ops3+=[op]
        
    # for i in range(Ly):
    #     mn=min(i,(i+1)%Ly)
    #     mx=max(i,(i+1)%Ly)
    #     s1="x"
    #     s2="y"
      
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
    #     ops1+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
    #     ops1+=[op]
        
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
    #     ops1+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
    #     ops1+=[op]
    #     s1="y"
    #     s2="z"
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
    #     ops2+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
    #     ops2+=[op]
        
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
    #     ops2+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
    #     ops2+=[op]
    #     s1="x"
    #     s2="z"
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
    #     ops3+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
    #     ops3+=[op]
        
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
    #     ops3+=[op]
    #     op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
    #     ops3+=[op]
    operators={(1,1): ops0}
    return operators
#     # ops contains spin list and list of symmetry sector
Lx=8
Ly=8
Delta=.0
J=1
Moms=2
if Moms==0:
    operators=get_moments_1(Lx,Ly)
elif Moms==1:
    operators=get_moments_2(Lx,Ly)
elif Moms==2:
    operators=get_moments_3(Lx,Ly)
#dn="2dHeisenberJ{0:.3f}Delta{1:.3f}Lx{2}Ly{3}Moms{4}".format(J,Delta,Lx,Ly,Moms)
#os.mkdir(dn)
    
P = picos.Problem()
#     # for o in ops0:
#     #     print(o.sym)
#     # for o in ops1:
#     #     print(o.sym)


kwargs={"Lx":Lx}
basis=mfuncs.momentum_basis( operators, Lx, P, spin_class_2d.spin_op_2d,**kwargs)


Is=Is=get_xxz2d(basis.sectors[(1,1)], J=J, Delta=Delta,Ly=Ly, Lx=Lx)
#print(Is)
Is_picos={}
#print(P)
Is_picos=picos.Constant("H{0}".format(0), Is)

P.set_objective("min",(basis.sectors[(1,1)].blocks[0] | Is_picos))

P.solve(solver="mosek")
E=np.array([(basis.sectors[(1,1)].blocks[0] | Is_picos).np/(Ly*np.sqrt(Lx))])
B0=(basis.sectors[(1,1)].blocks[0]).np
#np.save(dn+"/E"  ,E)
#np.save(dn+"/B0",B0 )
print((basis.sectors[(1,1)].blocks[0] | Is_picos).np/(Ly*np.sqrt(Lx)))