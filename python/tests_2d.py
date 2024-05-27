import numpy as np
from copy import deepcopy
import spin_class as spin_class
import spin_class_2d as spin_class_2d
import npa_funcs as npa_funcs
from npa_funcs import sym
import picos
import npa_hams_2d as npa_hams_2d
import sys
from marginal_problem import *

# def test_make_normalform_vector():
#     ops1=[]
#     ops2=[]
#     # ops contains spin list and list of symmetry sector
#     L=3
#     ops1+=[(spin_class.op_list([]),1)] # empty list gives unity
#     for i in range(L):
#         for s in sym:
#             op=spin_class.op_list([spin_class.spin_op(s,i)])
#             if s=="x":
#                 ops1+=[(op, 1)]
#             else:
#                 ops2+=[(op, -1)]
#     for i in range(L):
#         s1="z"
#         s2="x"
#         op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i)])
#         ops1+=[(op, 1)]
#    # for i in range(L):
     
#     # i=0
#     # s1="z"
#     # s2="x"
#     # s3="y"
#     # # op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i),spin_class.spin_op(s3,i)])
#     # # print(op.sym)
#     # # #ops1+=[(op, 1)]        
#     # # coeff,nf=spin_class.normal_form(op)
#     # # print(nf.sym, coeff)
    
#     # i=0
#     # s4="z"
#     # op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i),spin_class.spin_op(s3,i),spin_class.spin_op(s4,i)])
#     # print(op.sym)
#     # ops1+=[(op, 1)]        
#     # #coeff,nf=spin_class.normal_form(op)
#     # # print(nf.sym, coeff)
    
    
#     vecs=npa_funcs.get_normal_form_vector(ops1,1, spin_class.spin_op, {})
#     #vecs.update(npa_funcs.get_normal_form_vector(ops2, -1))
#     for v in vecs:
#         print(v.sym, vecs[v].get())
#     return
#test_make_normalform_vector()
def test_make_Totalform_vector():
    ops={1: [], -1: []}

    # ops contains spin list and list of symmetry sector
    Lx=3
    Ly=3
    ops[1]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(Lx):
        for j in range(Lx):
            for s in sym:
                op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,j],Lx= Lx)])
                if s=="x":
                    ops[1]+=[op]
                else:
                    ops[-1]+=[op]

    basis=npa_funcs.basis(ops, spin_class_2d.spin_op_2d, **{"Lx": Lx})
    for v in basis.total_nfs.keys():
        print("sym  sec:", v)
        curr_sec=basis.total_nfs[v]
        for w in curr_sec.keys():
            print(w.sym, curr_sec[w].get())
    return
#test_make_Totalform_vector()
def test_make_set_up_momenMatrix():
    ops={1: [], -1: []}

    # ops contains spin list and list of symmetry sector
    Lx=3
    Ly=3
    ops[1]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(Lx):
        for j in range(Lx):
            for s in sym:
                op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,j],Lx= Lx)])
                if s=="x":
                    ops[1]+=[op]
                else:
                    ops[-1]+=[op]

    basis=npa_funcs.basis(ops, spin_class_2d.spin_op_2d, **{"Lx": Lx})
    P = picos.Problem()
    M=basis.get_mom_matrix_M(P)
    return

#test_make_set_up_momenMatrix()
def test_tfi():
    ops={1: [], -1: []}

    # ops contains spin list and list of symmetry sector
    Lx=3
    Ly=3
    ops[1]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(Lx):
        for j in range(Ly):
            for s in sym:
                op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,j],Lx= Lx)])
                if s=="x":
                    ops[1]+=[op]
                else:
                    ops[-1]+=[op]

    basis=npa_funcs.basis(ops, spin_class_2d.spin_op_2d, **{"Lx": Lx})
    P = picos.Problem()
    Ms=basis.get_mom_matrix_M(P)
    shapes={}
    for k in Ms.keys():
        shapes[k]=Ms[k].shape
    PB=True
    J=4
    h=4
    Is=npa_hams_2d.get_tfi2d_npa(normal_forms=basis.total_nfs, J=J, h=h,shapes=shapes,Ls=[Lx,Ly], PB=PB)
    Is_picos={}
    for k in Is.keys():
        Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
        #print(Is[k])
    
    P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
    #print(P)
    P.solve(solver="mosek")
    #for k in Is.keys():
    #   print(Is[k])
    print(sum([(Ms[k] | Is_picos[k]).np/(Lx*Ly) for k in Ms.keys()]))
    return 

#test_tfi()

def test_xxz():
    
    ops={(1,1): [],(1,-1): [], (-1,1): [], (-1,-1): []}
    # ops contains spin list and list of symmetry sector
    Lx=3
    Ly=3
    ops[(1,1)]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(Lx):
        for j in range(Ly):
            s="z"
            op=spin_class.op_list([spin_class.spin_op(s,[i,j])])
            ops[(1,-1)]+=[op]
            s="x"
            op=spin_class.op_list([spin_class.spin_op(s,[i,j])])
            ops[(-1,1)]+=[op ]
            s="y"
            op=spin_class.op_list([spin_class.spin_op(s,[i,j])])
            ops[(-1,-1)]+=[op]
    for i in range(Lx):
        for j in range(Ly):
            mn_x=min(i, (i+1)%Lx)
            mx_x=max(i, (i+1)%Lx)
            mn_y=min(j, (j+1)%Ly)
            mx_y=max(j, (j+1)%Ly)
            s="z"
            op=spin_class.op_list([spin_class.spin_op(s,[mn_x, j]),spin_class.spin_op(s,[mx_x, j])])
            ops[(1,1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s,[i, mn_y]),spin_class.spin_op(s,[i, mx_y])])
            ops[(1,1)]+=[op]
            s="x"
            op=spin_class.op_list([spin_class.spin_op(s,[mn_x, j]),spin_class.spin_op(s,[mx_x, j])])
            ops[(1,1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s,[i, mn_y]),spin_class.spin_op(s,[i, mx_y])])
            ops[(1,1)]+=[op]
            s="y"
            op=spin_class.op_list([spin_class.spin_op(s,[mn_x, j]),spin_class.spin_op(s,[mx_x, j])])
            ops[(1,1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s,[i, mn_y]),spin_class.spin_op(s,[i, mx_y])])
            ops[(1,1)]+=[op]


    basis=npa_funcs.basis(ops, spin_class_2d.spin_op_2d, **{"Lx": Lx})
    P = picos.Problem()
    Ms=basis.get_mom_matrix_M(P)
    shapes={}
    
    for k in Ms.keys():
        shapes[k]=Ms[k].shape
    PB=True
    J=1
    # h=4
    Is=npa_hams_2d.get_xxz2d_npa(normal_forms=basis.total_nfs, J=J,Delta=1, shapes=shapes,Ls=[Lx,Ly], PB=PB)
    Is_picos={}
    for k in Is.keys():
        Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
        #print(Is[k])
    
    P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
    #print(P)
    P.solve(solver="mosek")
    # for k in Is.keys():
    #     print(Is[k])
   
    print(sum([(Ms[k] | Is_picos[k]).np/(Lx*Ly) for k in Ms.keys()]))
    return 
test_xxz()

def test_xxz_highO():
  
    ops={(1,1): [],(1,-1): [], (-1,1): [], (-1,-1): []}
    # ops contains spin list and list of symmetry sector
    Lx=3
    Ly=3
    ops[(1,1)]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(Lx):
        for j in range(Ly):
            s="z"
            op=spin_class.op_list([spin_class.spin_op(s,[i,j])])
            ops[(1,-1)]+=[op]
            s="x"
            op=spin_class.op_list([spin_class.spin_op(s,[i,j])])
            ops[(-1,1)]+=[op ]
            s="y"
            op=spin_class.op_list([spin_class.spin_op(s,[i,j])])
            ops[(-1,-1)]+=[op]
    for i in range(Lx):
        for j in range(Ly):
            mn_x=min(i, (i+1)%Lx)
            mx_x=max(i, (i+1)%Lx)
            mn_y=min(j, (j+1)%Ly)
            mx_y=max(j, (j+1)%Ly)
            s="z"
            op=spin_class.op_list([spin_class.spin_op(s,[mn_x, j]),spin_class.spin_op(s,[mx_x, j])])
            ops[(1,1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s,[i, mn_y]),spin_class.spin_op(s,[i, mx_y])])
            ops[(1,1)]+=[op]
            s="x"
            op=spin_class.op_list([spin_class.spin_op(s,[mn_x, j]),spin_class.spin_op(s,[mx_x, j])])
            ops[(1,1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s,[i, mn_y]),spin_class.spin_op(s,[i, mx_y])])
            ops[(1,1)]+=[op]
            s="y"
            op=spin_class.op_list([spin_class.spin_op(s,[mn_x, j]),spin_class.spin_op(s,[mx_x, j])])
            ops[(1,1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s,[i, mn_y]),spin_class.spin_op(s,[i, mx_y])])
            ops[(1,1)]+=[op]
    for i in range(Lx):
        for j in range(Ly):
            mn_x=min(i, (i+1)%Lx)
            mx_x=max(i, (i+1)%Lx)
            mn_y=min(j, (j+1)%Ly)
            mx_y=max(j, (j+1)%Ly)
            s1="x"
            s2="y"
            op=spin_class.op_list([spin_class.spin_op(s1,[mn_x, j]),spin_class.spin_op(s2,[mx_x, j])])
            ops[(1,-1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s1,[i, mn_y]),spin_class.spin_op(s2,[i, mx_y])])
            ops[(1,-1)]+=[op]
        
            op=spin_class.op_list([spin_class.spin_op(s2,[mn_x, j]),spin_class.spin_op(s1,[mx_x, j])])
            ops[(1,-1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s2,[i, mn_y]),spin_class.spin_op(s1,[i, mx_y])])
            ops[(1,-1)]+=[op]
            
            s1="y"
            s2="z"
            op=spin_class.op_list([spin_class.spin_op(s1,[mn_x, j]),spin_class.spin_op(s2,[mx_x, j])])
            ops[(1,-1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s1,[i, mn_y]),spin_class.spin_op(s2,[i, mx_y])])
            ops[(1,-1)]+=[op]
        
            op=spin_class.op_list([spin_class.spin_op(s2,[mn_x, j]),spin_class.spin_op(s1,[mx_x, j])])
            ops[(1,-1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s2,[i, mn_y]),spin_class.spin_op(s1,[i, mx_y])])
            ops[(1,-1)]+=[op]
            s1="x"
            s2="z"
            op=spin_class.op_list([spin_class.spin_op(s1,[mn_x, j]),spin_class.spin_op(s2,[mx_x, j])])
            ops[(1,-1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s1,[i, mn_y]),spin_class.spin_op(s2,[i, mx_y])])
            ops[(1,-1)]+=[op]
        
            op=spin_class.op_list([spin_class.spin_op(s2,[mn_x, j]),spin_class.spin_op(s1,[mx_x, j])])
            ops[(1,-1)]+=[op]
            op=spin_class.op_list([spin_class.spin_op(s2,[i, mn_y]),spin_class.spin_op(s1,[i, mx_y])])
            ops[(1,-1)]+=[op]
    basis=npa_funcs.basis(ops, spin_class_2d.spin_op_2d, **{"Lx": Lx})
    P = picos.Problem()
    Ms=basis.get_mom_matrix_M(P)
    shapes={}
    
    for k in Ms.keys():
        shapes[k]=Ms[k].shape
    PB=True
    J=1
    # h=4
    Is=npa_hams_2d.get_xxz2d_npa(normal_forms=basis.total_nfs, J=J,Delta=1, shapes=shapes,Ls=[Lx,Ly], PB=PB)
    Is_picos={}
    for k in Is.keys():
        Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
        #print(Is[k])
    
    P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
    #print(P)
    P.solve(solver="cvxopt")
    # for k in Is.keys():
    #     print(Is[k])
   
    print(sum([(Ms[k] | Is_picos[k]).np/(Lx*Ly) for k in Ms.keys()]))
        
    return
#test_make_normalform_vector()
#test_make_set_up_momenMatrix()
#test_tfi()
#test_xxz()
test_xxz_highO()#
#test_xxz_2()
