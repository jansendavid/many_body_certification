import numpy as np
from copy import deepcopy
import spin_class as spin_class
import npa_funcs as npa_funcs
from npa_funcs import sym
import picos
import npa_hams as npa_hams
import sys
from marginal_problem import *

def test_make_normalform_vector():
    ops1=[]
    ops2=[]
    # ops contains spin list and list of symmetry sector
    L=3
    ops1+=[(spin_class.op_list([]),1)] # empty list gives unity
    for i in range(L):
        for s in sym:
            op=spin_class.op_list([spin_class.spin_op(s,i)])
            if s=="x":
                ops1+=[(op, 1)]
            else:
                ops2+=[(op, -1)]
    for i in range(L):
        s1="z"
        s2="x"
        op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i)])
        ops1+=[(op, 1)]
   # for i in range(L):
     
    # i=0
    # s1="z"
    # s2="x"
    # s3="y"
    # # op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i),spin_class.spin_op(s3,i)])
    # # print(op.sym)
    # # #ops1+=[(op, 1)]        
    # # coeff,nf=spin_class.normal_form(op)
    # # print(nf.sym, coeff)
    
    # i=0
    # s4="z"
    # op=spin_class.op_list([spin_class.spin_op(s1,i),spin_class.spin_op(s2,i),spin_class.spin_op(s3,i),spin_class.spin_op(s4,i)])
    # print(op.sym)
    # ops1+=[(op, 1)]        
    # #coeff,nf=spin_class.normal_form(op)
    # # print(nf.sym, coeff)
    
    
    vecs=npa_funcs.get_normal_form_vector(ops1,1, spin_class.spin_op, {})
    #vecs.update(npa_funcs.get_normal_form_vector(ops2, -1))
    for v in vecs:
        print(v.sym, vecs[v].get())
    return
#test_make_normalform_vector()
def test_make_Totalform_vector():
    ops={1: [], -1: []}

    # ops contains spin list and list of symmetry sector
    L=3
    ops[1]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(L):
        sec=-1
        for s in sym:
            op=spin_class.op_list([spin_class.spin_op(s,i)])
            if s=="x":
                ops[1]+=[op]
            else:
                ops[-1]+=[op]

    basis=npa_funcs.basis(ops, spin_class.spin_op, {})
    # vecs1=npa_funcs.get_normal_form_vector(ops1, 1,spin_class.spin_op, {})
    # vecs2=npa_funcs.get_normal_form_vector(ops2, -1,spin_class.spin_op, {})
    # total_vecs={+1: vecs1, -1:vecs2}
    # total_nfs=npa_funcs.get_NF_matrix_elements(total_vecs,spin_class.spin_op, {})
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
    L=3
    ops[1]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(L):
        sec=-1
        for s in sym:
            op=spin_class.op_list([spin_class.spin_op(s,i)])
            if s=="x":
                ops[1]+=[op]
            else:
                ops[-1]+=[op]

    basis=npa_funcs.basis(ops, spin_class.spin_op, {})
    P = picos.Problem()
    M=basis.get_mom_matrix_M(P)
    return


def test_tfi():
    ops={1: [], -1: []}

    # ops contains spin list and list of symmetry sector
    L=8
    ops[1]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(L):
        for s in sym:
            op=spin_class.op_list([spin_class.spin_op(s,i)])
            if s=="x":
                ops[1]+=[op]
            else:
                ops[-1]+=[op]
        
            for s2 in ["z"]:
                mn=min(i,(i+1)%L)
                mx=max(i,(i+1)%L)
                op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s2,mx)])
                if s=="x" and  not s2=="x":
                    ops[-1]+=[op]
                else:
                    ops[1]+=[op]
    
    d={}
    basis=npa_funcs.basis(ops=ops, class_name=spin_class.spin_op, **{})
    # for l in basis.total_nfs.keys():
    #     print("sec", l)
    #     for g in basis.total_vectors[l]:
    #         print(g.sym, basis.total_vectors[l][g].get())
    # print("end")
    # for l in basis.total_nfs.keys():
    #     print("sec", l)
    #     for g in basis.total_nfs[l]:
    #         print(g.sym, basis.total_nfs[l][g].get())
    # print("end")
    P = picos.Problem()
    Ms=basis.get_mom_matrix_M(P)
    shapes={}
    for k in Ms.keys():
        shapes[k]=Ms[k].shape
    PB=True
    J=4
    h=4
    Is=npa_hams.get_tfi_npa(normal_forms=basis.total_nfs, J=J, h=h,shapes=shapes,L=L, PB=PB)
    Is_picos={}
    for k in Is.keys():
        Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
        #print(Is[k])
    
    P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
    #print(P)
    P.solve(solver="mosek")
    #for k in Is.keys():
    #   print(Is[k])
    print(sum([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()]))
    return 



def test_xxz():
    
    ops={(1,1): [],(1,-1): [], (-1,1): [], (-1,-1): []}
    # ops contains spin list and list of symmetry sector
    L=4
    ops[(1,1)]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(L):
        s="z"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops[(1,-1)]+=[op]
        s="x"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops[(-1,1)]+=[op ]
        s="y"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops[(-1,-1)]+=[op]
    for i in range(L):
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)
        s="z"
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        #print(op.sym)/np.sqrt(L)
        ops[(1,1)]+=[op]
        s="x"
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops[(1,1)]+=[op]
        s="y"
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops[(1,1)]+=[op]


    basis=npa_funcs.basis(ops, spin_class.spin_op, **{})
    P = picos.Problem()
    Ms=basis.get_mom_matrix_M(P)
    shapes={}
    
    for k in Ms.keys():
        shapes[k]=Ms[k].shape
    PB=True
    J=1
    # h=4
    Is=npa_hams.get_xxz_npa(normal_forms=basis.total_nfs, J=J,Delta=0, shapes=shapes,L=L, PB=PB)
    Is_picos={}
    for k in Is.keys():
        Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
        #print(Is[k])
    
    P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
    #print(P)
    P.solve(solver="mosek")
    # for k in Is.keys():
    #     print(Is[k])
   
    print(sum([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()]))
    return 
# def test_xxz_2():
#     ops1=[]
#     ops2=[]
#     ops3=[]
#     # ops contains spin list and list of symmetry sector
#     L=8
#     ops1+=[(spin_class.op_list([]),(1,1))] # empty list gives unity
#     for i in range(L):
#         s="z"
#         op=spin_class.op_list([spin_class.spin_op(s,i)])
#         ops1+=[(op, (1,-1))]
#         s="x"
#         op=spin_class.op_list([spin_class.spin_op(s,i)])
#         ops1+=[(op, (-1,1))]
#         s="y"
#         op=spin_class.op_list([spin_class.spin_op(s,i)])
#         ops1+=[(op, (-1,-1))]
#     for i in range(L):
#         mn=min(i, (i+1)%L)
#         mx=max(i, (i+1)%L)
#         s="z"
#         op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
#         #print(op.sym)/np.sqrt(L)
#         ops1+=[(op, (1,1))]
#         s="x"
#         op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
#         ops1+=[(op, (1,1))]
#         s="y"
#         op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
#         ops1+=[(op, (1,1))]

#     vecs1=npa_funcs.get_normal_form_vector(ops1, (1,-1))
#     #vecs2=npa_funcs.get_normal_form_vector(ops2, (-1,1))
#     #vecs3=npa_funcs.get_normal_form_vector(ops3, (-1,-1))
    
#     total_vecs={(1,-1): vecs1}#, (-1,1):vecs2, (-1,-1): vecs3}
#     total_nfs=npa_funcs.get_NF_matrix_elements(total_vecs)
#     P = picos.Problem()
#     Ms=npa_funcs.get_mom_matrix_M(P, total_vecs, total_nfs)
#     shapes={}
    
#     for k in Ms.keys():
#         shapes[k]=Ms[k].shape
#     PB=True
#     J=1
#     # h=4
#     Is=npa_hams.get_xxz_npa(normal_forms=total_nfs, J=J,Delta=1, shapes=shapes,L=L, PB=PB)
#     Is_picos={}
#     for k in Is.keys():
#         Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
#         #print(Is[k])
    
#     P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
#     #print(P)
#     P.solve(solver="mosek")
#     for k in Is.keys():
#         print(Is[k])
    
#     #print(sum([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()]))
#     #print(Ms[(1,-1)].np[0,:])
#     #print(Ms[k].np)
#     return 
def test_xxz_highO():
    ops={(1,1): [],(1,-1): [], (-1,1): [], (-1,-1): []}
    # ops contains spin list and list of symmetry sector
    L=5
    ops[(1,1)]+=[spin_class.op_list([])] # empty list gives unity
    for i in range(L):
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)
        s="z"
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        #print(op.sym)
        ops[(1,1)]+=[op]
        s="x"
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops[(1,1)]+=[op]
        s="y"
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops[(1,1)]+=[op]
    for i in range(L):
        s="z"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops[(1,-1)]+=[op]
        s="x"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops[(-1,1)]+=[op]
        s="y"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops[(-1,-1)]+=[op]
        
    # for i in range(L):
    #     mn=min(i,(i+1)%L)
    #     mx=max(i,(i+1)%L)
    #     s1="x"
    #     s2="y"
    #     op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
    #     ops[(1,-1)]+=[op]
        
    #     op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
    #     ops[(1,-1)]+=[op]
    #     s1="y"
    #     s2="z"
    #     op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
    #     ops[(-1,1)]+=[op]
    #     op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
    #     ops[(-1,1)]+=[op]
    #     s1="x"
    #     s2="z"
    #     op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
    #     ops[(-1,-1)]+=[op]
    #     op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
    #     ops[(-1,-1)]+=[op]
        

    basis=npa_funcs.basis(ops, spin_class.spin_op, **{})
    P = picos.Problem()
    Ms=basis.get_mom_matrix_M(P)
    shapes={}
    shapes={}
    for k in Ms.keys():
        shapes[k]=Ms[k].shape
    PB=True
    J=1
    Delta=0
    # h=4
    Is=npa_hams.get_xxz_npa(normal_forms=basis.total_nfs, J=J, Delta=Delta,shapes=shapes,L=L, PB=PB)
    Is_picos={}
    for k in Is.keys():
        Is_picos[k]=picos.Constant("H{0}".format(k), Is[k])
    #

    P.set_objective("min",sum([(Ms[k] | Is_picos[k]) for k in Ms.keys()]))
    #print(P)
    P.solve(solver="mosek")
    #print([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()])
    #print(Ms[(1,1)].np)
    print(sum([(Ms[k] | Is_picos[k]).np/L for k in Ms.keys()]))
    return
#test_make_normalform_vector()
#test_make_set_up_momenMatrix()
test_tfi()
#test_xxz()
#test_xxz_highO()
#test_xxz_2()
