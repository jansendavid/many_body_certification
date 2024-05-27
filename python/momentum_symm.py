import numpy as np
from copy import deepcopy
import spin_class as spin_class
import npa_funcs as npa_funcs
from npa_funcs import sym
import picos
import npa_hams as npa_hams
import sys
import networkx as nx
import momentum_funcs as mfuncs

def test_tranlation():
    L=4
    ops1=[]
    ops1=[]
    ops2=[]
    ops3=[]
    for i in range(L):
        s="z"
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)    
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops1+=[(op, (1,1))]
    for op in ops1:
        print("start")
        print(op[0].sym)
        for j in range(1,L):
            op_new=mfuncs.translate(op[0], j, L)
            print(op_new.sym)
        print("end")
  
def get_xxz(obj, J,Delta):
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
            if s=="z":
                Is[0][0,co.ind]=Delta/4 # plus one because it is in block [1,d,d
                Is[0][co.ind, 0]=Delta/4
            else:
                Is[0][0,co.ind]=J/4 # plus one because it is in block [1,d,d
                Is[0][co.ind, 0]=J/4
    return Is[0]/2

                
def test_find_blocks():
    L=20
    ops1=[]
    ops1=[]
    ops2=[]
    ops3=[]
    for i in range(1):
        s="z"
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)    
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[op]
        s="y"
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)    
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[op]
        s="x"
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)    
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[op]
    for i in range(1):
        s="z"
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)    
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops1+=[op]
        s="y"
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)    
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops1+=[op]
        s="x"
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)    
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops1+=[op]
    for i in range(1):
        mn=min(i,(i+1)%L)
        mx=max(i,(i+1)%L)
        s1="x"
        s2="y"
        op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
        ops1+=[op]
        
        op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
        ops1+=[op]
        s1="y"
        s2="z"
        op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
        ops1+=[op]
        op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
        ops1+=[op]
        s1="x"
        s2="z"
        op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
        ops1+=[op]
        op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
        ops1+=[op]
    #for s in ops1:
    #    print(s.sym)
    #print("start")
    P = picos.Problem()
    #print(ops1)
    M=mfuncs.mom_symm_block_zero(ops1, L,P, class_name=spin_class.spin_op,**{})
    # for n in M.map_TI_to_elements.keys():
    #     print(n.sym, "->", M.map_TI_to_elements[n])
    # print("start")
    # for n in M.map_TI.keys():
    #     print(n.sym, "->", M.map_TI[n].sym)
    # print("d's")
    # for n in M.d_vector_map.keys():
    #     print(n.sym, "->", M.d_vector_map[n].get())
    Is=get_xxz(M, J=1)
    #print(Is)
    
    Is_picos={}
    #print(P)
    Is_picos=picos.Constant("H{0}".format(0), Is)

    P.set_objective("min",(M.blocks[0] | Is_picos))
    print("start solve")
    P.solve(solver="mosek")
    print((M.blocks[0] | Is_picos).np/np.sqrt(L))
    print(M.blocks[0].np[0,:])
#test_find_blocks()    
def test_find_many_blocks():
    ops0=[] # (1,1)
    ops1=[] #(1,-1)
    ops2=[] #(-1,1)
    ops3=[] #(-1,-1)
    # ops contains spin list and list of symmetry sector
    L=4
    for i in range(1):
        mn=min(i, (i+1)%L)
        mx=max(i, (i+1)%L)
        s="z"
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        #print(op.sym)
        ops0+=[op]
        s="x"
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops0+=[op]
        s="y"
        op=spin_class.op_list([spin_class.spin_op(s,mn),spin_class.spin_op(s,mx)])
        ops0+=[op]
    for i in range(1):
        s="z"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[op]
        s="x"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops2+=[op]
        s="y"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops3+=[op]
        
    for i in range(1):
        mn=min(i,(i+1)%L)
        mx=max(i,(i+1)%L)
        s1="x"
        s2="y"
        op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
        ops1+=[op]
        
        op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
        ops1+=[op]
        s1="y"
        s2="z"
        op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
        ops2+=[op]
        op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
        ops2+=[op]
        s1="x"
        s2="z"
        op=spin_class.op_list([spin_class.spin_op(s1,mn),spin_class.spin_op(s2,mx)])
        ops3+=[op]
        op=spin_class.op_list([spin_class.spin_op(s2,mn),spin_class.spin_op(s1,mx)])
        ops3+=[op]
        
        
    P = picos.Problem()
    # for o in ops0:
    #     print(o.sym)
    # for o in ops1:
    #     print(o.sym)
    operators={(1,1): ops0,(1,-1): ops1,(-1,1): ops2,(-1,-1): ops3}
    kwargs={}
    basis=mfuncs.momentum_basis( operators, L, P, spin_class.spin_op,**{})
    Is=get_xxz( basis.sectors[(1,1)], J=1, Delta=0)
    
    #print(Is)

    Is_picos={}
    #print(P)
    Is_picos=picos.Constant("H{0}".format(0), Is)

    P.set_objective("min",(basis.sectors[(1,1)].blocks[0] | Is_picos))

    P.solve(solver="mosek")
    print((basis.sectors[(1,1)].blocks[0] | Is_picos).np/np.sqrt(L))
    # print(M0.blocks[0].np[0,:])
#    # print(M.blocks[0].np)
#     #for n in M.blocks.keys():
#     #    print(M.blocks[n].np)
#     #for m in M.G_matrices.keys():
#     #    print(m)
#     #for v in M.D_vector_map.keys():
#     #    print(v.sym)
# #test_tranlation()

test_find_many_blocks()





