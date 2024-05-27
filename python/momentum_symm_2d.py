import numpy as np
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
import npa_hams_2d as npa_hams_2d
def test_tranlation():
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

    #basis=npa_funcs.basis(ops, spin_class_2d.spin_op_2d, {"Lx": Lx})
    for op in ops[1]:
        print("start")
        print(op.sym)
        for j in range(1,Lx):
            op_new=mfuncs.translate_2d(op, j, L=Lx, Lx=Lx,index=0)
            print(op_new.sym)
        print("end")
    print("fin e")
    for op in ops[1]:
        print("start")
        print(op.sym)
        for j in range(1,Ly):
            op_new=mfuncs.translate_2d(op, j, L=Ly,Lx=Lx, index=1)
            print(op_new.sym·)
        print("end")
def get_xxz2d(obj, J, Ly,Lx):
    Is={}
    Is[0]=np.zeros(obj.blocks[0].shape)
    sym=["x", "y", "z"]
    #sym=["z"]
    for j in range(Ly):
        
        for s in sym:          
            o_1=spin_class_2d.spin_op_2d(s,[min(j,(j+1)%Ly),0], Lx=Lx)
            o_2=spin_class_2d.spin_op_2d(s,[max(j,(j+1)%Ly),0], Lx=Lx)
            otot=spin_class.op_list([o_1,o_2])
            coeff,nf=spin_class.normal_form(otot, class_name=spin_class_2d.spin_op_2d, **{"Lx": Lx})
            if nf in obj.d_vector_map.keys():
                #print(nf.sym)
                co=obj.d_vector_map[nf]
                #print(nf.sym, co.get())
                Is[0][0,co.ind]=J/4 # plus one because it is in block [1,d,d
                Is[0][co.ind, 0]=J/4
    
        for s in sym:          
            o_1=spin_class_2d.spin_op_2d(s,[j,0], Lx=Lx)
            o_2=spin_class_2d.spin_op_2d(s,[j,1], Lx=Lx)
            otot=spin_class.op_list([o_1,o_2])
            coeff,nf=spin_class.normal_form(otot, class_name=spin_class_2d.spin_op_2d, **{"Lx": Lx})
            if nf in obj.d_vector_map.keys():
                #print(nf.sym)
                co=obj.d_vector_map[nf]
                #print(nf.sym, co.get())
                Is[0][0,co.ind]=J/4 # plus one because it is in block [1,d,d
                Is[0][co.ind, 0]=J/4
                
    return Is[0]/2

                
def test_find_blocks():
    Lx=6
    Ly=1
    ops=[]
    for i in range(Ly):
        s="z"   
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        ops+=[op]
        s="y" 
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        ops+=[op]
        s="x" 
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx)])
        ops+=[op]
    for i in range(Ly):
        s="z"  
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[min(i, (i+1)%Ly),0], Lx=Lx),spin_class_2d.spin_op_2d(s,[max(i, (i+1)%Ly),0], Lx=Lx)])
        #ops+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops+=[op]
        s="y"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[min(i, (i+1)%Ly),0], Lx=Lx),spin_class_2d.spin_op_2d(s,[max(i, (i+1)%Ly),0], Lx=Lx)])
        #ops+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops+=[op]
        s="x"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[min(i, (i+1)%Ly),0], Lx=Lx),spin_class_2d.spin_op_2d(s,[max(i, (i+1)%Ly),0], Lx=Lx)])
        #ops+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s,[i,1], Lx=Lx)])
        ops+=[op]
    for i in range(Ly):
        mn=min(i,(i+1)%Ly)
        mx=max(i,(i+1)%Ly)
        s1="x"
        s2="y"
  
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
        ops+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
        #ops1+=[op]
        
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
        ops+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
        #ops+=[op]
        s1="y"
        s2="z"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
        ops+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
        #ops2+=[op]
        
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
        ops+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
        #ops2+=[op]
        s1="x"
        s2="z"
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[i,1], Lx=Lx)])
        ops+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s1,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s2,[mx,0], Lx=Lx)])
        #ops3+=[op]
        
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[i,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[i,1], Lx=Lx)])
        ops+=[op]
        op=spin_class.op_list([spin_class_2d.spin_op_2d(s2,[mn,0], Lx=Lx),spin_class_2d.spin_op_2d(s1,[mx,0], Lx=Lx)])
        #ops3+=[op]
#     #for s in ops1:
#     #    print(s.sym)
#     #print("start")
    P = picos.Problem()
    #print(ops1)
    M=M=mfuncs.mom_symm_block_zero(ops, Lx,P, class_name=spin_class_2d.spin_op_2d,sec="",**{"Lx":Lx})
    #mfuncs.mom_symm_block_zero(ops, Ls=[Lx,Ly],P=P)
#     # for n in M.map_TI_to_elements.keys():
#     #     print(n.sym, "->", M.map_TI_to_elements[n])
#     # print("start")
#     # for n in M.map_TI.keys():
#     #     print(n.sym, "->", M.map_TI[n].sym)
#     # print("d's")
#     # for n in M.d_vector_map.keys():
#     #     print(n.sym, "->", M.d_vector_map[n].get())
    Is=get_xxz2d(M, J=1, Ly=Ly, Lx=Lx)
    print(Is)
    
    Is_picos={}
    #print(P)
    Is_picos=picos.Constant("H{0}".format(0), Is)

    P.set_objective("min",(M.blocks[0] | Is_picos))
    print("start solve")
    P.solve(solver="mosek")
    print((M.blocks[0] | Is_picos).np/(Ly*np.sqrt(Lx)))
#     print(M.blocks[0].np[0,:])
#test_find_blocks()    
def test_find_many_blocks():
    ops0=[] # (1,1)
    ops1=[] #(1,-1)
    ops2=[] #(-1,1)
    ops3=[] #(-1,-1)
#     # ops contains spin list and list of symmetry sector
    Lx=6
    Ly=6
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

        
    P = picos.Problem()
#     # for o in ops0:
#     #     print(o.sym)
#     # for o in ops1:
#     #     print(o.sym)

    operators={(1,1): ops0,(1,-1): ops1,(-1,1): ops2,(-1,-1): ops3}
    kwargs={"Lx":Lx}
    basis=mfuncs.momentum_basis( operators, Lx, P, spin_class_2d.spin_op_2d,**kwargs)
    Is=Is=get_xxz2d(basis.sectors[(1,1)], J=1, Ly=Ly, Lx=Lx)
    #print(Is)
    Is_picos={}
    #print(P)
    Is_picos=picos.Constant("H{0}".format(0), Is)

    P.set_objective("min",(basis.sectors[(1,1)].blocks[0] | Is_picos))

    P.solve(solver="mosek")
    print((basis.sectors[(1,1)].blocks[0] | Is_picos).np/(Ly*np.sqrt(Lx)))
#test_tranlation()

test_find_many_blocks()





