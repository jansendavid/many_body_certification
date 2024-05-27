import numpy as np
from copy import deepcopy
import spin_class as spin_class
import npa_funcs as npa_funcs
from npa_funcs import sym
import picos
import npa_hams as npa_hams
import sys
from marginal_problem import *
import marginal_problem_hamiltonians as mag_ham

def test_rhos_gen():
    ops1=[]
    ops2=[]
    ops3=[]
    # ops contains spin list and list of symmetry sector
    L=3
    
    for i in range(L):
        s="z"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops1+=[(op, (1,-1))]
        s="x"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops2+=[(op, (-1,1))]
        s="y"
        op=spin_class.op_list([spin_class.spin_op(s,i)])
        ops3+=[(op, (-1,-1))]

    vecs1=npa_funcs.get_normal_form_vector(ops1, (1,-1), index=0)
    vecs2=npa_funcs.get_normal_form_vector(ops2, (-1,1), index=0)
    vecs3=npa_funcs.get_normal_form_vector(ops3, (-1,-1), index=0)
    
    total_vecs={(1,-1): vecs1, (-1,1):vecs2, (-1,-1): vecs3}
    total_nfs=npa_funcs.get_NF_matrix_elements(total_vecs)
    P = picos.Problem()
    Ms=npa_funcs.get_mom_matrix_M(P, total_vecs, total_nfs)
    shapes={}
    
    for k in Ms.keys():
        shapes[k]=Ms[k].shape
    general_rho([0,2],Ms=Ms,normal_forms=total_nfs)

    return 

def test_tfi_marginal():
    J=4

    h=4
    m=2
    H=mag_ham.TFI(J=J,h=h).reshape(4,4)
    print(H)
    #print(H)
    H_picos=picos.Constant("H", H)
    #mag_ham.XXZ_chain(J=J, Delta=Delta,h=h).reshape(4,4)
    P = picos.Problem()
    rho_list=solve_marginal_problem(P,  m)
    
    P.set_objective("min",(H_picos| rho_list[2]))
    print(P)
    P.solve(solver="mosek")
    print((H| rho_list[2]).np)
    return



def test_xxz_marginal():
    J=1
    Delta=1
    h=0
    m=5
    H=mag_ham.XXZ_chain(J=J,Delta=Delta,h=h).reshape(4,4)
    #print(H)
    H_picos=picos.Constant("H", H)
    #mag_ham.XXZ_chain(J=J, Delta=Delta,h=h).reshape(4,4)
    P = picos.Problem()
    rho_list=solve_marginal_problem(P,  m)
    
    P.set_objective("min",(H_picos| rho_list[2]))
    print(P)
    P.solve(solver="mosek")
    print((H| rho_list[2]).np)
    return

#test_rhos_gen()
#test_tfi_marginal()
#test_xxz_marginal()