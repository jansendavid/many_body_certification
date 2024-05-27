#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 17:28:26 2024

@author: david
"""

import tensor as tensor
import tensornetwork as tn
import numpy as np
import pickle
from marginal_problem_hamiltonians import *
from marginal_problem import *
def random_isometry(D1, D2, dtype):
  """Generate a random isometric matrix of dimension D1 x D2.

  We require D1 <= D2.

  Args:
    D1: Left dimension.
    D2: Right dimension.
    dtype: Element type.

  Returns:
    V: An isometry.
  """
  if not D1 <= D2:
    raise ValueError("The left dimension must be <= the right dimension.")
  A = np.random.rand(D2, D1)
  Q, R = np.linalg.qr(A)
  r = np.diagonal(R)
  L = np.diag(r / np.abs(r))
  return np.transpose(Q @ L, (1,0))
def test_vectorization():
    A=np.random.rand(3,2,2,3)+1j*np.random.rand(3,2,2,3)
    rho=np.random.rand(8,8)
    rho=rho+np.transpose(rho)
    rho=rho.reshape(2,2,2,2,2,2)
    rho_fin=tn.ncon([A, rho, np.conj(A)], [[-1,1,2,-2],[1,2,-3,3,4,-6], [-4,3,4,-5]])
    Add=A.transpose(0,3,1,2)
    Add=Add.reshape(9,4)
    Add=np.kron(Add, np.eye(2)) 
    Add.shape
    rho_fin_2=Add@rho.reshape(8,8)@np.conj(np.transpose(Add))
    np.allclose(rho_fin_2.reshape(3,3,2,3,3,2),rho_fin)
    rho_3=tensor.apply_vetorization(rho, Add)
    rho_3=rho_3.reshape(3,3,2,3,3,2)
    np.allclose(rho_fin_2.reshape(3,3,2,3,3,2),rho_fin),np.allclose(rho_fin_2.reshape(3,3,2,3,3,2),rho_3)
    return
#test_vectorization()
def test_apply_ttn():
    iso=random_isometry(3, 4, float)
    iso=iso.reshape(3,2,2)
    print(iso.shape)
    rho=np.random.randn(4,4)
    rho=rho+np.transpose(rho)
    rho/=np.trace(rho)
    print(np.trace(rho))
    rho_new=tensor.apply_CPTPmap(rho, iso)
    print(np.allclose(np.trace(rho_new), np.trace(rho)))
#test_apply_ttn() 
def test_marginal_prob_ttn():
    D=3
    J=4
    h=4
    itr=800
    ttn_dir_name="Tensors/DATA_TTN/TFI/TTN_TFI_J{0:.2f}h{1:.2f}D{2}itr{3}DIR".format(J,h,D, itr)   
    with open(ttn_dir_name+'/isos.p', 'rb') as handle:
        isos = pickle.load(handle)
    for t in isos:
        print(t.shape)
    H=TFI(J, h).reshape(4,4)
    Eye_2=np.eye(4,4)
    H_tot=np.kron(Eye_2, H)+np.kron( H, Eye_2)
    m=6
    P = picos.Problem()
    H_picos=picos.Constant("H", H_tot)
    rho_list=solve_marginal_problem_TTN(P=P,  m=m,isos=isos, name="")
    P.set_objective("min",(H_picos| rho_list[4]))

    print(P)
    P.solve(solver="mosek")
    print((H_picos| rho_list[4]).np/2)

#test_marginal_prob_ttn()


def test_marginal_prob_ttn():
    D=3
    J=4
    h=4
    itr=800
    ttn_dir_name="Tensors/DATA_TTN/TFI/TTN_TFI_J{0:.2f}h{1:.2f}D{2}itr{3}DIR".format(J,h,D, itr)   
    with open(ttn_dir_name+'/isos.p', 'rb') as handle:
        isos = pickle.load(handle)
    for t in isos:
        print(t.shape)
    H=TFI(J, h).reshape(4,4)
    Eye_2=np.eye(4,4)
    H_tot=np.kron(Eye_2, H)+np.kron( H, Eye_2)
    m=6
    sx=np.array([[0,1],[1,0]])
    O=np.kron(sx, np.eye(2))
    O=np.kron(O, np.eye(2))
    O=np.kron(O, sx)
    P = picos.Problem()
    H_picos=picos.Constant("H", H_tot)
    O_picos=picos.Constant("H",O)
    rho_list=solve_marginal_problem_TTN(P=P,  m=m,isos=isos, name="")
    P.add_constraint((H_picos| rho_list[4])/2<=-1.27)
    P.add_constraint((H_picos| rho_list[4])/2>=-1.30)
    P.set_objective("min",(O_picos| rho_list[4]))

    print(P)
    P.solve(solver="mosek")
    print((H_picos| rho_list[4]).np/2)
    print((O_picos| rho_list[4]))
    
test_marginal_prob_ttn()