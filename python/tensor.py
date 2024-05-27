#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 17:15:45 2024

@author: david
"""

import tensornetwork as tn
import numpy as np
from scipy.linalg import null_space

def apply_vetorization(rho, A):
    rho_new=rho.reshape(-1,1)
    op=(np.kron(A, np.conj(A)))
    return op@rho_new


def apply_CPTPmap(rho, iso):
    # expecting tree tensor of shape (D,d,d)
    D=iso.shape[0]
    d=iso.shape[1]
    iso_extended=tn.ncon([iso,np.eye(D,D+1)],[[1,-2,-3],[1,-1]])
    iso_extended=iso_extended.reshape(D+1,d*d)
    isos_ns=null_space(iso_extended.reshape(D+1,d*d))
    gvec=np.zeros(D+1)
    gvec[-1]=1
    kraus_ops=[]
    for n in np.arange(0,isos_ns.shape[1],1):
        kraus_ops+=[tn.ncon([gvec.reshape(-1,1),np.transpose(np.conj(isos_ns[:,n])).reshape(1,-1) ], [[-1,1],[1,-2]])]
    # apply to rho
    print("here")
    OP=(np.kron(iso_extended, np.conj(iso_extended)))
    for k_op in kraus_ops:
        OP+=(np.kron(k_op, np.conj(k_op)))
    rho_new=rho.reshaped((OP.shape[1],1),order="C")
            
    rho_new=OP*rho_new # picos mat amt multiplication

    rho_new=rho_new.reshaped((D+1,D+1),order="C")
        
        
    return rho_new