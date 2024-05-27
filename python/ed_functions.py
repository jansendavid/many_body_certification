#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 10:37:46 2023

@author: david
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

"""
To do generate sparse Hamiltonians, using Lanczos to only get the groundstate
"""

import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse import lil_matrix
def make_spinbasis(L, conserve_M=False, M=0):
    """
    Function to generate the spin basis

    Parameters
    ----------
    L : int
        Length of chain.
    conserve_M : (bool), optional
        DESCRIPTION. The default is False. Set True if magnetization 
        is conserved.
    M : int, optional
        DESCRIPTION. The default is 0. Magnetization. Only matters if 
        conserve_M=True
    Returns
    -------
    states : (dict)
        Dictionary containing the spin basis states as values and the decimal 
        representation.

    """
    
    
    states={}
    # the position of the state in the basis
    state_nr=0
    for i in np.arange(0, 2**L,1):
        state=spin_state(L, i)
        if(conserve_M):
            if(state.mag==M):
                state.set_basis_pos(state_nr)
                states[state.id]=state
                state_nr+=1
        else:
            state.set_basis_pos(state_nr)
            states[state.id]=state
            state_nr+=1
    return states

def convert_array_to_number(arr):
    """
    This function converts an array of spins e.g. [0,1,1,0]
    into their corresponding integer representation (each array is encoded as the binary representation
                                                     of an integer)
    Parameters
    ----------
    arr : numpy array
    Takes a numpy array with the elements 0,1 

    Returns
    -------
    nr : int
        The integer from which the arr is the binary representation

    """
    nr=0
    j=0
    for i in np.flip(arr):
        nr+=i*2**(j)
        j+=1
    return nr

########### Observables #############################
def generate_XXZ(states,J,Delta,h,L, PB=True, use_sparse=False):
    """
    This function generates the XXZ-Hamiltonian with periodic boundary condictions.
    
    H=\sum\limits_{i=0}^{L-1}( J (S_i^{x}S^{x}_{i+1} +S_i^{y}S^{y}_{i+1})+
                              \Delta (S_i^{z}\S^{z}_{i+1})+hS^{z}_i)

    Parameters
    ----------
    states : TYPE
        DESCRIPTION.
    J : float
        The spin-spin interaction.
    h : float
        Magnetic field.
    L : int
        Length of chain.
    PB : bool
        Indicating whether to use periodic boundary conditions

    Returns
    -------
    H : numpy array
        Hamiltonian of the system.

    """
    end_val=L

    if PB==False:
        end_val=L-1
    if use_sparse:
        H=lil_matrix((len(states),len(states)))
    else:
        H=np.zeros((len(states),len(states)))       
    for state in states.values():
        for i in np.arange(0,L,1):
            H[state.basis_pos,state.basis_pos]+=(float(-state.state[i])+0.5)*h
        for i in np.arange(0,end_val,1):
            j=(i+1)%L
            
            if(state.state[i]==state.state[j]):
                #print(None)
                H[state.basis_pos,state.basis_pos]+=Delta*1./4
            else:
                H[state.basis_pos,state.basis_pos]+=-Delta*1./4
                arr=np.copy(state.state)
                arr_cp=np.copy(state.state)
                arr[i]=arr_cp[j]
                arr[j]=arr_cp[i]
                H[state.basis_pos,states[convert_array_to_number(arr)].basis_pos]+=J*1./2
                        
                
    return H


def generate_fermi_ekin(states,t0,L, PB=True, use_sparse=False):
    """
    This function generates the XXZ-Hamiltonian with periodic boundary condictions.
    
    H=\sum\limits_{i=0}^{L-1}( J (S_i^{x}S^{x}_{i+1} +S_i^{y}S^{y}_{i+1})+
                              \Delta (S_i^{z}\S^{z}_{i+1})+hS^{z}_i)

    Parameters
    ----------
    states : TYPE
        DESCRIPTION.
    J : float
        The spin-spin interaction.
    h : float
        Magnetic field.
    L : int
        Length of chain.
    PB : bool
        Indicating whether to use periodic boundary conditions

    Returns
    -------
    H : numpy array
        Hamiltonian of the system.

    """
    end_val=L
    if PB==False:
        end_val=L-1
    if use_sparse:
        H=lil_matrix((len(states),len(states)))
    else:
        H=np.zeros((len(states),len(states)))        
    for state in states.values():
        for i in np.arange(0,end_val,1):
            j=(i+1)%L
            if(state.state[i]!=state.state[j]):
                arr=np.copy(state.state)
                arr_cp=np.copy(state.state)
                arr[i]=arr_cp[j]
                arr[j]=arr_cp[i]
                sign=1.
                nr_fermions=np.sum(state.state)
                if j==0 and nr_fermions%2==0:
                    sign=-1
                    
                H[states[convert_array_to_number(arr)].basis_pos,state.basis_pos]+=-t0*sign
                        
                
    return H
def generate_fermi_int(states,U,L, PB=True, use_sparse=False):
    """
    This function generates the XXZ-Hamiltonian with periodic boundary condictions.
    
    H=\sum\limits_{i=0}^{L-1}( J (S_i^{x}S^{x}_{i+1} +S_i^{y}S^{y}_{i+1})+
                              \Delta (S_i^{z}\S^{z}_{i+1})+hS^{z}_i)

    Parameters
    ----------
    states : TYPE
        DESCRIPTION.
    J : float
        The spin-spin interaction.
    h : float
        Magnetic field.
    L : int
        Length of chain.
    PB : bool
        Indicating whether to use periodic boundary conditions

    Returns
    -------
    H : numpy array
        Hamiltonian of the system.

    """
    end_val=L
    if PB==False:
        end_val=L-1
    if use_sparse:
        H=lil_matrix((len(states),len(states)))
    else:
        H=np.zeros((len(states),len(states)))        
    for state in states.values():
        for i in np.arange(0,end_val,1):
            j=(i+1)%L
            if(state.state[i]==1 and state.state[j]==1):
                H[state.basis_pos,state.basis_pos]+=U
                        
                
    return H
########### Observables #############################
def generate_TFI(states,J,h,L, PB=True, use_sparse=False):
    """
    This function generates the Hamiltonian with periodic boundary condictions.
    
    H=\sum\limits_{i=0}^{L-1}( J (sigma_i^{z}S^{z}_{i+1}  +\frac{h}{2}sigma^{x}_i)

    Parameters
    ----------
    states : TYPE
        DESCRIPTION.
    J : float
        The spin-spin interaction.
    h : float
        Magnetic field.
    L : int
        Length of chain.
    PB : bool
        Indicating whether to use periodic boundary conditions

    Returns
    -------
    H : (numpy array)
        Hamiltonian of the system.

    """
    end_val=L
    if PB==False:
        end_val=L-1
    if use_sparse:
        H=lil_matrix((len(states),len(states)))
    else:
        H=np.zeros((len(states),len(states))) 

    for state in states.values():
        for i in np.arange(0,end_val,1):                       
            arr=np.copy(state.state)
            arr[i]=(arr[i]+1)%2
            H[states[convert_array_to_number(arr)].basis_pos,state.basis_pos]+=h/4
            
            j=(i+1)%L
            if(state.state[i]==state.state[j]):
                #print(None)
                H[state.basis_pos,state.basis_pos]+=J/4
            else:
                H[state.basis_pos,state.basis_pos]+=-J/4
        if(PB==False):
            arr=np.copy(state.state)
            arr[L-1]=(arr[L-1]+1)%2
            H[states[convert_array_to_number(arr)].basis_pos,state.basis_pos]+=h/4
            
                
                        
                
    return H


# def generate_TFI_forsdp(states,J,h,L, PB=True, use_sparse=False):
#     """
#     This function generates the Hamiltonian with periodic boundary condictions.
    
#     H=\sum\limits_{i=0}^{L-1}( J (sigma_i^{z}S^{z}_{i+1}  +\frac{h}{2}sigma^{x}_i)

#     Parameters
#     ----------
#     states : TYPE
#         DESCRIPTION.
#     J : float
#         The spin-spin interaction.
#     h : float
#         Magnetic field.
#     L : int
#         Length of chain.
#     PB : bool
#         Indicating whether to use periodic boundary conditions

#     Returns
#     -------
#     H : (numpy array)
#         Hamiltonian of the system.

#     """
#     end_val=L
#     if PB==False:
#         end_val=L-1
#     if use_sparse:
#         H=lil_matrix((len(states),len(states)))
#     else:
#         H=np.zeros((len(states),len(states))) 

#     for state in states.values():

#         for i in np.arange(0,end_val,1):                       
#             arr=np.copy(state.state)
#             arr[i]=(arr[i]+1)%2
#             if i ==0:
#                 H[states[convert_array_to_number(arr)].basis_pos,state.basis_pos]+=0.5*h/4
#             elif i==end_val-1:
#                 H[states[convert_array_to_number(arr)].basis_pos,state.basis_pos]+=0.5*h/4
#             else:
#                 H[states[convert_array_to_number(arr)].basis_pos,state.basis_pos]+=h/4
            
#             j=(i+1)%L
#             if(state.state[i]==state.state[j]):
#                 #print(None)
#                 H[state.basis_pos,state.basis_pos]+=J/4
#             else:
#                 H[state.basis_pos,state.basis_pos]+=-J/4
#         if(PB==False):
#             arr=np.copy(state.state)
#             arr[L-1]=(arr[L-1]+1)%2
#             H[states[convert_array_to_number(arr)].basis_pos,state.basis_pos]+=0.5*h/4
            
                
                        
                
#     return H




def generate_mag(states,L,use_sparse=False):
    """
    Genereates the magnetization 
    \sum\limits_{i=0}^{L-1}S^{z}_i
    Parameters
    ----------
    states : dict{state id, state}
    Basis of the system
        .
    L : int
        Length of the system.

    Returns
    -------
    M : numpy array
        Matrix containing the magnetization.

    """
    if use_sparse:
        M=lil_matrix((len(states),len(states)))
    else:
        M=np.zeros((len(states),len(states)))
    for state in states.values():
        for i in np.arange(0,L,1):

            M[state.basis_pos,state.basis_pos]+=(float(-state.state[i])+0.5)

                        
                
    return M

def generate_density(states,i,use_sparse=False):
    """
    Genereates the magnetization 
    \sum\limits_{i=0}^{L-1}S^{z}_i
    Parameters
    ----------
    states : dict{state id, state}
    Basis of the system
        .
    L : int
        Length of the system.

    Returns
    -------
    M : numpy array
        Matrix containing the magnetization.

    """
    if use_sparse:
        M=lil_matrix((len(states),len(states)))
    else:
        M=np.zeros((len(states),len(states)))
    for state in states.values():
        M[state.basis_pos,state.basis_pos]+=(state.state[i])

                        
                
    return M



def generate_Szcorrel(states,L,i,j, use_sparse=False):
    """
    Genereates the S^z_i S^z_{i+r} correlation function with periodic boundary conditions.
    \sum\limits_{i=0}^{L-1}S^{z}_iS^{z}_{i+r}
    Parameters
    ----------
    states : dict{state id, state}
    Basis of the system.
        
    L : int
        Length of the system.
    r : int 
        The distance between the operators

    Returns
    -------
    M : numpy array
        Matrix containing the correlation function.

    """
    if use_sparse:
        M=lil_matrix((len(states),len(states)))
    else:
        M=np.zeros((len(states),len(states))) 
    for state in states.values():
            M[state.basis_pos,state.basis_pos]+=(float(-state.state[i])+0.5)*(float(-state.state[j])+0.5)

                        
                
    return M
def generate_Sxcorrel(states,L,i,j,use_sparse=False):
    """
    Genereates the S^z_i S^z_{i+r} correlation function with periodic boundary conditions.
    \sum\limits_{i=0}^{L-1}S^{z}_iS^{z}_{i+r}
    Parameters
    ----------
    states : dict{state id, state}
    Basis of the system.
        
    L : int
        Length of the system.
    r : int 
        The distance between the operators

    Returns
    -------
    M : numpy array
        Matrix containing the correlation function.

    """
    if use_sparse:
        M=lil_matrix((len(states),len(states)))
    else:
        M=np.zeros((len(states),len(states)))   
    for state in states.values():
        if(state.state[i]!=state.state[j]):
              arr=np.copy(state.state)
              arr_cp=np.copy(state.state)
              arr[i]=(arr_cp[i]+1)%2
              arr[j]=(arr_cp[j]+1)%2
              M[states[convert_array_to_number(arr)].basis_pos,state.basis_pos]+=0.25
                        
                
    return M




########### classes ####################
class spin_state():
    """
    This is the spin state class
    """
    def __init__(self, L,i):
        """

        Parameters
        ----------
        L : int 
            Length of chain.
        i : int
            The number that shall be encoded into binary representation 
            as a spin state

        Returns
        -------
        None.

        """
        self.state=np.zeros(L, dtype=int)
        self.maxval=0
        self.id=i
        self.basis_pos=-1

        for j in np.arange(0,L,1):
            self.maxval+=2**(j)
        self.nr=i
        if(i>self.maxval):
            print("i to large")
            return
        i_cp=i
        index=0
        while i_cp>0:
            if(i_cp%2!=0):
                self.state[L-1-index]=1
            i_cp//=2
            index+=1
        self.mag=np.sum(self.state)
        return

    def set_basis_pos(self, x):
        """
        

        Parameters
        ----------
        x : int
            The position of the state in a basis.

        Returns
        -------
        None.

        """
        self.basis_pos = x 
        return
    
    
    
