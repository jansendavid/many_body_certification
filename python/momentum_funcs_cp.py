#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 14:08:21 2024

@author: david
"""
import numpy as np
from copy import deepcopy
import spin_class as spin_class
import npa_funcs as npa_funcs
from npa_funcs import sym
import picos
import npa_hams as npa_hams
import sys
import networkx as nx
import npa_funcs as npa_funcs

# new idea in G_matrices is object of length L were 
# we store ([operator], coefficient), the index will be the postion
# convention, I store the operator and the coefficient in front
class nf_monomial_op:
    """
    Class containing the monomial operator position in the vector u and and their prefactor
    """
    def __init__(self, op,   prefac):
        self.op=op
        self.prefac=prefac

def translate(op, j, L):
    new_ops=[]
    for o in op.ops:
        new_ops+=[spin_class.spin_op(o.cor,(o.site+j)%L)]
    return spin_class.op_list(new_ops)
def generate_all_translations(op,L):
    # returns list of all translations
    if op.sym=="1":
        return [op]
    all_T=[]
    for i in range(1,L):
        new_state=translate(op, i, L)
        if new_state in all_T:
            break 
        else:
            all_T+=[new_state]
    return [op]+all_T
        
class block():
    def __init__(self, start, size):
        self.diags=[]
        for i in range(size):
            self.diags+=[(start[0]+i, start[1]+i)]
        return
class mom_symm_block_base:
    ops=[]    
    map_TI={} # maps all states that occur to their translation invariant part
    map_TI_to_elements={} # kety are a translation invariant state, contains map with each element
    
    G_matrices={} # contains vector [1,s^x_1s_2s^x_2s_3,s^x_1s_2s^x_3s_4,...]
    G_variables={} # contains picos variables 
    D_matrices={}
    G_representatives={} # contains on element representative wt TI
    blocks={} # final blocks
    def __init__(self, ops,L,P, sec=""):
        self.sec=sec
        self.ops=ops
        self.L=L
        self.FT=np.zeros((self.L,self.L), dtype=complex)
        for n in range(L):
            for m in range(L):
                self.FT[n,m]=np.exp(-2j*np.pi*n*m/self.L)/np.sqrt(self.L)
        # this can be optimized to only store one vector
        #for i in self.map_TI.keys():
         #   print(i.sym , "-> ", self.map_TI[i].sym)

             

        return
    def generate_G_elements(self, op1,op2, g_key, i):
        if i>0:
            op2_new=translate(op2, i, self.L)
        else:
            op2_new=op2
        op_new=spin_class.op_list(op1.ops+op2_new.ops)
        coeff,nf=spin_class.normal_form(op_new)
        if nf.sym=="1":
            self.map_TI[nf]=nf
            if nf not in self.map_TI_to_elements.keys():
                self.map_TI_to_elements[nf]=[]
            self.map_TI_to_elements[nf]+=[(self.G_variables[g_key][i],coeff, g_key,i)]
                
        else:
            ind_now=nf.ops[0].site
            if ind_now!=0:
                op_new_index=translate(nf, self.L-ind_now, self.L)
            else:
                op_new_index=nf
            all_T=generate_all_translations(op_new_index,self.L)
            ops=list(set(all_T)& set(self.map_TI.values()))
            if len(ops)==0:
                self.map_TI[nf]=nf
                if nf not in self.map_TI_to_elements.keys():
                    self.map_TI_to_elements[nf]=[]
                self.map_TI_to_elements[nf]+=[(self.G_variables[g_key][i],coeff, g_key,i)]
            else:
                self.map_TI[nf]=ops[0]
                if ops[0] not in self.map_TI_to_elements.keys():
                    self.map_TI_to_elements[ops[0]]=[]
                self.map_TI_to_elements[ops[0]]+=[(self.G_variables[g_key][i],coeff, g_key,i)]
            if len(ops)>1:
                print("errpr")
            #         vec_rep[op_new_index]=i
        return 
class mom_symm_block_zero(mom_symm_block_base): 
    d_vectors=[] # the vector with the values d_1, d_2 ..., d_1=c_1L
    d_vector_map={}

    def __init__(self, ops, L, P, sec=""):    
        super().__init__(ops, L, P, sec=sec)
        self.generate_blocks(P)
        for operator_key in self.map_TI_to_elements.keys():
            for i in range(0,len(self.map_TI_to_elements[operator_key])-1):
                P.add_constraint(self.map_TI_to_elements[operator_key][i][1]*self.map_TI_to_elements[operator_key][i][0] == \
                                  self.map_TI_to_elements[operator_key][i+1][1]*self.map_TI_to_elements[operator_key][i+1][0])
        for key in self.d_vector_map.keys():
            # maybe unnecisary
            all_T=generate_all_translations(key,self.L)
            ops=list(set(all_T)& set(self.map_TI_to_elements.keys()))
            if len(ops)>0:
                if len(ops)>1:
                    print("error")                
                P.add_constraint(self.map_TI_to_elements[ops[0]][0][1]*self.map_TI_to_elements[operator_key][0][0] == \
                                  self.d_vector_map[key].prefac*self.d_vectors[self.d_vector_map[key].ind-1,0])
        if spin_class.op_list([]) in self.map_TI_to_elements.keys():
            print("PPP", self.map_TI_to_elements[spin_class.op_list([])][0][0])
            P.add_constraint(self.map_TI_to_elements[spin_class.op_list([])][0][1]*self.map_TI_to_elements[spin_class.op_list([])][0][0]==1)
                
                
            
                
                
             
        
    def generate_blocks(self, P):
        self.blocks[0]=picos.HermitianVariable("block_nr{1}{0}".format(0, self.sec), (len(self.ops)+1, len(self.ops)+1))
        P.add_constraint(self.blocks[0] >> 0)
        P.add_constraint(self.blocks[0][0,0]==1)
        for i in range(1,self.L,1):
            self.blocks[i]=picos.HermitianVariable("block_nr{1}{0}".format(i, self.sec), (len(self.ops), len(self.ops)))
            P.add_constraint(self.blocks[i] >> 0)

    
        shift=0 # shift for the zeroth block due to zero vectors
        self.d_vectors=picos.ComplexVariable("d_vec", (len(self.ops),1))
        for i in range(len(self.ops)):  
            coeff,nf=spin_class.normal_form(self.ops[i])
            self.d_vector_map[nf]=npa_funcs.nf_monomial(ind=i+1,  prefac=coeff*np.sqrt(self.L))
            P.add_constraint(self.blocks[0][0,i+1]==self.d_vector_map[nf].prefac*self.d_vectors[self.d_vector_map[nf].ind-1,0])
            P.add_constraint(self.blocks[0][i+1,0]==np.conj(self.d_vector_map[nf].prefac)*self.d_vectors[self.d_vector_map[nf].ind-1,0])
            for j in range(i,len(self.ops)):
                g_key="G_{2}{0}/{1}".format(i,j,self.sec)
                self.G_variables[g_key]=picos.ComplexVariable("G_{2}{0}/{1}".format(i,j,self.sec), (self.L))
                expression_1=[0]*self.L
                expression_2=[0]*self.L
                for mat_pos in range(0,self.L, 1):
                    if mat_pos==0:
                        shift=1
                    else:
                        shift=0
                    self.generate_G_elements(self.ops[i],self.ops[j], g_key, mat_pos)
                    for pos in range(0,self.L):
                        position_in_G=(self.L-pos)%self.L
                        expression_1[mat_pos]+=self.FT[pos,mat_pos]*self.G_variables[g_key][position_in_G,0]
                        expression_2[mat_pos]+=np.conj(self.FT[pos,mat_pos])*self.G_variables[g_key][position_in_G,0]                        
                    P.add_constraint(self.blocks[mat_pos][shift+i,shift+j]==expression_1[mat_pos])
                    if i!=j:
                        P.add_constraint(self.blocks[mat_pos][shift+j,shift+i]==expression_2[mat_pos])
        return
