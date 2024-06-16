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
        new_ops+=[o.get_translated(j, L)]
    return spin_class.op_list(new_ops)
def generate_all_translations(op,L,class_name, **kwargs):
    # returns list of all translations
    if op.sym=="1":
        return [op]
    all_T=[]
    for i in range(1,L):
        new_state=translate(op, i, L)
        coeff,nf=spin_class.normal_form(new_state, class_name=class_name,**kwargs)
        if np.abs(coeff-1.)>1e-09:
            print("translation error")
        all_T+=[nf]
    return [op]+all_T
        
class block():
    def __init__(self, start, size):
        self.diags=[]
        for i in range(size):
            self.diags+=[(start[0]+i, start[1]+i)]
        return
class mom_symm_block_base:

    def __init__(self, ops,L,P, sec="",class_name=spin_class.spin_op,**kwargs):   
        self.map_TI={} # maps all states that occur to their translation invariant part
        self.map_TI_to_elements={} # kety are a translation invariant state, contains map with each element
        self.G_variables={} # contains picos variables
        self.blocks={} # final blocks
        self.sec=sec
        self.ops=ops
        self.L=L
        self.class_name=class_name
        self.FT=np.zeros((self.L,self.L), dtype=complex)
        self.kwargs=kwargs
        for n in range(L):
            for m in range(L):
                self.FT[n,m]=np.exp(-2j*np.pi*n*m/self.L)#/np.sqrt(self.L)
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
        coeff,nf=spin_class.normal_form(op_new,class_name=self.class_name,**self.kwargs)
        if nf.sym=="1":
            self.map_TI[nf]=nf
            if nf not in self.map_TI_to_elements.keys():
                self.map_TI_to_elements[nf]=[]
            self.map_TI_to_elements[nf]+=[(self.G_variables[g_key][i],coeff, g_key,i)]
                
        else:
            ind_now=nf.ops[0].get_xsite()
            if ind_now!=0:
                op_new_index=translate(nf, self.L-ind_now, self.L)
            else:
                op_new_index=nf
            all_T=generate_all_translations(op_new_index,self.L,class_name=self.class_name,**self.kwargs)
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
        return 
class mom_symm_block_zero(mom_symm_block_base): 

    def __init__(self, ops, L, P, class_name,sec="",**kwargs):    
        super().__init__(ops, L, P, sec=sec,class_name=class_name,**kwargs)
        self.d_vectors=[] # the vector with the values d_1, d_2 ..., d_1=c_1L
        self.d_vector_map={}
        self.generate_blocks(P)
        for operator_key in self.map_TI_to_elements.keys():
            for i in range(0,len(self.map_TI_to_elements[operator_key])-1):
                P.add_constraint(self.map_TI_to_elements[operator_key][i][1]*self.map_TI_to_elements[operator_key][i][0] == \
                                  self.map_TI_to_elements[operator_key][i+1][1]*self.map_TI_to_elements[operator_key][i+1][0])

        if spin_class.op_list([]) in self.map_TI_to_elements.keys():
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
        # can maybe be merged into next loop
        for i in range(len(self.ops)):  
            
            coeff,nf=spin_class.normal_form(self.ops[i],self.class_name, **self.kwargs)
            all_T=generate_all_translations(nf,self.L,self.class_name, **self.kwargs)
            ops=list(set(all_T)& set(self.map_TI.values()))
            if len(ops)==0:
                self.map_TI[nf]=nf
                if nf not in self.map_TI_to_elements.keys():
                    self.map_TI_to_elements[nf]=[]
                self.map_TI_to_elements[nf]+=[(self.d_vectors[i],coeff, "d_{0}".format(i),i)]
            self.d_vector_map[nf]=npa_funcs.nf_monomial(ind=i+1,  prefac=coeff*np.sqrt(self.L))
            P.add_constraint(self.blocks[0][0,i+1]==self.d_vector_map[nf].prefac*self.d_vectors[self.d_vector_map[nf].ind-1,0])
            P.add_constraint(self.blocks[0][i+1,0]==np.conj(self.d_vector_map[nf].prefac)*self.d_vectors[self.d_vector_map[nf].ind-1,0])
            
            
        for i in range(len(self.ops)):  
            for j in range(i,len(self.ops)):
                g_key="G_{2}{0}/{1}".format(i,j,self.sec)
                g_key_2="G_{2}{0}/{1}".format(j,i,self.sec)
                self.G_variables[g_key]=picos.ComplexVariable("G_{2}{0}/{1}".format(i,j,self.sec), (self.L))
                expression_1=[0]*self.L
                expression_2=[0]*self.L
                if j!=i:
                    self.G_variables[g_key_2]=picos.ComplexVariable(g_key_2, (self.L))
              
                for mat_pos in range(0,self.L, 1):
                    if mat_pos==0:
                        shift=1
                    else:
                        shift=0
                    self.generate_G_elements(self.ops[i],self.ops[j], g_key, mat_pos)
                    if i!=j:
                        self.generate_G_elements(self.ops[j],self.ops[i], g_key_2, mat_pos)
                    for pos in range(0,self.L):
                        position_in_G_1=(self.L-pos)%self.L
                        position_in_G_2=(self.L-pos)%self.L
                        #print("pos_1",position_in_G_1,"pos_2",position_in_G_2  )
                        expression_1[mat_pos]+=self.FT[pos,mat_pos]*self.G_variables[g_key][position_in_G_1,0]
                        if i!=j:
                            expression_2[mat_pos]+=self.FT[pos,mat_pos]*self.G_variables[g_key_2][position_in_G_2,0]                       
                    P.add_constraint(self.blocks[mat_pos][shift+i,shift+j]==expression_1[mat_pos])
                    if i!=j:
                        P.add_constraint(self.blocks[mat_pos][shift+j,shift+i]==expression_2[mat_pos])
        return
    
class mom_symm_block(mom_symm_block_base): 
    def __init__(self, ops, L, P,class_name,  sec="",**kwargs):    
        super().__init__(ops, L, P, sec=sec,class_name=class_name,**kwargs)
        self.generate_blocks(P)
        for operator_key in self.map_TI_to_elements.keys():
            for i in range(0,len(self.map_TI_to_elements[operator_key])-1):
                P.add_constraint(self.map_TI_to_elements[operator_key][i][1]*self.map_TI_to_elements[operator_key][i][0] == \
                                  self.map_TI_to_elements[operator_key][i+1][1]*self.map_TI_to_elements[operator_key][i+1][0])

        if spin_class.op_list([]) in self.map_TI_to_elements.keys():
            P.add_constraint(self.map_TI_to_elements[spin_class.op_list([])][0][1]*self.map_TI_to_elements[spin_class.op_list([])][0][0]==1)
    def generate_blocks(self, P):

        for i in range(0,self.L,1):
            self.blocks[i]=picos.HermitianVariable("block_nr{1}{0}".format(i, self.sec), (len(self.ops), len(self.ops)))
            P.add_constraint(self.blocks[i] >> 0)

        for i in range(len(self.ops)):  
            for j in range(i,len(self.ops)):
                g_key="G_{2}{0}/{1}".format(i,j,self.sec)
                g_key_2="G_{2}{0}/{1}".format(j,i,self.sec)
                self.G_variables[g_key]=picos.ComplexVariable("G_{2}{0}/{1}".format(i,j,self.sec), (self.L))
                expression_1=[0]*self.L
                expression_2=[0]*self.L
                if j!=i:
                    self.G_variables[g_key_2]=picos.ComplexVariable(g_key_2, (self.L))
              
                for mat_pos in range(0,self.L, 1):
                    self.generate_G_elements(self.ops[i],self.ops[j], g_key, mat_pos)
                    if i!=j:
                        self.generate_G_elements(self.ops[j],self.ops[i], g_key_2, mat_pos)
                    for pos in range(0,self.L):
                        position_in_G_1=(self.L-pos)%self.L
                        position_in_G_2=(self.L-pos)%self.L
                        #print("pos_1",position_in_G_1,"pos_2",position_in_G_2  )
                        expression_1[mat_pos]+=self.FT[pos,mat_pos]*self.G_variables[g_key][position_in_G_1,0]
                        if i!=j:
                            expression_2[mat_pos]+=self.FT[pos,mat_pos]*self.G_variables[g_key_2][position_in_G_2,0]                       
                    P.add_constraint(self.blocks[mat_pos][i,j]==expression_1[mat_pos])
                    if i!=j:
                        P.add_constraint(self.blocks[mat_pos][j,i]==expression_2[mat_pos])
        return
class momentum_basis:
    def __init__(self, operators, L, P,class_name,**kwargs):
        self.sectors={}
        self.L=L
        self.class_name=class_name
        self.kwargs=kwargs

        for i in range(0, len(operators.keys())):
            key=list(operators.keys())[i]
            if i==0:
                M=mom_symm_block_zero(ops=operators[key], L=self.L,P=P, sec="{0}_".format(i), class_name=self.class_name,**self.kwargs)
            else:
                M=mom_symm_block(ops=operators[key], L=self.L,P=P, sec="{0}_".format(i), class_name=self.class_name,**self.kwargs)
            self.sectors[key]=M
            
        # relate elements
        for i in range(0,len(self.sectors.keys())):
            key_i=list(self.sectors.keys())[i]
            for j in range(i+1,len(self.sectors.keys())):
                for key_1 in self.sectors[key_i].map_TI_to_elements.keys():
                    key_j=list(self.sectors.keys())[j]
                    all_T=generate_all_translations(key_1,L,class_name=self.class_name,**self.kwargs)
                    ops=list(set(all_T)& set(self.sectors[key_j].map_TI_to_elements.keys()))
                    #print("jjeje", ops)
                    if len(ops)!=0:
                        #print(tot[j].map_TI_to_elements[ops[0]])
                        P.add_constraint(self.sectors[key_j].map_TI_to_elements[ops[0]][0][0]*np.conj(self.sectors[key_j].map_TI_to_elements[ops[0]][0][1])\
                                          ==self.sectors[key_i].map_TI_to_elements[key_1][0][0]*np.conj(self.sectors[key_i].map_TI_to_elements[key_1][0][1]))
                        
    
            
            