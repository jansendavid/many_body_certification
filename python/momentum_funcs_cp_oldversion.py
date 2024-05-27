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


def translate(op, j, L):
    new_ops=[]
    for o in op.ops:
        new_ops+=[spin_class.spin_op(o.cor,(o.site+j)%L)]
    return spin_class.op_list(new_ops)
def generate_all_translations(op,L):
    # returns list of all translations
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
    G_matrices={} # contains vector [1,s^x_1s_2s^x_2s_3,s^x_1s_2s^x_3s_4,...]
    G_variables={} # contains picos variables 
    D_matrices={}
    blocks={} # final blocks
    ind1=1 # how the G blocks will be initializes
    L=0
    sec=""# extra label for symmetrie sector
    def __init__(self, ops,L,P, sec=""):
        self.sec=sec
        self.ops=ops
        self.L=L
        self.FT=np.zeros((self.L,self.L), dtype=complex)
        for n in range(L):
            for m in range(L):
                self.FT[n,m]=np.exp(-2j*np.pi*n*m/self.L)/np.sqrt(self.L)
        #self.generate_G_matrices()
        #self.generate_D_matrices(P)
        return
    def fix_ones_inGmatrices(self,P):
        # # fix all the ones
        for key_g in  self.G_matrices.keys():
            vec=self.G_matrices[key_g]
            for n in vec.keys():
                if n.sym=="1":
                    G_element_1=self.G_matrices[key_g][n]
                    P.add_constraint(np.conj(G_element_1.prefac)*self.G_variables[key_g][G_element_1.ind]==1)
        
        return
    def generate_G_matrices(self):
        i1=self.ind1
        # Generate G matrix elements
        for o1 in self.ops:  
            i2=self.ind1
            for o2 in self.ops:
                #print("G_{0}/{1}".format(i1,i2),"elements",o1.sym, o2.sym)
                vec=self.generate_G_elements( o1,o2)
                #print("G_{0}{1}".format(i1,i2))
                #for n in vec.keys():
                #    print(n.sym, vec[n].get())
                print("G_{2}{0}/{1}".format(i1,i2,self.sec))
                self.G_matrices["G_{2}{0}/{1}".format(i1,i2,self.sec)]=vec
                self.G_variables["G_{2}{0}/{1}".format(i1,i2,self.sec)]=picos.ComplexVariable("G_{2}{0}/{1}".format(i1,i2,self.sec), (len(vec)))
                i2+=1
            i1+=1 
        return
    def associate(self,P):
        print("hehe")
        for i in range(1,len(self.ops)+1):
            for j in range(i,len(self.ops)+1): 
                for l in self.G_matrices["G_{2}{0}/{1}".format(i,j,self.sec)].keys():
                    G_sec_monomial=self.G_matrices["G_{2}{0}/{1}".format(i,j,self.sec)][l]
                    print("G_{2}{0}/{1}".format(i,j,self.sec), "l", l.sym)
                    print("end")
                    P.add_constraint(G_sec_monomial.prefac*self.G_variables["G_{2}{0}/{1}".format(i,j,self.sec)][G_sec_monomial.ind,0]\
                                      ==np.conj(G_sec_monomial.prefac)*self.G_variables["G_{2}{0}/{1}".format(j,i,self.sec)][G_sec_monomial.ind,0])
        return
            
    def generate_D_matrices(self,P):
        # # Generate D matrices, diagonal where each entry is sum of G matrices elements
        # TODO, optimize, not keepng L^2 terms
        #print("START")
        i1=self.ind1
        for o1 in self.ops:  
            i2=self.ind1
            for o2 in self.ops:
                self.D_matrices["D_{2}{0}/{1}".format(i1,i2, self.sec)]=picos.ComplexVariable("D_{2}{0}/{1}".format(i1,i2, self.sec), (self.L))
                #print("D_{0}{1}".format(i1,i2))
                for mat_pos in range(0,self.L, 1):
                    shift=0
                    expression=0
                    #print("MAT", mat_pos)
                    for pos_1 in range(0,self.L, 1):
                        for pos_2 in range(0,self.L, 1):
                            coeff,el=self.generate_D_elements(o1,o2, pos_2+shift)
                            #print(el.sym, "G_{0}/{1}".format(i1,i2)," at ", pos_1, pos_2)
                            G_sec_monomial=self.G_matrices["G_{2}{0}/{1}".format(i1,i2,self.sec)][el]
                            # multiply G_sec_monomial.prefac*? I GUESS NOT
                            expression+=np.transpose(np.conj(self.FT))[mat_pos,pos_1]*self.G_variables["G_{2}{0}/{1}".format(i1,i2,self.sec)][G_sec_monomial.ind,0]*self.FT[pos_2, mat_pos]
                        shift=self.L-1-pos_1
                    P.add_constraint(self.D_matrices["D_{2}{0}/{1}".format(i1,i2, self.sec)][mat_pos]==expression)
                i2+=1
            i1+=1
        return
    
    def associate_Gelements(self,P,already_associated_elements):
        # associates all elements in the G matrices with each other
        for i in range(len(self.G_matrices.keys())):
            keys_i=list(self.G_matrices.keys())[i]
            operators_i=list(self.G_matrices[keys_i].keys())
            for j in range(i+1,len(self.G_matrices.keys())):
                keys_j=list(self.G_matrices.keys())[j]
                operators_j=list(self.G_matrices[keys_j].keys())
                
                
                overlap=list(set(operators_i) & set(operators_j))
                for o in overlap:
                #    print( o)
                    # If two sectors contain the same elemet, we relate them to each other
                    #print("overlap", first_sec, second_sec,o.sym, (first_set[o].ind1,first_set[o].ind2),(sec_set[o].ind1,sec_set[o].ind2))
                    #if (keys_j,o) not in already_associated_elements:
                        G_element_1=self.G_matrices[keys_i][o]
                        G_element_2=self.G_matrices[keys_j][o]
                        #
                        P.add_constraint(np.conj(G_element_1.prefac)*self.G_variables[keys_i][G_element_1.ind]==np.conj(G_element_2.prefac)*self.G_variables[keys_j][G_element_2.ind])
                        already_associated_elements+=[(keys_j,o)]
            return
        # old function, problem distinguishing s^z_0s^x_2 and s^x_2s^z_0
    # def generate_G_elements(self, op1,op2):
    #     vec={}
    #     ind=0 # dummy index to ensure position of element in G_{i,j} array
    #     #print("rnu")
    #     for i in range(0,self.L, 1):
    #         if i>0:
    #             op2_new=translate(op2, i, self.L)
    #         else:
    #             op2_new=op2
    #         op_new=spin_class.op_list(op1.ops+op2_new.ops)
    #         coeff,nf=spin_class.normal_form(op_new)
    #         if nf.sym!="1":
    #             ind_now=nf.ops[0].site
    #             if ind_now!=0:
    #                 op_new_index=translate(nf, self.L-ind_now, self.L)
    #             else:
    #                 op_new_index=nf
    #         else:
    #             op_new_index=nf
    #         if op_new_index not in vec.keys():
    #             vec[op_new_index]=npa_funcs.nf_monomial(ind=ind,  prefac=coeff)
    #             ind+=1
    #     return vec
    def generate_G_elements(self, op1,op2):
        vec={}
        ind=0 # dummy index to ensure position of element in G_{i,j} array
        #print("rnu")
        for i in range(0,self.L, 1):
            if i>0:
                op2_new=translate(op2, i, self.L)
            else:
                op2_new=op2
            op_new=spin_class.op_list(op1.ops+op2_new.ops)
            coeff,nf=spin_class.normal_form(op_new)
            if nf.sym=="1":
                if nf not in vec.keys():
                    vec[nf]=npa_funcs.nf_monomial(ind=ind,  prefac=coeff)
                    ind+=1
            else:
                ind_now=nf.ops[0].site
                if ind_now!=0:
                    op_new_index=translate(nf, self.L-ind_now, self.L)
                else:
                    op_new_index=nf
                # overview
                print("overview", op_new_index.sym)
                
                all_translations=generate_all_translations(op_new_index,self.L)
                for a in all_translations:
                    print(a.sym)
                print("end")
                overlap=list(set(all_translations) & set(vec.keys()))
                if len(overlap)==0:
                    vec[op_new_index]=npa_funcs.nf_monomial(ind=ind,  prefac=coeff)
                    ind+=1
                        
        return vec


    
    def generate_D_elements(self,op1,op2, translation):
        if translation>0:
            op2_new=translate(op2, translation, self.L)
        else:
            op2_new=op2
        op_new=spin_class.op_list(op1.ops+op2_new.ops)
        coeff,nf=spin_class.normal_form(op_new)
        # translate to first site having index 0
        op_new_index=nf
        if nf.sym!="1":
            ind_now=nf.ops[0].site
            if ind_now!=0:
                op_new_index=translate(nf, self.L-ind_now, self.L)
        return (coeff, op_new_index)
    def generate_blocks(self,P):
        for i in range(1,self.L,1):
            self.blocks[i]=picos.HermitianVariable("block_nr{1}{0}".format(i, self.sec), (len(self.ops), len(self.ops)))
            P.add_constraint(self.blocks[i] >> 0)
            i1=1
            for o1 in self.ops:
                i2=1
                for o2 in self.ops:
        #             #print(i, self.blocks[i].shape, self.D_matrices["D_{0}{1}".format(i1,i2)].shape)
                    P.add_constraint(self.blocks[i][i1-1,i2-1]==self.D_matrices["D_{2}{0}/{1}".format(i1,i2, self.sec)][i,0])
                    i2+=1
                i1+=1 
class mom_symm_block_zero(mom_symm_block_base): 

    d_vectors=[] # the vector with the values d_1, d_2 ..., d_1=c_1L
    d_vector_map={}

    def __init__(self, ops, L, P, sec=""):    
        super().__init__(ops, L, P, sec=sec)
        # initialize the 
        # generate G matrices
        self.generate_G_matrices()
        for i in self.G_matrices.keys():
            print(i)
            for j in self.G_matrices[i].keys():
                print(j.sym,self.G_matrices[i][j].get() )
        #generate D matrices
        # self.generate_D_matrices(P)
        # # # # the Dvector are arranged (c_1sqrt(L), c_2sqrt(L)...)        
        # self.d_vectors=picos.ComplexVariable("d_vec", (len(ops),1))
        # ind=0
        # for o1 in ops:  
        #     coeff,nf=spin_class.normal_form(o1)
        #     self.d_vector_map[nf]=npa_funcs.nf_monomial(ind=ind,  prefac=coeff*np.sqrt(L)) # times  c_1 sqrt(L)
        #     ind+=1
        # # # associate all variables with each other   
        # ind=0
        # already_associated_elements=[]
        # for key in self.d_vector_map.keys():
        #     for G_key in self.G_matrices.keys():
        #         if key in self.G_matrices[G_key].keys():
        #             d_element=self.d_vector_map[key]
        #             G_element=self.G_matrices[G_key][key]
        #             P.add_constraint(self.d_vectors[d_element.ind,0]*d_element.prefac/np.sqrt(L)==G_element.prefac*self.G_variables[G_key][G_element.ind])
        # self.fix_ones_inGmatrices(P)
        # self.associate_Gelements(P,already_associated_elements) 
        # # # generate blocks
        # self.blocks[0]=picos.HermitianVariable("block_nr{1}{0}".format(0, self.sec), (len(ops)+1, len(ops)+1))
        # P.add_constraint(self.blocks[0] >> 0)
        # P.add_constraint(self.blocks[0][0,0]==1)
        # for i in range(1,len(ops)+1):
        # #for d_key in self.D_vector_map.keys():
        #     d_key=list(self.d_vector_map.keys())[i-1]
        #     d_element=self.d_vector_map[d_key]
        #     P.add_constraint(self.blocks[0][0,i]==d_element.prefac*self.d_vectors[d_element.ind,0])
        #     P.add_constraint(self.blocks[0][i,0]==np.conj(d_element.prefac)*self.d_vectors[d_element.ind,0])

        # i1=1
        # for o1 in self.ops:  
        #     i2=1
        #     for o2 in self.ops:               
        #         P.add_constraint(self.blocks[0][i1,i2]==self.D_matrices["D_{2}{0}/{1}".format(i1,i2, self.sec)][0,0])
        #         i2+=1
        #     i1+=1 
        
        # if len(self.ops)>1:
        #     self.generate_blocks(P)
        # self.associate(P)
        return
    
    
class mom_symm_block(mom_symm_block_base): 

    def __init__(self, ops, L, P, sec=""):    
        super().__init__(ops, L, P, sec=sec)
        # initialize the 
        # generate G matrices
        self.generate_G_matrices()
        #for i in self.G_matrices.keys():
        # generate D matrices
        self.generate_D_matrices(P)

        # # associate all variables with each other   
        ind=0
        already_associated_elements=[]
        self.fix_ones_inGmatrices(P)
        self.associate_Gelements(P,already_associated_elements) 
        # # generate blocks
        self.blocks[0]=picos.HermitianVariable("block_nr{1}{0}".format(0, self.sec), (len(ops), len(ops)))
        P.add_constraint(self.blocks[0] >> 0)
        P.add_constraint(self.blocks[0][0,0]==1)

        i1=0
        for o1 in self.ops:  
            i2=0
            for o2 in self.ops:               
                P.add_constraint(self.blocks[0][i1,i2]==self.D_matrices["D_{2}{0}/{1}".format(i1+1,i2+1, self.sec)][0,0])
                i2+=1
            i1+=1 
        if len(self.ops)>1:
            self.generate_blocks(P)
        return
            




















################### might not be needed#######################################
#     def sort(self, ops):
#           block_sizes={(0,0):1,}
#           i1=1
#           i2=1
#           # Generate G matrix elements
#           for o1 in ops:  
#               key="G_{0}{1}".format(i1,i2)
#               #print(len(self.G_matrices[key]))
#               block_sizes[(0,i1)]=len(self.G_matrices[key])
#               i1+=1
#           #mat_size=sum(block_sizes.values())
#           G = nx.Graph()
#           nr_node=0
#           for key in block_sizes.keys():
#               #print("key",key)
#               for j in range(block_sizes[key]):
#                   G.add_node(nr_node, matrix_ind=key, position_in_matrix=j)
#                   #print(G.nodes[nr_node])
                  
#                   #G.add_node(nr_node)
#                   nr_node+=1
#           G.add_edges_from([(0, 1)])
#           start=1
#           for i in range(1,len(block_sizes.keys())):
#               G.add_edges_from([(0, start)])
#               start+=block_sizes[list(block_sizes.keys())[i]]
#           #for n in G.nodes:
#            #   print(G.nodes[n])
#           #print(G.edges)
#           blocks=[]
#           bs_1=1 
#           for b1 in range(1,len(block_sizes.keys())):  
#               bs_2=1
#               for b2 in range(1,len(block_sizes.keys())):
#                   #print("added blocks", (bs_1,bs_2),block_sizes[b2])
#                   blocks+=[block((bs_1,bs_2),block_sizes[list(block_sizes.keys())[b2]])]
#                   bs_2+=block_sizes[list(block_sizes.keys())[b2]]
#               bs_1+=block_sizes[list(block_sizes.keys())[b1]] 
#           #print("blocks")
#           for bl in blocks:
#               #print(bl.diags)
#               for d in bl.diags:
#                   G.add_edge(d[0],d[1])
        
#           #print(G.nodes)
#           #print(G.edges)
#           S = [G.subgraph(c).copy() for c in nx.connected_components(G)]
#           # print(type(S[0]))
#           # print(S[0].nodes[0])
# #          for s in S:
#  #             print(s.nodes)
#   #            for l in s.nodes:
#    #               print(s.nodes[l])
#           return S
            