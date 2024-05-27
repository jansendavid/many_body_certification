
import spin_class as spin_class
import picos
import numpy as np
from copy import deepcopy
import npa_hams as npa_hams
import networkx as nx
from utils import sym

# convention, I store the operator and the coefficient in front
class nf_monomial:
    """
    Class containing the monomial position in the vector u and and their prefactor
    """
    def __init__(self, ind,   prefac):
        self.ind=ind
        self.prefac=prefac
    def get(self):
        return (self.ind,  self.prefac)


class matrix_element:
    """
    Class containing the matrix element position and their prefactor.
    So that if u[1]=\sigma^x_0 \sigma^y_1 , u[2]=\sigma^z_0->M[1,2]=\sigma^x_0 \sigma^y_1 \sigma^z_0, 
    NF(M[1,2])=-i \sigma^y_0 \sigma^y_1=prefac*\sigma^y_0 \sigma^y_1
    (-i) is store as prefactor: so \sigma^y_0 \sigma^y_0=conj(prefac)M[1,2]
    """
    def __init__(self, ind1, ind2,   prefac):
        self.ind1=ind1
        self.ind2=ind2
        self.prefac=prefac
    def get(self):
        return (self.ind1,self.ind2, self.prefac)
        
class basis():
    def __init__(self,ops, class_name, **kwargs):
        self.ops=ops
        
        self.class_name=class_name
        self.kwargs=kwargs
        self.total_vectors={}
        
        for symm in ops.keys():
            vecs=self.get_normal_form_vector( ops=self.ops[symm],class_name=self.class_name,kwargs=self.kwargs)
            self.total_vectors[symm]=vecs
        self.total_nfs=self.get_NF_matrix_elements()
            
        
    
    def get_mom_matrix_M(self,P):
        """
        Makes moment matrix M for different symmetrie sectors
    
        Parameters
        ----------
        P : Picos Problem
            DESCRIPTION.
        vector : Dict, keys: symmetrie sector, value Dics: key=monomial, value=nf_monomial class
            Dictionary containing all symmetri sectors which dictionary containing all monomials in the block.
        normal_forms : Dict, keys: symmetrie sector, value Dics: key=monomial, value=matrix_element
            A dictionary containing all symmetrie sectors, each with a dictionary containing all the normal
            forms in the block matrix.
    
        Returns
        -------
        total_Ms : Dict
            Dictiobary containing symmetri sectors as keys and momentmatrix as values.
    
        """
    
        total_Ms={}
        for symm_sector in self.total_vectors.keys():
            vector_symm=self.total_vectors[symm_sector]
            sector_size=len(vector_symm)
            M = picos.HermitianVariable("M_{0}".format(symm_sector), (sector_size,sector_size)) 
            P.add_constraint(M >> 0)
            for i in vector_symm.keys():
                for j in vector_symm.keys():
                    inds_i=vector_symm[i].ind
                    inds_j=vector_symm[j].ind
                    # going through operators in reverse since we  have u^{dagger}u
                    op_j=spin_class.op_list(j.ops[::-1])
                    element=op_j.ops+i.ops
                    total_operator=spin_class.op_list(element)
                    # finding normal form NF(u) and the coefficient so that 
                    #  NF(\sigma^x_0 \sigma^y_1 \sigma^z_0)=coeff* \sigma^y_0 \sigma^y_1
                    coeff_tot,nf_tot=spin_class.normal_form(total_operator, self.class_name, **self.kwargs)
                    if i==j:
                        # the diagonal term is always 1 for spin systems
                        P.add_constraint(M[inds_i,inds_i]==1)
                    else:
                        # finding matrix element of the normal form
                        new_obj=self.total_nfs[symm_sector][nf_tot]
                        
                        if (inds_j,inds_i)!=(new_obj.ind1,new_obj.ind2):
                            #print("Now:",symm_sec, np.conj(coeff_tot), total_el.sym , "->", nf_tot.sym )
                            # relating the two matrix elements
                            P.add_constraint(coeff_tot*M[inds_j,inds_i]==new_obj.prefac*M[new_obj.ind1,new_obj.ind2])
                        
            total_Ms[symm_sector]=M
            # fix ones
            for nf_sym in self.total_nfs.keys():
                for nf_key in self.total_nfs[nf_sym].keys():
                    if nf_key.sym=="1":
                        obj=self.total_nfs[nf_sym][nf_key]
                        P.add_constraint((obj.prefac)*M[obj.ind1,obj.ind2]==1)
                    
            
            
        # relate elements different moment sectors, probably a bottle neck which is not needed often
        # if a alament is already associated, we do not add it again
        already_associated_elements=[]
        for i in range(len(self.total_nfs.keys())-1):
            first_sector=list(self.total_nfs.keys())[i]
            first_set=self.total_nfs[first_sector]
            # j will go through all remaining moment sectors
            for j in range(i+1,len(self.total_nfs.keys())):
                second_sector=list(self.total_nfs.keys())[j]
                sec_set=self.total_nfs[second_sector]
                overlap=list(set(first_set.keys()) & set(sec_set.keys()))
                #print("secor")
                #print(first_sector, second_sector)
                for o in overlap:
                    #print( o.sym)
                    # If two sectors contain the same elemet, we relate them to each other
                    #print("relate")
                   # print("overlap", first_sector, second_sector,o.sym, (first_set[o].ind1,first_set[o].ind2),(sec_set[o].ind1,sec_set[o].ind2))
                    if (j,o) not in already_associated_elements:
                    #    print("relate", first_set[o].prefac,"tp", sec_set[o].prefac,(first_set[o].ind1,first_set[o].ind2),(sec_set[o].ind1,sec_set[o].ind2))
                        P.add_constraint(first_set[o].prefac*total_Ms[first_sector][first_set[o].ind1,first_set[o].ind2 ]==sec_set[o].prefac*total_Ms[second_sector][sec_set[o].ind1,sec_set[o].ind2 ])
                        already_associated_elements+=[(j,o)]
        return total_Ms
    
    def get_NF_matrix_elements(self):
        """
        Generates a dictionary with all the symmetrie sectors and the normal form matrix elements
    
        Parameters
        ----------
        vector : TYPE
            DESCRIPTION.
    
        Returns
        -------
        total_nfs : TYPE
            DESCRIPTION.
    
        """
        #finds the normal form an sign of all the monomials in states
        total_nfs={}
        unit_sec=True
        for symmsec in self.total_vectors.keys():
            symm_vector=self.total_vectors[symmsec]
            dict_normal_forms={}
            for i in symm_vector.keys():
                
                if unit_sec:
                    dict_normal_forms[i]=matrix_element(ind1=0,ind2=symm_vector[i].ind,  prefac=symm_vector[i].prefac)
                for j in self.total_vectors[symmsec].keys():
                    op_j=spin_class.op_list(j.ops[::-1])
                    #print(j.sym, op_j.sym)
                    element=op_j.ops+i.ops
                    total_el=spin_class.op_list(element)
                    coeff_tot,nf_tot=spin_class.normal_form(total_el, self.class_name, **self.kwargs)
                    if nf_tot.ops!=[] and nf_tot not in dict_normal_forms.keys():
                        dict_normal_forms[nf_tot]=matrix_element(ind1=symm_vector[j].ind,ind2=symm_vector[i].ind,  prefac=coeff_tot)
            unit_sec=False
            total_nfs[symmsec]=dict_normal_forms
        return total_nfs
    
    
    
                
    def get_normal_form_vector(self,ops,  class_name, **kwargs):
        """
        Generates a vector with the normal forms of the monomials.
    
        Parameters
        ----------
        ops : TYPE
            DESCRIPTION.
        sym_sec : TYPE
            DESCRIPTION.
        index : TYPE, optional
            DESCRIPTION. The default is 0.
    
        Returns
        -------
        dict_normal_forms : TYPE
            DESCRIPTION.
    
        """
        index=0
        dict_normal_forms={}
        for op in ops:
            coeff,nf=spin_class.normal_form(op,class_name, **kwargs)
            if nf not in dict_normal_forms.keys():
                dict_normal_forms[nf]=nf_monomial(ind=index,  prefac=coeff)
                index+=1
        return dict_normal_forms
    
    
    
    # # def apply_virial_theorem(P,M,nf_vector, total_nfs, J, h, shape,L, PB):
    # #     VT=npa_hams.get_tfi_npa_virial_sigy(states=nf_vector,dict_normal_forms=total_nfs, J=J, h=h, shape=M.shape,L=L, PB=PB)
    # #     VT_2=npa_hams.get_tfi_npa_virial_sigz(states=nf_vector,dict_normal_forms=total_nfs, J=J, h=h, shape=M.shape,L=L, PB=PB)
    
    # #     VY_picos = picos.Constant("VTy", VT)
    # #     VY_2_picos = picos.Constant("VTz", VT)
    # #     P.add_constraint((M | VY_picos )==0)
    # #     P.add_constraint((M | VY_2_picos )==0)
    # #     return
    
    
                    
