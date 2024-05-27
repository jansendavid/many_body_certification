import spin_class as spin_class
import picos
import numpy as np
import momentum_funcs as mfuncs
from copy import deepcopy
import npa_hams as npa_hams
import tensor as tensor
import itertools
from utils import sym, sx, sy,sz, I
def solve_marginal_problem(P,  m):
    """
    Implements marginal problem rho_2=tr_L[rho_3]=tr_R[rho_3] etc
    up to order m
    Note: for now, we use Picos hermitian variables to make compatible with 
    momen matrix. For better efficiency, they can be defined to be real and symmetric
    Parameters
    ----------
    P : Picos Problem
        DESCRIPTION.
    m : Level of hirarchy
        DESCRIPTION.

    Returns
    -------
    list_rhos : Dictioary 
        The reduced density matrix of i sites corresponds to entry .list_rhos[i]

    """
    list_rhos={}
    for i in range(1,m+1,1):
        rhonew_picos = picos.HermitianVariable("rho{0}".format(i), (2**i,2**i)) 
        list_rhos[i]=rhonew_picos
        P.add_constraint(list_rhos[i]>>0)

    P.add_constraint(picos.trace(list_rhos[2])==1)
    for i in range(2,m,1):
  
        P.add_constraint(picos.partial_trace(list_rhos[i+1], 0,2)==list_rhos[i])
        P.add_constraint(picos.partial_trace(list_rhos[i+1], i,2)==list_rhos[i])

    return list_rhos
def solve_marginal_problem_noTI(P,  m, name=""):
    """
    Solving marginal problem without translation symmetry.
    
    Implements marginal problem with only rho_2==tr_R[rho_3] etc
    up to order m
    Note: for now, we use Picos hermitian variables to make compatible with 
    momen matrix. For better efficiency, they can be defined to be real and symmetric
    Parameters
    ----------
    P : Picos Problem
        DESCRIPTION.
    m : Level of hirarchy
        DESCRIPTION.

    Returns
    -------
    list_rhos : Dictioary 
        The reduced density matrix of i sites corresponds to entry .list_rhos[i]

    """
    list_rhos={}
    for i in range(2,m+1,1):
        rhonew_picos = picos.HermitianVariable("rho{0}".format(i)+name, (2**i,2**i)) 
        list_rhos[i]=rhonew_picos
        P.add_constraint(list_rhos[i]>>0)

    P.add_constraint(picos.trace(list_rhos[2])==1)
    for i in range(2,m,1):
  
        #P.add_constraint(picos.partial_trace(list_rhos[i+1], 0,2)==list_rhos[i])
        P.add_constraint(picos.partial_trace(list_rhos[i+1], i,2)==list_rhos[i])

    return list_rhos

# TODO marginal problem with tree tensors

# def solve_marginal_problem_TTN(P,  m,isos, name=""):
#     """
#     Solve marginal problem using tree tensor CPTP map

#     Parameters
#     ----------
#     P : TYPE
#         DESCRIPTION.
#     m : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     list_rhos : TYPE
#         DESCRIPTION.

#     """
#     list_rhos={}
#     rhonew_picos = picos.HermitianVariable("rho{0}".format(4)+name, (2**4,2**4))
#     list_rhos[4]=rhonew_picos
#     P.add_constraint(picos.trace(list_rhos[4])==1)
#     P.add_constraint(picos.partial_trace(list_rhos[4], 0,2)==picos.partial_trace(list_rhos[4], 3,2))
#     P.add_constraint(list_rhos[4]>>0)
#     D=isos[0].shape[0]
#     d=isos[0].shape[1]
#     iso_kron=np.kron(isos[0].reshape(D,d*d),isos[0].reshape(D,d*d)).reshape(D**2,d**2,d**2)
#     rho_x=list_rhos[4].reshaped((2**8,1), order="C")
#     #print(iso_kron.shape)
#     rho_new=tensor.apply_CPTPmap(list_rhos[4], iso_kron)
#     rhonew_picos = picos.HermitianVariable("rho{0}".format(8)+name, ((D+1)**4,(D+1)**4)) 
#     list_rhos[8]=rhonew_picos
#     P.add_constraint(list_rhos[4]==picos.partial_trace(list_rhos[8], [0,1],(D+1)))
#     P.add_constraint(picos.partial_trace(list_rhos[8], [0],(D+1))==picos.partial_trace(list_rhos[8], [3],(D+1)))
#     #P.add_constraint(list_rhos[8]>>0)
    
#     # for i in range(2,m+1,1):
#     #     rhonew_picos = picos.HermitianVariable("rho{0}".format(i)+name, (2**i,2**i)) 
#     #     list_rhos[i]=rhonew_picos
#     #     P.add_constraint(list_rhos[i]>>0)

#     # P.add_constraint(picos.trace(list_rhos[2])==1)
#     # for i in range(2,m,1):
  
#     #     #P.add_constraint(picos.partial_trace(list_rhos[i+1], 0,2)==list_rhos[i])
#     #     P.add_constraint(picos.partial_trace(list_rhos[i+1], i,2)==list_rhos[i])

#     return list_rhos

# def get_rho_1(i, Ms,normal_forms):
#     # for now, only considering real rho's

#     rho=np.eye(2, dtype=complex)
#     matrix_dic={"z": sz, "x": sx,"y": sy}
    
#     for s in sym:
#         op1=spin_class.op_list([spin_class.spin_op(s,i)])
#         for symm_sector in normal_forms.keys():
#             if op1 in normal_forms[symm_sector].keys():
#                 new_obj=normal_forms[symm_sector][op1]
#                 rho=rho+np.conj(new_obj.prefac)*(Ms[symm_sector][new_obj.ind1,new_obj.ind2])*matrix_dic[s]
#                 break
#     # normalization
#     return rho/2


# def get_rho_2(i,j, Ms,normal_forms):
#     # for now, only considering real rho's

#     d=4
#     rho=np.eye(d, dtype=complex)
#     sxx=np.kron(sx, sx)
#     szz=np.kron(sz, sz)
#     syy=np.kron(sy, sy)
#     sx1=np.kron(sx, I)
#     sz1=np.kron(sz, I)
#     sy1=np.kron(sy,I)
#     s1x=np.kron( I,sx)
#     s1z=np.kron(I,sz)
#     s1y=np.kron(I,sy)
#     szx=np.kron(sz, sx)
#     sxz=np.kron(sx, sz)
#     sxy=np.kron(sx, sy)   
#     syx=np.kron(sy, sx)
#     szy=np.kron(sz, sy)
#     syz=np.kron(sy, sz)


    
#     matrix_dic={"zz": szz, "xx": sxx, "yy": syy,"zx": szx,"xz": sxz,\
#                 "zy": szy,"yz": syz, "yx": syx,"xy": sxy, "x1": sx1, "y1":sy1, "z1": sz1, "1z":s1z, "1x": s1x, "1y": s1y}
    
   
#     for s in sym:
#         for s_2 in sym:
#             for symm_sector in normal_forms.keys():
#                 op1=spin_class.op_list([spin_class.spin_op(s,i),spin_class.spin_op(s_2,j)])
#                 if op1 in normal_forms[symm_sector].keys():
#                     new_obj=normal_forms[symm_sector][op1]
#                     rho=rho+np.conj(new_obj.prefac)*(Ms[symm_sector][new_obj.ind1,new_obj.ind2])*matrix_dic[s+s_2]
#                     break
#     for s in sym:
#         for symm_sector in normal_forms.keys():
#             op1=spin_class.op_list([spin_class.spin_op(s,i)])
#             if op1 in normal_forms[symm_sector].keys():
#                 new_obj=normal_forms[symm_sector][op1]  
#                 #print(op1.sym, symm_sector, np.conj(new_obj.prefac))
#                 rho=rho+np.conj(new_obj.prefac)*(Ms[symm_sector][new_obj.ind1,new_obj.ind2])*matrix_dic[s+"1"]
#                 break
#     for s in sym:
#         for symm_sector in normal_forms.keys():        
#             op1=spin_class.op_list([spin_class.spin_op(s,j)])
#             if op1 in normal_forms[symm_sector].keys():
#                new_obj=normal_forms[symm_sector][op1]  
#                #print(op1.sym, symm_sector, np.conj(new_obj.prefac))
#                rho=rho+np.conj(new_obj.prefac)*(Ms[symm_sector][new_obj.ind1,new_obj.ind2])*matrix_dic["1"+s]
#                break

#     return rho/d                 
def general_rho(sites, Ms,normal_forms, class_name, **kwargs):
    """
    Generates a general reduced density matrix based on the moment matrix.
    

    Parameters
    ----------
    sites :  List
        Entries are sites [i,j,k] would generate rho_{ijk}.
    Ms : Picos Variable
        Moment matrix.
    normal_forms : Dict
        Keys are symmetrie sector, each entry is a dict where the normal form of the operator is stored.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    d=2**len(sites)
    rho=np.eye(d)
    sym_tot=["x","y", "z", "1"]
    total_comb=list(itertools.product(sym_tot, repeat=len(sites)))
    #for i in range(1,len(sites)):
    #    total_comb=list(itertools.product(total_comb, sym_tot))
    matrix_dic={"z": sz, "x": sx,"y": sy, "1": I}
    #print(total_comb)
    total_comb.pop(-1)
    #print(total_comb)
    for data in total_comb:
        Mat=matrix_dic[data[0]]
        #print(matrix_dic[data[0]])
        operator=[]
        if data[0]!="1":
            operator+=[spin_class.spin_op(data[0],sites[0])]
        for j in range(1,len(data)):
            Mat=np.kron(Mat, matrix_dic[data[j]])
            if data[j]!="1":
                operator+=[spin_class.spin_op(data[j],sites[j], **kwargs)]
        op_tot=spin_class.op_list(operator)
        #print(op_tot.sym)
        for symm_sector in normal_forms.keys(): 
            if op_tot in normal_forms[symm_sector].keys():
                #print("found",op_tot.sym )
                new_obj=normal_forms[symm_sector][op_tot]
                rho=rho+np.conj(new_obj.prefac)*(Ms[symm_sector][new_obj.ind1,new_obj.ind2])*Mat
                break

    return rho/d



def general_rho_TI(sites,basis):
    """
    Generates a general reduced density matrix based on the moment matrix.
    

    Parameters
    ----------
    sites :  List
        Entries are sites [i,j,k] would generate rho_{ijk}.
    Ms : Picos Variable
        Moment matrix.
    normal_forms : Dict
        Keys are symmetrie sector, each entry is a dict where the normal form of the operator is stored.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    d=2**len(sites)
    rho=np.eye(d)
    sym_tot=["x","y", "z", "1"]
    total_comb=list(itertools.product(sym_tot, repeat=len(sites)))
    #for i in range(1,len(sites)):
    #    total_comb=list(itertools.product(total_comb, sym_tot))
    matrix_dic={"z": sz, "x": sx,"y": sy, "1": I}
    #print(total_comb)
    total_comb.pop(-1)
    #print(total_comb)
    for data in total_comb:
        Mat=matrix_dic[data[0]]
        #print(matrix_dic[data[0]])
        operator=[]
        if data[0]!="1":
            operator+=[basis.class_name(data[0],sites[0],**basis.kwargs )]
        for j in range(1,len(data)):
            Mat=np.kron(Mat, matrix_dic[data[j]])
            if data[j]!="1":
                operator+=[basis.class_name(data[j],sites[j], **basis.kwargs)]
        op_tot=spin_class.op_list(operator)
        coeff,nf=spin_class.normal_form(op_tot,class_name=basis.class_name,**basis.kwargs)
        all_T=mfuncs.generate_all_translations(nf,basis.L,class_name=basis.class_name,**basis.kwargs)
        
        #print(op_tot.sym)
        for symm_sector in basis.sectors.keys(): 
            ops=list(set(all_T)& set(basis.sectors[symm_sector].map_TI.values()))
            if len(ops)>0:
                rho=rho+basis.sectors[symm_sector].map_TI_to_elements[ops[0]][0][0]*np.conj(basis.sectors[symm_sector].map_TI_to_elements[ops[0]][0][1])*Mat
                break

    return rho/d
