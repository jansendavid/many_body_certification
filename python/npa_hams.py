
import spin_class as spin_class
import picos
import numpy as np
from utils import sym

def get_xxz_npa(normal_forms, J, Delta,shapes,L, PB=False):
    """
    XXZ Hamiltonian H=\sum_i J(S^{x}_i_S^{x}_{i+1}+S^{y}_i_S^{y}_{i+1})+ Delta S^{z}_i_S^{z}_{i+1}

    Parameters
    ----------
    normal_forms : TYPE
        DESCRIPTION.
    J : TYPE
        DESCRIPTION.
    Delta : TYPE
        DESCRIPTION.
    shapes : TYPE
        DESCRIPTION.
    L : TYPE
        DESCRIPTION.
    PB : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    Is : TYPE
        DESCRIPTION.

    """
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k], dtype=complex)
    end=L-1
    if PB:
        end=L
    for i in np.arange(0,end,1):
        
        for s in sym:          
            o_1=spin_class.spin_op(s,min(i, (i+1)%L))
            o_2=spin_class.spin_op(s,max(i, (i+1)%L))
            otot=spin_class.op_list([o_1,o_2])
            add_term=True
            for symmetrie in normal_forms.keys():
                if otot in normal_forms[symmetrie].keys() and add_term:
                    matel=normal_forms[symmetrie][otot]
                    if s=="z":
                        Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*Delta/4
                        Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*Delta/4
                        
                    else:
                        Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*J/4
                        Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*J/4
                    break
    for k in Is.keys():
        #print("np", k, Is[k])
        Is[k]=0.5*Is[k]
    return Is
def get_tfi_npa(normal_forms, J, h,shapes,L, PB=False):
    """
    Transverse field ising model Hamiltonian H=\sum_i J(S^{x}_i_S^{x}_{i+1}+S^{y}_i_S^{y}_{i+1})+ Delta S^{z}_i_S^{z}_{i+1}

    """
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k])
    end=L-1
    if PB:
        end=L
    for i in np.arange(0,L,1):
        s="x"
        o=spin_class.spin_op(s,i)
        for symm in normal_forms.keys():
            if spin_class.op_list([o]) in normal_forms[symm].keys():
                matel=normal_forms[symm][spin_class.op_list([o])]
                #print(o.sym)
                #print(ind)
                Is[symm][matel.ind1,matel.ind2]+=1.*h/4
                Is[symm][matel.ind2,matel.ind1]+=1.*h/4
                break
    for i in np.arange(0,end,1):
        s="z"
        o_1=spin_class.spin_op(s,min(i, (i+1)%L))
        o_2=spin_class.spin_op(s,max(i, (i+1)%L))
        otot=spin_class.op_list([o_1,o_2])
        for symm in normal_forms.keys():
            if otot in normal_forms[symm].keys():
                matel=normal_forms[symm][spin_class.op_list([o_1,o_2])]
                Is[symm][matel.ind1,matel.ind2]+=1.*J/4
                Is[symm][matel.ind2,matel.ind1]+=1.*J/4
                break
    for k in Is.keys():
        Is[k]=0.5*Is[k]
    return Is
# MBL Hamiltonian
def get_mbl_npa(normal_forms, J, hx,hz,shapes,L, PB=False):
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k])
    end=L-1
    if PB:
        end=L
    for i in np.arange(0,L,1):
        s="x"
        o=spin_class.op_list([spin_class.spin_op(s,i)])
        for symm in normal_forms.keys(): 
            if o in normal_forms[symm].keys():
                matel=normal_forms[symm][spin_class.op_list([o])]
                Is[symm][matel.ind1,matel.ind2]+=1.*hx[i]/2
                Is[symm][matel.ind2,matel.ind1]+=1.*hx[i]/2
            break
        # s="z"
        # o=spin_class.spin_op(s,i)
        # for symm in normal_forms.keys():     
        #     matel=normal_forms[symm][spin_class.op_list([o])]
        #     Is[symm][matel.ind1,matel.ind2]+=1.*hz[i]/2
        #     Is[symm][matel.ind2,matel.ind1]+=1.*hz[i]/2
        #     break

    for i in np.arange(0,end,1):
        
        for s in sym:          
            o_1=spin_class.spin_op(s,min(i, (i+1)%L))
            o_2=spin_class.spin_op(s,max(i, (i+1)%L))
            otot=spin_class.op_list([o_1,o_2])
            add_term=True
            for symmetrie in normal_forms.keys():
                if otot in normal_forms[symmetrie].keys() and add_term:
                    matel=normal_forms[symmetrie][otot]
                    Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*J/4
                    Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*J/4
                    break
                    # if s=="z":
                    #     Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*Delta/4
                    #     Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*Delta/4
                        
                    # else:
                    #     Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*J/4
                    #     Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*J/4
                    # break

    for k in Is.keys():
        Is[k]=0.5*Is[k]
    return Is
def get_mbl_commleft_npa(normal_forms, J, hx,hz,shapes,L, OP,PB=False):
    # computed [H,x]=Hx
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k])
    end=L-1
    if PB:
        end=L
    for i in np.arange(0,L,1):
        s="x"
        op_new=spin_class.op_list([spin_class.spin_op(s,i)]+OP.ops)
        for symm in normal_forms.keys():   
            coeff,nf=spin_class.normal_form(op_new)
            if nf in normal_forms[symm].keys():
                #print("hehe")
                matel=normal_forms[symm][nf]
                Is[symm][matel.ind1,matel.ind2]+=1.*hx[i]*coeff/2
                Is[symm][matel.ind2,matel.ind1]+=1.*hx[i]*np.conj(coeff)/2
            break
        # s="z"
        # o=spin_class.spin_op(s,i)
        # for symm in normal_forms.keys():     
        #     matel=normal_forms[symm][spin_class.op_list([o])]
        #     Is[symm][matel.ind1,matel.ind2]+=1.*hz[i]/2
        #     Is[symm][matel.ind2,matel.ind1]+=1.*hz[i]/2
        #     break

    for i in np.arange(0,end,1):
        
        for s in sym:          
            o_1=spin_class.spin_op(s,min(i, (i+1)%L))
            o_2=spin_class.spin_op(s,max(i, (i+1)%L))
            otot=spin_class.op_list([o_1,o_2]+OP.ops)
            coeff,nf=spin_class.normal_form(otot)
            
            add_term=True
            for symmetrie in normal_forms.keys():
                if nf in normal_forms[symmetrie].keys() and add_term:
                    #print("haha")
                    matel=normal_forms[symmetrie][nf]
                    Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*J/4
                    Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*J/4
                    break
                    # if s=="z":
                    #     Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*Delta/4
                    #     Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*Delta/4
                        
                    # else:
                    #     Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*J/4
                    #     Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*J/4
                    # break

    for k in Is.keys():
        Is[k]=0.5*Is[k]
    return Is
def get_mbl_commright_npa(normal_forms, J, hx,hz,shapes,L, OP,PB=False):
    # computed [H,x]=Hx
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k])
    end=L-1
    if PB:
        end=L
    for i in np.arange(0,L,1):
        s="x"
        op_new=spin_class.op_list(OP.ops+[spin_class.spin_op(s,i)])
        for symm in normal_forms.keys():   
            coeff,nf=spin_class.normal_form(op_new)
            if nf in normal_forms[symm].keys():
                #print("hehe")
                matel=normal_forms[symm][nf]
                Is[symm][matel.ind1,matel.ind2]+=1.*hx[i]*coeff/2
                Is[symm][matel.ind2,matel.ind1]+=1.*hx[i]*np.conj(coeff)/2
            break
        # s="z"
        # o=spin_class.spin_op(s,i)
        # for symm in normal_forms.keys():     
        #     matel=normal_forms[symm][spin_class.op_list([o])]
        #     Is[symm][matel.ind1,matel.ind2]+=1.*hz[i]/2
        #     Is[symm][matel.ind2,matel.ind1]+=1.*hz[i]/2
        #     break

    for i in np.arange(0,end,1):
        
        for s in sym:          
            o_1=spin_class.spin_op(s,min(i, (i+1)%L))
            o_2=spin_class.spin_op(s,max(i, (i+1)%L))
            otot=spin_class.op_list(OP.ops+[o_1,o_2])
            coeff,nf=spin_class.normal_form(otot)
            
            add_term=True
            for symmetrie in normal_forms.keys():
                if nf in normal_forms[symmetrie].keys() and add_term:
                    #print("haha")
                    matel=normal_forms[symmetrie][nf]
                    Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*J/4
                    Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*J/4
                    break
                    # if s=="z":
                    #     Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*Delta/4
                    #     Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*Delta/4
                        
                    # else:
                    #     Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*J/4
                    #     Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*J/4
                    # break

    for k in Is.keys():
        Is[k]=0.5*Is[k]
    return Is
def get_mbl_squared_npa(normal_forms, J, hx,hz,shapes,L, PB=False):
    # Squared MBL Hamiltonian
    Is={}
    for k in normal_forms.keys():
        Is[k]=np.zeros(shapes[k], dtype=complex)
    end=L-1
    if PB:
        end=L
    const=0
    # sx_i sx_j term
    for i in np.arange(0,L,1):
        for j in np.arange(0,L,1):
            if i!=j:
                s="x"
                
                op_tot=spin_class.op_list([spin_class.spin_op(s,min(i,j)),spin_class.spin_op(s,max(i,j))])
                for symm in normal_forms.keys():     
                    matel=normal_forms[symm][op_tot]
                    Is[symm][matel.ind1,matel.ind2]+=1.*hx[i]*hx[j]/4
                    Is[symm][matel.ind2,matel.ind1]+=1.*hx[i]*hx[j]/4
                    break
            else:
                Is[list(normal_forms.keys())[0]][0,0]+=2*1.*hx[i]*hx[i]/4
                
        # s="z"
        # o=spin_class.spin_op(s,i)
        # for symm in normal_forms.keys():     
        #     matel=normal_forms[symm][spin_class.op_list([o])]
        #     Is[symm][matel.ind1,matel.ind2]+=1.*hz[i]/2
        #     Is[symm][matel.ind2,matel.ind1]+=1.*hz[i]/2
        #     break

    for i in np.arange(0,end,1):
        for j in np.arange(0,end,1):
        
            for s in sym: 
                for s2 in sym: 
                    o_1=spin_class.spin_op(s,min(i, (i+1)%L))
                    o_2=spin_class.spin_op(s,max(i, (i+1)%L))
                    o_3=spin_class.spin_op(s2,min(j, (j+1)%L))
                    o_4=spin_class.spin_op(s2,max(j, (j+1)%L))
                    total_operator=spin_class.op_list([o_1,o_2, o_3, o_4])
                    coeff_tot,nf_tot=spin_class.normal_form(total_operator)     
                    for symmetrie in normal_forms.keys():
                        if nf_tot in normal_forms[symmetrie].keys():
                            matel=normal_forms[symmetrie][nf_tot]
                            Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(coeff_tot)*np.conj(matel.prefac)*1.*J*J/8
                            Is[symmetrie][matel.ind2,matel.ind1]+=coeff_tot*(matel.prefac)*1.*J*J/8
                            break
                    if len(nf_tot.ops)==0:
                        Is[list(normal_forms.keys())[0]][0,0]+=2*1.*J*J/8
                        
    for i in np.arange(0,end,1):
        for j in np.arange(0,end,1):
            for s in sym: 
                o_1=spin_class.spin_op(s,min(i, (i+1)%L))
                o_2=spin_class.spin_op(s,max(i, (i+1)%L))
                o_3=spin_class.spin_op("x",j)
    
                total_operator=spin_class.op_list([o_1,o_2, o_3])
                coeff_tot,nf_tot=spin_class.normal_form(total_operator)     
                for symmetrie in normal_forms.keys():
                    if nf_tot in normal_forms[symmetrie].keys():
                        matel=normal_forms[symmetrie][nf_tot]
                        Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(coeff_tot)*np.conj(matel.prefac)*1.*J*J/8
                        Is[symmetrie][matel.ind2,matel.ind1]+=coeff_tot*(matel.prefac)*1.*J*J/8
                        break                       
    for i in np.arange(0,end,1):
        for j in np.arange(0,end,1):
            for s in sym: 
                o_1=spin_class.spin_op("x",j)
                o_2=spin_class.spin_op(s,min(i, (i+1)%L))
                o_3=spin_class.spin_op(s,max(i, (i+1)%L))
                
    
                total_operator=spin_class.op_list([o_1,o_2, o_3])
                coeff_tot,nf_tot=spin_class.normal_form(total_operator)     
                for symmetrie in normal_forms.keys():
                    if nf_tot in normal_forms[symmetrie].keys():
                        matel=normal_forms[symmetrie][nf_tot]
                        Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(coeff_tot)*np.conj(matel.prefac)*1.*J*J/8
                        Is[symmetrie][matel.ind2,matel.ind1]+=coeff_tot*(matel.prefac)*1.*J*J/8
                        break  
    for k in Is.keys():
        Is[k]=0.5*Is[k]
    return Is



# def get_mbl_npa_local(normal_forms,i,j, J,Delta, hx,hz,shapes,L, PB=False):
#     Is={}
#     for k in normal_forms.keys():
#         Is[k]=np.zeros(shapes[k])
#     end=L-1
#     if PB:
#         end=L
#     s="x"
#     o=spin_class.spin_op(s,i)
#     for symm in normal_forms.keys():     
#         matel=normal_forms[symm][spin_class.op_list([o])]
#         Is[symm][matel.ind1,matel.ind2]+=1.*hx[i]/2
#         Is[symm][matel.ind2,matel.ind1]+=1.*hx[i]/2
#         break
#     o=spin_class.spin_op(s,j)
#     for symm in normal_forms.keys():     
#         matel=normal_forms[symm][spin_class.op_list([o])]
#         Is[symm][matel.ind1,matel.ind2]+=1.*hx[i]/2
#         Is[symm][matel.ind2,matel.ind1]+=1.*hx[i]/2
#         break
#         # s="z"
#         # o=spin_class.spin_op(s,i)
#         # for symm in normal_forms.keys():     
#         #     matel=normal_forms[symm][spin_class.op_list([o])]
#         #     Is[symm][matel.ind1,matel.ind2]+=1.*hz[i]/2
#         #     Is[symm][matel.ind2,matel.ind1]+=1.*hz[i]/2
#         #     break
#     for s in sym:          
#         o_1=spin_class.spin_op(s,min(i, (j)))
#         o_2=spin_class.spin_op(s,max(i, (j)))
#         otot=spin_class.op_list([o_1,o_2])
#         add_term=True
#         for symmetrie in normal_forms.keys():
#             if otot in normal_forms[symmetrie].keys() and add_term:
#                 matel=normal_forms[symmetrie][otot]
#                 if s=="z":
#                     Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*Delta/4
#                     Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*Delta/4
#                 else:
#                     Is[symmetrie][matel.ind1,matel.ind2]+=np.conj(matel.prefac)*1.*J/4
#                     Is[symmetrie][matel.ind2,matel.ind1]+=(matel.prefac)*1.*J/4
#                 break

#     for k in Is.keys():
#         Is[k]=0.5*Is[k]
#     return Is

# TODO!
# def get_tfi_squared_npa(states,dict_normal_forms, J, h,shape,L, PB=False):
#     I=np.zeros((shape))
#     for i in np.arange(0,L,1):
#             for j in np.arange(0,L,1):
#                 s="x"
#                 o1=spin_class.spin_op(s,i)
#                 o2=spin_class.spin_op(s,i)
        
#                 ind=states[spin_class.op_list([o])]
#                 #print(o.sym)
#                 #print(ind)
#                 I[ind[0],ind[1]]+=1.*h/4
#                 I[ind[1],ind[0]]+=1.*h/4
#     for i in np.arange(0,L-1,1):
#         s="z"
#         o_1=spin_class.spin_op(s,i)
#         o_2=spin_class.spin_op(s,(i+1)%L)
#         coeff,nf=spin_class.normal_form(spin_class.op_list([o_2,o_1]))
#         ind=dict_normal_forms[spin_class.op_list([o_1,o_2])]
#         I[ind[0],ind[1]]+=1.*J/4
#         I[ind[1],ind[0]]+=1.*J/4
#     if PB:
#         s="z"
#         o_1=spin_class.spin_op(s,0)
#         o_2=spin_class.spin_op(s,L-1)
#         #print(op_list([o_1,o_2]).sym)
#         coeff,nf=spin_class.normal_form(spin_class.op_list([o_2,o_1]))
#         ind=dict_normal_forms[spin_class.op_list([o_1,o_2])]    
#         I[ind[0],ind[1]]+=1.*J/4
#         I[ind[1],ind[0]]+=1.*J/4
#     I=0.5*I
#     return I
