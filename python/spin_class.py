import numpy as np
from copy import deepcopy

class spin_op:
    """
    Class containing spin operators
    """
    def __init__(self, cor, site):
        self.cor=cor # coordinate, x,y,z
        self.site=site
        self.sym="S^"+self.cor+"_"+str(self.site)
    def pos(self):
        return self.site # default
        
    def __lt__(self, other):
        # on sites
        return self.site < other.site
    def __le__(self, other):
        # on sites
        return self.site <= other.site
    def __eq__(self, other):
        # on sites
        return self.site == other.site
    def __gt__(self, other):
        # on sites
        return self.site > other.site
    def __ge__(self, other):
        # on sites
        return self.site >= other.site
    def get_translated(self,i, L):
        return spin_op(self.cor,(self.site+i)%L)
    def get_xsite(self):
        return self.site
        

def get_sym(a,b):
    if a=="x" and b=="y":
        return (1j, "z")
    elif a=="y" and b=="x":
        return (-1j, "z")
    elif a=="z" and b=="x":
        return (1j, "y")
    elif a=="x" and b=="z":
        return (-1j, "y")
    elif a=="y" and b=="z":
        return (1j, "x")
    elif a=="z" and b=="y":
        return (-1j, "x")
    else:
        return (0,0)
class op_list:
    """
    List of operators, e.g., S^{z}_0 S^{x}_1 
    """
    def __init__(self, ops):
        self.ops=ops
        self.sym=""
        if len(ops)==0:
            self.sym="1"
        else:
            for o in self.ops:
                self.sym+=o.sym
    def __hash__(self):
        return hash(self.sym)
    def __eq__(self, other):
        # on sites    
        return self.sym == other.sym

def commute(op_1, op_2):
    if op_1.site>op_2.site:
        return (1,[op_2, op_1])
    elif op_1.site==op_2.site:
        coeff, op=get_sym(op_1.cor,op_2.cor)
        #print(coeff, op)
        return (coeff, [spin_op(op,op_1.site)])
    else:
        return (1, [op_1, op_2])
def normal_form(ops, class_name, **kwargs):
    """
    Generates normal for of operator ops.
    The convetion is S^x_0S^y_0 -> 1j,S^z_0 


    Parameters
    ----------
    ops : op_list
        DESCRIPTION.
    class_name: name of class to be called e.g., (spin 1d, spin 2d etc)
    kwargs: extra parameters for the spin class
    Returns
    -------
    TYPE
        DPreafactor.
    TYPE
        op_list of normal form.

    """
    list_cp=deepcopy(ops.ops)
    list_cp=op_list(sorted(list_cp, key=lambda x: x.pos()))

    def run_loop(list_cp,class_name, **kwargs):
        pref=1
        list_cp_final={}
        for i in range((len(list_cp.ops))):
        #print("sym",list_cp.ops[i].sym)
            list_cp_final[i]=list_cp.ops[i]
        for i in range(len(list_cp.ops)-1):
            if i in list_cp_final.keys():
                o1=list_cp_final[i]
                o2=list_cp_final[i+1]
                if o1.site!=o2.site:
                    continue
                elif o1.sym==o2.sym:
                    del list_cp_final[i]
                    del list_cp_final[i+1]
                else:
                    coeff, op=get_sym(o1.cor,o2.cor)
                    del list_cp_final[i]
                    del list_cp_final[i+1]
                    list_cp_final[i]=class_name(op,o1.site, **kwargs)
                    pref*=coeff

        final_op=[]
        for o in sorted(list_cp_final.keys()):
            final_op+=[list_cp_final[o]]
    
        return (pref, op_list(final_op))
    final_op_new=op_list([])
    new_list=list_cp
    pref=1
    #dummy variable that is not new_list
    new_op=op_list([]) 
    # iterate until normal form is found
    # check if state changes
    change=False
    while (not change):
        new_pref, new_op=run_loop(new_list,class_name, **kwargs)
        pref*=new_pref
        if new_op==new_list:
            # state is in normal form
            change=True
        else:
            new_list=new_op
        #    break
        
    new_list=new_op
    return (pref, new_list)



