import numpy as np
from copy import deepcopy
import spin_class 
class spin_op_2d(spin_class.spin_op):
    """
    Class containing spin operators
    Convention 0,1,2,3
               4,5,6,7
               8,9,10,11
    """
    def __init__(self, cor, site, Lx):
        self.cor=cor # coordinate, x,y,z
        self.site=site # list with indices []
        self.Lx=Lx
        self.sym="S^"+self.cor+"_"+str(self.site)
    def pos(self):
        return self.site[0]*self.Lx+self.site[1]
    def get_translated(self,i, L):
        return spin_op_2d(self.cor,[self.site[0],(self.site[1]+i)%L], Lx=self.Lx)
    def get_xsite(self):
        return self.site[1]
