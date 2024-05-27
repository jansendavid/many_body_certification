from ncon import ncon
import numpy as np
from utils import sx, sy,sz,I
#fGoing from paulis to S operators
sX = sx/2
sY = sy/2
sZ = sz/2
sI = I
# cup=np.array([[0,0,1,0],
#               [0,0,0,1],
#               [0,0,0,0],
#               [0,0,0,0]])

# cdo=np.array([[0,1,0,0],
#               [0,0,0,0],
#               [0,0,0,-1],
#               [0,0,0,0]])

# F=np.diag([1,-1,-1,1])
# nup=np.diag([0,0,1,1])
# ndo=np.diag([0,1,0,1])

# nup2=nup-0.5*np.eye(4);
# ndo2=ndo-0.5*np.eye(4);
def TFI(J, h):
    """
    # Note that we transform the first site by rotating around z-axis by pi 
    for better convergence 
    Construct the spin-1 Heisenberg Hamiltonian for given couplings.
    
        Parameters
        ----------
        Jx : float
            Coupling strength in x direction
        Jy : float
            Coupling strength in y direction
        Jy : float
            Coupling strength in z direction
        hz : float
            Coupling for Sz terms

        Returns
        -------
        h : np.array (3, 3, 3, 3)
            Spin-1 Heisenberg Hamiltonian.
    """

    I = np.eye(2)

    return J*ncon((sZ, sZ), ([-1, -3], [-2, -4])) \
            + 0.25*h*ncon((sI, sX), ([-1, -3], [-2, -4])) + 0.25*h*ncon((sX, sI), ([-1, -3], [-2, -4]))
# def Mag():

#     Sz = 0.5*np.array([[1, 0], [ 0, -1]])
#     I = np.eye(2)

#     return 0.5*ncon((sI,sZ), ([-1, -3], [-2, -4])) + 0.5*ncon((sZ, sI), ([-1, -3], [-2, -4]))
             

# def Mag_XXZ_d4():

#     d=4
#     I = np.eye(2)
#     m_temp=0.5*ncon((sI,sZ), ([-1, -3], [-2, -4])) + 0.5*ncon((sZ, sI), ([-1, -3], [-2, -4]))
#     return (0.5 * np.kron(m_temp, np.eye(4)) + 0.5 * np.kron(np.eye(4), m_temp) +
#         np.kron(np.eye(2), np.kron(m_temp, np.eye(2)))).reshape(d, d, d, d)
             



# def XXZ_chain_d4(J, Delta,h):
#     """
#     Construct the spin-1 Heisenberg Hamiltonian for given couplings.
    
#         Parameters
#         ----------
#         Jx : float
#             Coupling strength in x direction
#         Jy : float
#             Coupling strength in y direction
#         Jy : float
#             Coupling strength in z direction
#         hz : float
#             Coupling for Sz terms

#         Returns
#         -------
#         h : np.array (3, 3, 3, 3)
#             Spin-1 Heisenberg Hamiltonian.
#     """
#     d = 4
#     # diviging by factor 2/4 going from pauli to spin
#     h_temp = (np.real(J*np.kron(sX, sX) + J*np.kron(sY, sY) + Delta*np.kron(sZ, sZ)+0.5*h*np.kron(sI, sZ)+0.5*h*np.kron(sZ, sI)))
#     h0 = (0.5 * np.kron(h_temp, np.eye(4)) + 0.5 * np.kron(np.eye(4), h_temp) +
#         np.kron(np.eye(2), np.kron(h_temp, np.eye(2)))).reshape(d, d, d, d)
#     # making local dimension two sites=4 for better convergenceSB7WJy
#     return h0

def XXZ_chain_spinflip(J, Delta,h):
    """
    WITH SPIN FLIP!
    Construct the spin-1 Heisenberg Hamiltonian for given couplings.
    Used for better convergence of VUMPS, see original paper.
    
        Parameters
        ----------
        Jx : float
            Coupling strength in x direction
        Jy : float
            Coupling strength in y direction
        Jy : float
            Coupling strength in z direction
        hz : float
            Coupling for Sz terms

        Returns
        -------
        h : np.array (3, 3, 3, 3)
            Spin-1 Heisenberg Hamiltonian.
    """
    d = 2
    # diviging by factor 2/4 going from pauli to spin
    h0 = (np.real(J*np.kron(-sX, sX) + J*np.kron(-sY, sY) + Delta*np.kron(sZ, sZ)+0.5*h*np.kron(sI, sZ)+0.5*h*np.kron(sZ, sI)))

    return h0.reshape(d, d, d, d).reshape(d, d, d, d)

def XXZ_chain(J, Delta,h):
    """
    Construct the spin-1 Heisenberg Hamiltonian for given couplings.
    
        Parameters
        ----------
        Jx : float
            Coupling strength in x direction
        Jy : float
            Coupling strength in y direction
        Jy : float
            Coupling strength in z direction
        hz : float
            Coupling for Sz terms

        Returns
        -------
        h : np.array (3, 3, 3, 3)
            Spin-1 Heisenberg Hamiltonian.
    """
    d = 2
    # diviging by factor 2/4 going from pauli to spin
    h0 = (np.real(J*np.kron(sX, sX) + J*np.kron(sY, sY) + Delta*np.kron(sZ, sZ)+0.5*h*np.kron(sI, sZ)+0.5*h*np.kron(sZ, sI)))

    return h0.reshape(d, d, d, d).reshape(d, d, d, d)
def MBL_chain(J, Delta,h):
    """
    WITH SPIN FLIP!
    Construct the spin-1 Heisenberg Hamiltonian for given couplings.
    
        Parameters
        ----------
        Jx : float
            Coupling strength in x direction
        Jy : float
            Coupling strength in y direction
        Jy : float
            Coupling strength in z direction
        hz : float
            Coupling for Sz terms

        Returns
        -------
        h : np.array (3, 3, 3, 3)
            Spin-1 Heisenberg Hamiltonian.
    """
    d = 2
    # diviging by factor 2/4 going from pauli to spin
    h0 = (np.real(J*np.kron(-sX, sX) + J*np.kron(-sY, sY) + Delta*np.kron(sZ, sZ)+0.5*h[0]*np.kron(sI, sZ)+0.5*h[1]*np.kron(sZ, sI)))

    return h0.reshape(d, d, d, d).reshape(d, d, d, d)


# def Fermi_Hubbard_chain(t0, U):
#     """
#     Construct the spin-1 Heisenberg Hamiltonian for given couplings.
    
#         Parameters
#         ----------
#         Jx : float
#             Coupling strength in x direction
#         Jy : float
#             Coupling strength in y direction
#         Jy : float
#             Coupling strength in z direction
#         hz : float
#             Coupling for Sz terms

#         Returns
#         -------
#         h : np.array (3, 3, 3, 3)
#             Spin-1 Heisenberg Hamiltonian.
            
            
#     """
#     d = 4


#     id_=np.eye(d)
#     h_0 = - t0*(np.kron(np.conj(cup).T@F,cup) - np.kron(cup@F,np.conj(cup).T) + np.kron(np.conj(cdo).T@F,cdo)\
#                 - np.kron(cdo@F,np.conj(cdo).T)) + 0.5*(U*np.kron(nup2@ndo2,id_) + U*np.kron(id_,nup2@ndo2)) 

#     # making local dimension two sites=4 for better convergence
#     #h_0=-t0*np.kron(Ekin_mat, np.eye(d))#-t0*np.kron( np.eye(d),Ekin_mat)+\
#         #0.5*U*np.kron((np.kron(density-0.5*sI,density-0.5*sI)),np.eye(4))+\
#         #    0.5*U*np.kron(np.eye(4),(np.kron(density-0.5*sI,density-0.5*sI)))

#     return h_0.reshape(d, d, d, d)

# def Fermi_Hubbard_chain_2(t0, U,mu,h):
#     """
#     Construct the spin-1 Heisenberg Hamiltonian for given couplings.
    
#         Parameters
#         ----------
#         Jx : float
#             Coupling strength in x direction
#         Jy : float
#             Coupling strength in y direction
#         Jy : float
#             Coupling strength in z direction
#         hz : float
#             Coupling for Sz terms

#         Returns
#         -------
#         h : np.array (3, 3, 3, 3)
#             Spin-1 Heisenberg Hamiltonian.
            
            
#     """
#     d = 4


#     id_=np.eye(d)
#     h_0 = - t0*(np.kron(np.conj(cup).T@F,cup) - np.kron(cup@F,np.conj(cup).T) + np.kron(np.conj(cdo).T@F,cdo)\
#                 - np.kron(cdo@F,np.conj(cdo).T)) + 0.5*(U*np.kron(nup@ndo,id_) + U*np.kron(id_,nup@ndo))+ \
#             0.5*(mu*np.kron(nup,id_) + mu*np.kron(id_,nup)+mu*np.kron(ndo,id_) + mu*np.kron(id_,ndo))+\
#                 0.5*(h*np.kron(nup-ndo,id_) + -0.5*h*np.kron(id_,nup-ndo))#+\
#               #  0.5*(h*np.kron(nup,id_)+h*np.kron(id_,nup) + -h*np.kron(id_,ndo)-h*np.kron(ndo,id_))

            
                                                      
                                                       

#     # making local dimension two sites=4 for better convergence
#     #h_0=-t0*np.kron(Ekin_mat, np.eye(d))#-t0*np.kron( np.eye(d),Ekin_mat)+\
#         #0.5*U*np.kron((np.kron(density-0.5*sI,density-0.5*sI)),np.eye(4))+\
#         #    0.5*U*np.kron(np.eye(4),(np.kron(density-0.5*sI,density-0.5*sI)))

#     return h_0.reshape(d, d, d, d)



# def CHSH():
#     """
#     Construct the spin-1 Heisenberg Hamiltonian for given couplings.
    
#         Parameters
#         ----------
   

#         Returns
#         -------
#         h : np.array (3, 3, 3, 3)
#             Spin-1 Heisenberg Hamiltonian.
#     """
#     d = 2
#     Q=np.kron(2*sZ, sI)
#     R=np.kron(sI, 2*sX)
#     S=(-np.kron(sI, 2*sZ)-np.kron(sI, 2*sX))/np.sqrt(2)
#     T=(np.kron(sI, 2*sZ)-np.kron(sI, 2*sX))/np.sqrt(2)
#     h0=Q@S+R@S+R@T-Q@T
#     # diviging by factor 2/4 going from pauli to spin
#     #h_temp = (np.real(J*np.kron(sX, sX) + J*np.kron(sY, sY) + Delta*np.kron(sZ, sZ)+0.5*h*np.kron(sI, sZ)+0.5*h*np.kron(sZ, sI)))
#    # h0 = (0.5 * np.kron(h_temp, np.eye(4)) + 0.5 * np.kron(np.eye(4), h_temp) +
#     #    np.kron(np.eye(2), np.kron(h_temp, np.eye(2)))).reshape(d, d, d, d)
#     # making local dimension two sites=4 for better convergenceSB7WJy
#     return h0.reshape(d,d,d,d)
