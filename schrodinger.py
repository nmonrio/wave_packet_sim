import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import constants as cts

class wave_packet(object):
    """
    Class which implements a numerical solution of the time-dependent
    Schrodinger equation for an arbitrary potential
    """
    def __init__(self, t_f, dt, dx, x_min, x_max, k_0, sigma_0, x_0, L, V_0, M):
        """
        Constructor method for the schrodinger wave packet. Stores the presets
        of step 1.
        """

        # Inputs from table 1
        self.t_final = t_f
        self.dt = dt
        self.dx = dx          
        self.x_min = x_min
        self.x_max = x_max
        self.k_not = k_0
        self.sigma_not = sigma_0
        self.x_0 = x_0
        self.L = L
        self.V_0 = V_0
        self.M = M

        # Necessary variables for potential matrices. In Python the imaginary unit is
        # denoted as "j" and needs to be multiplied by some number to act as such.
        self.alpha = self.dt / (4 * (self.dx**2))
        self.beta = 1 + 2 * 1j * self.alpha
        self.beta_conjugate = 1 - 2 * 1j * self.alpha
          
    
    def psi_initial_state(self, x):
        """
        Method for computing the schrodinger equation for a given x.
        """
        psi = (1/(np.pi * self.sigma_not**2))**(1/4) * np.exp(1j*self.k_not*x) * np.exp(-((x-self.x_0)**2) / (2*self.sigma_not**2))
        return psi
    
    def potential_barrier(self):
        """
        Method for creating the potential barrier
        """
        potential_barrier = []     

        for j in range(1, self.M + 1):
            
            x_j = self.x_min + (j - 1)*self.dx
            
            if 0 <= x_j <= 2:
                potential_barrier.append(2)
            else:
                potential_barrier.append(0)
        
        return potential_barrier

    def compute_rmat(self, rmat_diag, potential):
        """
        Method for computing the main diagonal of the potential matrix R.
        """
        for i in range(0, self.M):
            r_i = self.beta_conjugate - 1j * self.dt * (potential[i]/2)   # Each of the elements of the diagonal
            rmat_diag.append(r_i)
        
        return rmat_diag

    def compute_rmat_subdiag(self):
        """
        Method for computing the upper subdiagonal for M.
        """
        r_subdiag_entry = 1j * self.alpha
        
        return r_subdiag_entry

    def compute_lmat_(self, lmat_diag, potential):
        """
        Method for computing the main diagonal of the potential matrix L. 
        Similar to the method compute_rmat
        """
        for i in range(0, self.M):
            l_i = self.beta + 1j * self.dt * (potential[i]/2)   # Each of the elements of the diagonal
            lmat_diag.append(l_i)
        
        return lmat_diag

    def compute_lmat_subdiag(self):
        """
        Method for computing the upper subdiagonal for L.
        """
        l_subdiag_entry = - 1j * self.alpha

        return l_subdiag_entry

    def compute_strid(self, rmat_diag, rmat_subdiag, psi_current):
        """
        Method that computes STRID, the multiplication of R and PSI_CURRENT,
        as specified in step 5.
        """
        s_trid = [] 

        for i in range(0, self.M):
            if i == 0:
                entry = rmat_diag[i]*psi_current[i] + rmat_subdiag*psi_current[i+1]
                s_trid.append(entry)
            elif i == (self.M - 1):
                entry = rmat_subdiag*psi_current[i-1] + rmat_diag[i]*psi_current[i]
                s_trid.append(entry)
            else:
                entry = rmat_subdiag*psi_current[i-1] + rmat_diag[i]*psi_current[i] + rmat_subdiag*psi_current[i+1]
                s_trid.append(entry)
        
        return s_trid
    
    def compute_atrid(self, lmat_diag, lmat_subdiag):
        """
        Method for computing the ATRID matrix using method from guide. (Eq. 20)
        """
        a_trid = []

        for i in range(0, (self.M - 1)):
            if i == 0:   # Equivalent to i = 1 in Eq 20
                a_trid.append(lmat_subdiag / lmat_diag[i])
            else:
                a_trid.append(lmat_subdiag / (lmat_diag[i] - lmat_subdiag*a_trid[i-1]))
        return a_trid

    def overwrite_strid(self, s_trid, a_trid, lmat_diag, lmar_subdiag):
        """
        Method that implements equation 21 from guide, overwriting S_TRID into a clean version.
        """
        s_trid_prime = []

        for i in range(0, self.M):
            if i == 0:
                s_trid_prime.append(s_trid[i] / lmat_diag[i])
            else:
                s_trid_prime.append((s_trid[i] - lmar_subdiag*s_trid_prime[i-1])/(lmat_diag[i] - lmar_subdiag*a_trid[i-1]))

        return s_trid_prime

    def compute_psi_next(self, s_trid_prime, a_trid):
        """
        Method that computes PSI_NEXT array. Last part from step 5.
        """
        psi_next = []

        for i in range(self.M - 1, -1, -1):
            if (i == self.M - 1):
                psi_next.append(s_trid_prime[i])
            else: 
                psi_next.append(s_trid_prime[k] - a_trid[k]*psi_next[self.M - i -2])
        
        # overwrite psi_next:
        psi_next_good = []
        for k in range(self.M - 1,-1,-1):
            psi_next_good.append(psi_next[k])
    
        return psi_next
        

    


    













































# class Potential(object):
#     """
#     Class that implements any arbitrary potential
#     """
    
#     def theta(x):
#         """
#         theta function :
#         returns 0 if x<=0, and 1 if x>0
#         """
#         x = np.asarray(x)
#         y = np.zeros(x.shape)
#         y[x > 0] = 1.0
#         return y

#     def square_barrier(x, width, height):
#         """
#         Method that represents the square barrier.
#         """
#         return height * (theta(x) - theta(x - width))