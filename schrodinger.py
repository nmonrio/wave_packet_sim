import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

import constants as cts

class Schrodinger(object):
    """
    Class which implements a numerical solution of the time-dependent
    Schrodinger equation for an arbitrary potential
    """
    def __init__(self, x, psi_x0, V_x, k0 = None, hbar=1, m=1, t0=0.0):

         # Validation of array inputs
        self.x, psi_x0, self.V_x = map(np.asarray, (x, psi_x0, V_x))
        N = self.x.size
        assert self.x.shape == (N,)
        assert psi_x0.shape == (N,)
        assert self.V_x.shape == (N,)

        # Set internal parameters
        self.hbar = hbar
        self.m = m
        self.t = t0
        self.dt_ = None
        self.N = len(x)
        self.dx = self.x[1] - self.x[0]
        self.dk = 2 * np.pi / (self.N * self.dx)

        # set momentum scale
        if k0 == None:
            self.k0 = -0.5 * self.N * self.dk
        else:
            self.k0 = k0
        self.k = self.k0 + self.dk * np.arange(self.N)

        self.psi_x = psi_x0
        self.compute_k_from_x()

        # variables which hold steps in evolution of the
        self.x_evolve_half = None
        self.x_evolve = None
        self.k_evolve = None

        # attributes used for dynamic plotting
        self.psi_x_line = None
        self.psi_k_line = None
        self.V_x_line = None

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