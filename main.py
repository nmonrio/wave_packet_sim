import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

from methods import wave_packet
from methods import results 

if __name__ == "__main__":
    
    dt = 0.005   #1.- define parameters from table 1
    t_final = 50
    dx = 0.05
    x_min = -60
    x_max = 60
    k_not = 1
    sigma_not = 3
    x_not = -10
    l = 2
    v_initial = 2
    m = 2401

    wave = wave_packet(t_final, dt, x_min, x_max, k_not, sigma_not, x_not, l, v_initial, m)

    psi_current = []                    
    psi_next = []   
    lmat_diag = []
    lmat_subdiag = []
    rmat_diag = []
    rmat_subdiag = []
    potential = wave.potential_barrier()   
    
    for i in range(1, m+1):
        x_j = x_min + (j - 1)*dx
        entry = wave.psi_initial_state(x_j)
        psi_current.append(entry)

    lmat_diag = wave.compute_lmat(lmat_diag, potential)
    lmat_subdiag = wave.compute_lmat_subdiag()
    rmat_diag =  wave.compute_rmat(rmat_diag, potential)
    rmat_subdiag = wave.compute_rmat_subdiag()
    a_trid = wave.compute_atrid(lmat_diag, lmat_subdiag)
    psi_squared = []
    integral_list = []
    integral_list.append(wave.integral(psi_current))

    # Initialise time loop
    for t in np.arange(0, t_final, dt):
        