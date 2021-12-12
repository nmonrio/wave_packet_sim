import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Helvetica"]})
from src.methods import wave_packet         # Related to computations and algorythm
from src.methods import results             # Related to the outputting of results

if __name__ == "__main__":
    
    # ALGORYTHM IMPLEMENTATION
    # We define parametres from table and create a wave packet with such characteristics.
    dt = 0.005   
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
    wave = wave_packet(t_final, dt, dx, x_min, x_max, k_not, sigma_not, x_not, l, v_initial, m)

    x_values = wave.x_values()
    t_values = wave.t_values()
    potential = wave.potential_barrier()

    psi_current = []                    
    psi_next = []   
    lmat_diag = []
    lmat_subdiag = []
    rmat_diag = []
    rmat_subdiag = []

    
    for j in range(1, m+1):
        x_j = x_min + (j - 1)*dx
        entry = wave.psi_initial_state(x_j)
        psi_current.append(entry)


    # Creation of L & R matrices, constant through time
    lmat_diag = wave.compute_lmat(lmat_diag, potential)
    lmat_subdiag = wave.compute_lmat_subdiag()
    rmat_diag =  wave.compute_rmat(rmat_diag, potential)
    rmat_subdiag = wave.compute_rmat_subdiag()
    
    a_trid = wave.compute_atrid(lmat_diag, lmat_subdiag)
    psi_squared = []

    # Integral data for plotting.
    integral_data = []
    integral_data.append(wave.integral(psi_current))

    # Probability data for plotting
    psi_squared_data = []
    psi_squared_data.append(wave.compute_psi_squared(psi_current))

    # Initialise transmission coefficient array for plotting
    trans_coeff = []

    # Initialise time loop
    for t in np.arange(0, t_final, dt):
        s_trid = wave.compute_strid(rmat_diag, rmat_subdiag, psi_current)
        s_trid_prime = wave.overwrite_strid(s_trid, a_trid, lmat_diag, lmat_subdiag)
        psi_next = wave.compute_pre_psi_next(s_trid_prime, a_trid)
        ow_psi_next = wave.overwrite_psi_next(psi_next)
        psi_current = ow_psi_next
        trans_coeff.append(wave.compute_transmission_coeff(psi_current))
        # integral_data.append(wave.integral(psi_current))

        if t == 500*dt or t == 1000*dt or t == 1500*dt or t == 2000*dt:
            psi_squared = wave.compute_psi_squared(psi_current)
            psi_squared_data.append(psi_squared)

        psi_squared = []
        ow_psi_next = []
        psi_next = []


    # ANALYSIS OF RESULTS

    # Storing as text files. (Uncomment to use)
    # results.trans_coeff_2_txt(trans_coeff)
    # results.integral_2_txt(integral_data)
    # results.psi_squared_2_txt(psi_squared_data)

    # Exercise 1
    results.psi_squared_2_plot(psi_squared_data, x_values, potential)

    # Exercise 2
    asymptote = wave.trans_coef_asymptote(trans_coeff)
    print(f'Transmission coefficient\'s asymptote: {asymptote}')
    results.plot_trans_coefficient(trans_coeff, t_values, asymptote)
    
    # Exercise 3

    