import numpy as np
from src.schrodinger import wave_packet 
from src.analysis import results

if __name__ == "__main__":

    # Initialize values
    dt = 0.005
    t_final = 50
    dx = 0.05
    x_min = -60
    x_max = 60
    k_not = 1
    sigma_not = 7 
    x_not = -10
    l = 2
    v_initial = 2
    m = 2401
    wave = wave_packet(t_final, dt, dx, x_min, x_max, k_not, sigma_not, x_not, l, v_initial, m)

    # THEORETICAL VALUES
    # Energy that will be inputted inside of theoretical values 
    # for a smoother line.
    energy = np.linspace(0, 35, 200) 
    T_theoretical = wave.compute_theoretical_transc(energy)

    # WAVE PACKET VALUES
    e_expe = []
    t_exp = [0.0014848704537606275, 0.01849684533559284, 0.037395512991021486, 0.24361592929484127, 0.9562382873022899, 0.925730146377701,0.9911664168224473, 0.9894604146892002, 0.9899280002164821, 0.9995428657034344, 0.9961962340020798, 0.997532999909798, 0.9998478239348638, 0.9984310795077675, 0.9990986857580242, 0.9999216717027525 ]
    for k in np.linspace(0, 8, 16):
        exp_wave = wave_packet(t_final, dt, dx, x_min, x_max, k, sigma_not, x_not, l, v_initial, m)
        e_expe.append(exp_wave.E_0)

    # Graphs for part ii
    results.plot_trans_vs_energy(T_theoretical, energy, t_exp, e_expe)





    