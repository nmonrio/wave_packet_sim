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
    sigma_not = 3
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
    t_exp = [0.0005184753635453, 0.0037387418098775, 0.0320168043340280, 0.2961256362742885, 0.8571273717728469, 0.9393059097791561, 0.9827213815087718, 0.9903886275799959, 0.9905519280703766, 0.9981519629612955, 0.9967316372243031, 0.9975524449641112, 0.9994550887034639, 0.9986515619552828, 0.9990866774819243, 0.9997635488213144]              # Array of simptotic values obtained in last excercise.
    for k in np.linspace(0, 8, 16):
        exp_wave = wave_packet(t_final, dt, dx, x_min, x_max, k, sigma_not, x_not, l, v_initial, m)
        e_expe.append(exp_wave.E_0)

    results.plot_trans_vs_energy(T_theoretical, energy, t_exp, e_expe)


    