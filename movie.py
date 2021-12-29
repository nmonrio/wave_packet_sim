import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from src.schrodinger import wave_packet  
plt.rcParams['animation.ffmpeg_path'] = r'C:\Users\Usuario\Desktop\ffmpeg-master-latest-win64-gpl\bin\ffmpeg.exe'

# Time related contrains
dt = 0.005
t_final = 50

# Spatial related contrains
dx = 0.05
x_min = -60
x_max = 60

# To illustrate hypothesis we are variating parametres.
k_not = 1
k_not2 = 8
sigma_not = 3
sigma_not2 = 7

# The rest of necessary parametres
x_not = -10
l = 2
v_initial = 2
m = 2401

# Create the three wave packets
wave = wave_packet(t_final, dt, dx, x_min, x_max, k_not, sigma_not, x_not, l, v_initial, m)
wave2 = wave_packet(t_final, dt, dx, x_min, x_max, k_not2, sigma_not, x_not, l, v_initial, m)
wave3 = wave_packet(t_final, dt, dx, x_min, x_max, k_not, sigma_not2, x_not, l, v_initial, m)


x_values = wave.x_values()
t_values = wave.t_values()
potential = wave.potential_barrier()

psi_current, psi_current2, psi_current3 = [], [], []
psi_next, psi_next2, psi_next3 = [], [], []

lmat_diag, lmat_diag2, lmat_diag3 = [], [], []
lmat_subdiag, lmat_subdiag2, lmat_subdiag3 = [], [], []
rmat_diag, rmat_diag2, rmat_diag3 = [], [], []
rmat_subdiag, rmat_subdiag2, rmat_subdiag3 = [], [], []

for j in range(1, m + 1):
    x_j = x_min + (j - 1) * dx
    entry = wave.psi_initial_state(x_j)
    entry2 = wave2.psi_initial_state(x_j)
    entry3 = wave3.psi_initial_state(x_j)
    psi_current.append(entry)
    psi_current2.append(entry2)
    psi_current3.append(entry3)

# Creation of L & R matrices, constant through time for each of the cases
lmat_diag = wave.compute_lmat(lmat_diag, potential)
lmat_subdiag = wave.compute_lmat_subdiag()
rmat_diag = wave.compute_rmat(rmat_diag, potential)
rmat_subdiag = wave.compute_rmat_subdiag()

lmat_diag2 = wave2.compute_lmat(lmat_diag2, potential)
lmat_subdiag2 = wave2.compute_lmat_subdiag()
rmat_diag2 = wave2.compute_rmat(rmat_diag2, potential)
rmat_subdiag2 = wave2.compute_rmat_subdiag()

lmat_diag3 = wave3.compute_lmat(lmat_diag3, potential)
lmat_subdiag3 = wave3.compute_lmat_subdiag()
rmat_diag3 = wave3.compute_rmat(rmat_diag3, potential)
rmat_subdiag3 = wave3.compute_rmat_subdiag()

a_trid = wave.compute_atrid(lmat_diag, lmat_subdiag)
a_trid2 = wave2.compute_atrid(lmat_diag2, lmat_subdiag2)
a_trid3 = wave3.compute_atrid(lmat_diag3, lmat_subdiag3)

psi_squared, psi_squared2, psi_squared3  = [], [], []

# Probability data for plotting
psi_squared_data = [wave.compute_psi_squared(psi_current)]
psi_squared_data2 = [wave2.compute_psi_squared(psi_current2)]
psi_squared_data3 = [wave3.compute_psi_squared(psi_current3)]

# Initialise time loop
for t in np.arange(0, t_final, dt):
    s_trid = wave.compute_strid(rmat_diag, rmat_subdiag, psi_current)
    s_trid_prime = wave.overwrite_strid(s_trid, a_trid, lmat_diag, lmat_subdiag)
    s_trid2 = wave2.compute_strid(rmat_diag2, rmat_subdiag2, psi_current2)
    s_trid_prime2 = wave2.overwrite_strid(s_trid2, a_trid2, lmat_diag2, lmat_subdiag2)
    s_trid3 = wave3.compute_strid(rmat_diag3, rmat_subdiag3, psi_current3)
    s_trid_prime3 = wave.overwrite_strid(s_trid3, a_trid3, lmat_diag3, lmat_subdiag3)

    psi_next = wave.compute_pre_psi_next(s_trid_prime, a_trid)
    ow_psi_next = wave.overwrite_psi_next(psi_next)
    psi_next2 = wave2.compute_pre_psi_next(s_trid_prime2, a_trid2)
    ow_psi_next2 = wave2.overwrite_psi_next(psi_next2)
    psi_next3 = wave3.compute_pre_psi_next(s_trid_prime3, a_trid3)
    ow_psi_next3 = wave.overwrite_psi_next(psi_next3)

    psi_current = ow_psi_next
    psi_current2 = ow_psi_next2
    psi_current3 = ow_psi_next3

    psi_squared = wave.compute_psi_squared(psi_current)
    psi_squared_data.append(psi_squared)
    psi_squared2 = wave2.compute_psi_squared(psi_current2)
    psi_squared_data2.append(psi_squared2)
    psi_squared3 = wave3.compute_psi_squared(psi_current3)
    psi_squared_data3.append(psi_squared3)

    psi_squared = []
    ow_psi_next = []
    psi_next = []

    psi_squared2 = []
    ow_psi_next2 = []
    psi_next2 = []

    psi_squared3 = []
    ow_psi_next3 = []
    psi_next3 = []


fig = plt.figure()
# Plotting limits
xlim = (-25, 25)
ylim = (0, 0.5)

ax1 = fig.add_subplot(111, xlim=xlim, ylim=ylim)
psi_squared_line, = ax1.plot([], [], c='r', label='$k_0 = 1, \sigma_0 = 3$')
psi_squared_line2, = ax1.plot([], [], c='b', label='$k_0 = 8, \sigma_0 = 3$')
psi_squared_line3, = ax1.plot([], [], c='g', label='$k_0 = 1, \sigma_0 = 7$')
title = ax1.set_title("Quantum Physics: Numerical Lab")
ax1.legend(prop=dict(size=12))
ax1.set_xlabel('$x$')
ax1.set_ylabel('$|\psi(x)|^2$')

# Plot alongside the potential, stays the same all the time
vx = plt.twinx()
vx.plot(x_values, potential, color='black')
vx.set_ylim(0, 3)
vx.set_ylabel('$V(x)$')


def animate(i):

    psi_squared_line.set_data(x_values, psi_squared_data[i])
    psi_squared_line2.set_data(x_values, psi_squared_data2[i])
    psi_squared_line3.set_data(x_values, psi_squared_data3[i])

    return (psi_squared_line, psi_squared_line2, psi_squared_line3,)

anim = animation.FuncAnimation(fig, animate, frames=5000, interval=1, blit=True)



plt.show()

anim.save('animation.mp4', fps=120)

