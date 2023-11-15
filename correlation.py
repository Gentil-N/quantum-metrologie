import qutip as qp
import cavity as cty
from matplotlib import pyplot as plt
import scipy.constants as const
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

ATOM_COUNT = 8
PHOTON_CAPACITY = 8

#simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 1.0, 1.0, 1.0, False, False)
#simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 1.0, 1.0, 1.0, True, True) # freq 7000, g 73, k 4, r 4
simple_cav = cty.Cavity(434.0 * const.tera, 434.0 * const.tera, ATOM_COUNT, PHOTON_CAPACITY, 820.0, np.sqrt(800.0 * const.kilo), np.sqrt(7.5 * const.kilo), True, True)

psi_init = qp.tensor(qp.spin_state(ATOM_COUNT/2, ATOM_COUNT/2), qp.fock(PHOTON_CAPACITY + 1, 0))

#time_range = simple_cav.run_simulation(psi_init, 0, 5, 100)
time_range = simple_cav.run_simulation(psi_init, 0, 0.002, 100)
exit()
#time_range = simple_cav.run_simulation(psi_init, 0, 0.01, 10000000)
#last_state = simple_cav.get_state()[-1]
#time_range = simple_cav.run_simulation(last_state, 10 * const.milli * const.nano, 10 * const.milli * const.nano, 10000)

print("computing g2...")

taulist = time_range

g2_array = simple_cav.compute_g2_array(psi_init, time_range, taulist)

#g2_array_qutip = simple_cav.compute_g2_array_qutip(psi_init, time_range, taulist)

###############
### PLOTING ###
###############

common_title = "Cavity (Dicke Model - Master Equation)"

fig0 = plt.figure(num=0)
ax0 = fig0.subplots(nrows=1, ncols=1)
ax0.set_title(common_title)
ax0.contourf(time_range, taulist, g2_array, cmap='inferno', levels=70)
ax0.set(xlabel="time (s)", ylabel="tau (s)")
fig0.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='inferno'), ax=ax0)

fig1 = plt.figure(num=1)
ax1 = fig1.subplots(nrows=1, ncols=1)
ax1.set_title(label=common_title)
#ax1.plot(time_range, g2_array[0], label="r0")
ax1.plot(time_range, g2_array[1], label="r1")
ax1.plot(time_range, g2_array[50], label="r2")
ax1.plot(time_range, g2_array[70], label="r3")
ax1.set(xlabel="time", ylabel="")
ax1.legend()

np.save('data_correlation_8_8.npy', g2_array)
#arr_loaded = np.load('myarray.npy')

#fig2 = plt.figure(num=2)
#ax2 = fig2.subplots(nrows=1, ncols=1)
#ax2.set_title(common_title + " Qutip")
#ax2.contourf(time_range, taulist, np.real(g2_array_qutip), cmap='inferno', levels=70)
#ax2.set(xlabel="time (s)", ylabel="tau (s)")
#fig2.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='inferno'), ax=ax2)
#
#fig3 = plt.figure(num=3)
#ax3 = fig3.subplots(nrows=1, ncols=1)
#ax3.set_title(label=common_title + " Qutip")
#ax1.plot(time_range, g2_array[0], label="r0")
#ax3.plot(time_range, np.real(g2_array_qutip[1]), label="r1")
#ax3.plot(time_range, np.real(g2_array_qutip[50]), label="r2")
#ax3.plot(time_range, np.real(g2_array_qutip[70]), label="r3")
#ax3.set(xlabel="time", ylabel="")
#ax3.legend()

plt.show()