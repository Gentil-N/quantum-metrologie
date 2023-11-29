import qutip as qp
import cavity as cty
from matplotlib import pyplot as plt
import scipy.constants as const
import numpy as np
from cumulant.cavity_global import *

simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, g * np.sqrt(ATOM_COUNT), kappa, gamma, nu, True, True, True)

#simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 1.0, 1.0, 1.0, 0.0, False, False, False)
#simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 1.0, 1.0, 1.0, 1.0, True, True, True) # freq 7000, g 73, k 4, r 4
#simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 10.0, 40.0, 9.0, 1.0, True, True, True) # Superradiance from julia tutorial
#simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 10.0, 40.0, 9.0, 1.0, True, False, True) # Superradirance burst
#simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 1.5, 1.0, 0.25, 4.0, True, True, True)
#simple_cav = cty.Cavity(434.0 * const.tera, 434.0 * const.tera, ATOM_COUNT, PHOTON_CAPACITY, 820.0, 800.0 * const.kilo, 7.5 * const.kilo, 7.5 * const.kilo, True, True, True)

#TOTAL_TIME = 100 * 10 * const.milli * const.nano #20 * const.milli
#NUM_DIV = 200 # = const.giga
#STEP = TOTAL_TIME / NUM_DIV
#
#last_state = init_state
#total_time_range = []
#total_time_range.extend(simple_cav.run_simulation(init_state, 0, STEP, 10000))
#res = simple_cav.get_all_probabilities()
#probs = []
#for ind_prob in res:
#    probs.append([ind_prob[0].tolist(), ind_prob[1]])
#last_state = simple_cav.get_state()[-1]
#
#for i in range(1, int(NUM_DIV)):
#    total_time_range.extend(simple_cav.run_simulation(last_state, i * STEP, (i + 1) * STEP, 10000))
#    res = simple_cav.get_all_probabilities()
#    last_state = simple_cav.get_state()[-1]
#    for i in range(len(res)):
#        probs[i][0].extend(res[i][0].tolist())

#time_range = simple_cav.run_simulation(init_state, 0, 50, 1000)
time_range = simple_cav.run_simulation(init_state, 0, 10, 1000) # Superradiance from julia tutorial
#time_range = simple_cav.run_simulation(init_state, 0, 0.0001, 1000)
#time_range = simple_cav.run_simulation(init_state, 0, 0.01, 10000000)
#last_state = simple_cav.get_state()[-1]
#time_range = simple_cav.run_simulation(last_state, 10 * const.milli * const.nano, 10 * const.milli * const.nano, 10000)

steady_state = simple_cav.compute_steady_state()

### compute <n> <Sz>

#expect = simple_cav.compute_expectation_values(0)
#
#fig1 = plt.figure(num=4)
#ax1 = fig1.subplots(nrows=1, ncols=1)
#ax1.set_title(label="Expectation values")
#ax1.plot(time_range, expect[0], label="<n>")
#ax1.plot(time_range, expect[1], label="<Sz>")
#ax1.plot(time_range, expect[2], label="<S+>")
#ax1.set(xlabel="time", ylabel="")
#ax1.legend()
#
#plt.show()
#
#exit()

### compute <n> [end]

print("computing g2...")

#time_range = np.linspace(0.0, 0.0001, 1000)
time_range = np.linspace(0.0, 10, 1000)

#g2 = simple_cav.compute_g2()
g2_qutip = simple_cav.compute_g2_qutip(steady_state, time_range)
g2_perso = simple_cav.compute_g2_perso(steady_state, time_range)
#print(g2_perso)

print("computing g1...")

#g1 = simple_cav.compute_g1()
g1_qutip = simple_cav.compute_g1_qutip(steady_state, time_range)
g1_perso = simple_cav.compute_g1_perso(steady_state, time_range)
#print(g1_perso)

print("computing probabilities...")

#probs = simple_cav.get_all_probabilities()

###############
### PLOTING ###
###############

common_title = "Cavity (Dicke Model - Master Equation)"

#fig0 = plt.figure(num=0)
#ax0 = fig0.subplots(nrows=1, ncols=1)
#ax0.set_title(label=common_title)
#for prob in probs:
#    ax0.plot(total_time_range, prob[0], label=prob[1])
#ax0.set(xlabel="time", ylabel="probabilities")
#ax0.legend()

#fig0 = plt.figure(num=0)
#ax0 = fig0.subplots(nrows=1, ncols=1)
#ax0.set_title(label=common_title)
#for prob in probs:
#    ax0.plot(time_range, prob[0], label= prob[1])
#ax0.set(xlabel="time", ylabel="probabilities")
#ax0.legend()

fig4 = plt.figure(num=4)
ax4 = fig4.subplots(nrows=1, ncols=1)
ax4.set_title(label=common_title)
#ax4.plot(time_range[1::], np.real(g2)[1::], label="g²(t) repro real part")
#ax4.plot(time_range[1::], np.imag(g2)[1::], label="g²(t) repro imaginary part")
ax4.plot(time_range, np.real(g2_qutip), label="g²(t) real part")
ax4.plot(time_range, np.imag(g2_qutip), label="g²(t) imaginary part")
ax4.plot(time_range, np.real(g2_perso), label="g²(t) perso real part", linestyle="dashed", color="gold")
ax4.plot(time_range, np.imag(g2_perso), label="g²(t) perso imaginary part", linestyle="dashed")
ax4.set(xlabel="time", ylabel="")
ax4.legend()

fig2 = plt.figure(num=2)
ax2 = fig2.subplots(nrows=1, ncols=1)
ax2.set_title(label=common_title)
#ax2.plot(time_range[1::], np.real(g1)[1::], label="g¹(t) repro real part")
#ax2.plot(time_range[1::], np.imag(g1)[1::], label="g¹(t) repro imaginary part")
ax2.plot(time_range, np.real(g1_qutip), label="g¹(t) real part")
ax2.plot(time_range, np.imag(g1_qutip), label="g¹(t) imaginary part")
ax2.plot(time_range, np.real(g1_perso), label="g¹(t) perso real part", linestyle="dashed", color="gold")
ax2.plot(time_range, np.imag(g1_perso), label="g¹(t) perso imaginary part", linestyle="dashed")
ax2.set(xlabel="time", ylabel="")
ax2.legend()

plt.show()