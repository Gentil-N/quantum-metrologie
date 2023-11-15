import qutip as qp
import cavity as cty
from matplotlib import pyplot as plt
import scipy.constants as const
import numpy as np

ATOM_COUNT = 1
PHOTON_CAPACITY = 1

#simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 1.0, 1.0, 1.0, False, False)
simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 1.0, 1.0, 1.0, True, True) # freq 7000, g 73, k 4, r 4
#simple_cav = cty.Cavity(434.0 * const.tera, 434.0 * const.tera, ATOM_COUNT, PHOTON_CAPACITY, 820.0, np.sqrt(800.0 * const.kilo), np.sqrt(7.5 * const.kilo), True, True)

psi_init = qp.tensor(qp.spin_state(ATOM_COUNT/2, ATOM_COUNT/2), qp.fock(PHOTON_CAPACITY + 1, 0))
#psi_init = qp.tensor(qp.spin_state(ATOM_COUNT/2, -ATOM_COUNT/2), qp.coherent(PHOTON_CAPACITY + 1, 1))

#TOTAL_TIME = 100 * 10 * const.milli * const.nano #20 * const.milli
#NUM_DIV = 200 # = const.giga
#STEP = TOTAL_TIME / NUM_DIV
#
#last_state = psi_init
#total_time_range = []
#total_time_range.extend(simple_cav.run_simulation(psi_init, 0, STEP, 10000))
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

time_range = simple_cav.run_simulation(psi_init, 0, 50, 1000)
#time_range = simple_cav.run_simulation(psi_init, 0, 0.02, 1000)
#time_range = simple_cav.run_simulation(psi_init, 0, 0.01, 10000000)
#last_state = simple_cav.get_state()[-1]
#time_range = simple_cav.run_simulation(last_state, 10 * const.milli * const.nano, 10 * const.milli * const.nano, 10000)

print("computing g2...")

g2 = simple_cav.compute_g2()
#g2_qutip = simple_cav.compute_g2_qutip(psi_init, time_range)
g2_perso = simple_cav.compute_g2_perso(psi_init, time_range)
#print(g2_perso)

print("computing g1...")

g1 = simple_cav.compute_g1()
#g1_qutip = simple_cav.compute_g1_qutip(psi_init, time_range)
g1_perso = simple_cav.compute_g1_perso(psi_init, time_range)
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
ax4.plot(time_range[1::], np.real(g2)[1::], label="g²(t) real part")
ax4.plot(time_range[1::], np.imag(g2)[1::], label="g²(t) imaginary part")
#ax4.plot(time_range, np.real(g2_qutip), label="g²(t) qutip real part")
#ax4.plot(time_range, np.imag(g2_qutip), label="g²(t) qutip imaginary part")
ax4.plot(time_range, np.real(g2_perso), label="g²(t) perso real part", linestyle="dashed", color="gold")
ax4.plot(time_range, np.imag(g2_perso), label="g²(t) perso imaginary part", linestyle="dashed")
ax4.set(xlabel="time", ylabel="")
ax4.legend()

fig2 = plt.figure(num=2)
ax2 = fig2.subplots(nrows=1, ncols=1)
ax2.set_title(label=common_title)
ax2.plot(time_range[1::], np.real(g1)[1::], label="g¹(t) real part")
ax2.plot(time_range[1::], np.imag(g1)[1::], label="g¹(t) imaginary part")
#ax2.plot(time_range, np.real(g1_qutip), label="g¹(t) qutip real part")
#ax2.plot(time_range, np.imag(g1_qutip), label="g¹(t) qutip imaginary part")
ax2.plot(time_range, np.real(g1_perso), label="g¹(t) perso real part", linestyle="dashed", color="gold")
ax2.plot(time_range, np.imag(g1_perso), label="g¹(t) perso imaginary part", linestyle="dashed")
ax2.set(xlabel="time", ylabel="")
ax2.legend()

plt.show()