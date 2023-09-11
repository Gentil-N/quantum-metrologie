import qutip as qp
import cavity as cty
from matplotlib import pyplot as plt
import numpy as np

ATOM_COUNT = 1
PHOTON_CAPACITY = 30

simple_cav = cty.Cavity(1.0, 1.0, ATOM_COUNT, PHOTON_CAPACITY, 10.0, 1.0, 4.0, True, True) # freq 7000, g 73, k 4, r 4

psi_init = qp.tensor(qp.spin_state(ATOM_COUNT/2, -ATOM_COUNT/2), qp.fock(PHOTON_CAPACITY + 1, 0))
#psi_init = qp.tensor(qp.spin_state(ATOM_COUNT/2, -ATOM_COUNT/2), qp.coherent(PHOTON_CAPACITY + 1, 1))

time_range = simple_cav.run_simulation(psi_init, 0, 5, 100)

print("computing g2...")

#probs = simple_cav.get_all_probabilities()
g2 = simple_cav.compute_g2()

print("computing g1...")

g1 = simple_cav.compute_g1()

###############
### PLOTING ###
###############

common_title = "Simple Cavity (1 atom)"

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
ax4.set(xlabel="time", ylabel="probabilities")
ax4.legend()

fig2 = plt.figure(num=2)
ax2 = fig2.subplots(nrows=1, ncols=1)
ax2.set_title(label=common_title)
ax2.plot(time_range[1::], np.real(g1)[1::], label="g¹(t) real part")
ax2.plot(time_range[1::], np.imag(g1)[1::], label="g¹(t) imaginary part")
ax2.set(xlabel="time", ylabel="probabilities")
ax2.legend()

plt.show()