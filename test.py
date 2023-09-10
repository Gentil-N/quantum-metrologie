import qutip as qp
import cavity as cty
from matplotlib import pyplot as plt
import numpy as np

PHOTON_CAPACITY = 10

simple_cav = cty.Cavity(1.0, 1.0, 5, PHOTON_CAPACITY, 10.0, 3.0, 1.0, True, True)

psi_init = qp.tensor(qp.spin_state(5/2, 5/2), qp.fock(PHOTON_CAPACITY + 1, 0))

time_range = simple_cav.run_simulation(psi_init, 0, 5, 1000)

#probs = simple_cav.get_all_probabilities()
g2 = simple_cav.compute_g2()

###############
### PLOTING ###
###############

common_title = "Simple Cavity (1 atom - 1 excitation) start |e, 0>"

#fig0 = plt.figure(num=0)
#ax0 = fig0.subplots(nrows=1, ncols=1)
#ax0.set_title(label=common_title)
#for prob in probs:
#    ax0.plot(time_range, prob[0], label= prob[1])
#ax0.set(xlabel="time", ylabel="probabilities")
#ax0.legend()

fig4 = plt.figure(num=4)
ax4 = fig4.subplots(nrows=1, ncols=1)
ax4.set_title(label="Simple Cavity (1 atom - 1 photon) start |e, 0>")
ax4.plot(time_range, np.real(g2), label="g²(t) real part")
ax4.plot(time_range, np.imag(g2), label="g²(t) imaginary part")
ax4.set(xlabel="time", ylabel="probabilities")
ax4.legend()

plt.show()