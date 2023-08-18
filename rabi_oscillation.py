import qutip as qx
from matplotlib import pyplot as plt
import numpy as np
import math

GAMMA_B = 1.0

psi_init = qx.Qobj([[1], [1]])
psi_init *= 1 / math.sqrt(2)

dens_plusx = psi_init * psi_init.dag()

hamiltonian = - GAMMA_B * qx.sigmaz()

time_range = np.linspace(0, 2 * math.pi, 100)

result = qx.mesolve(hamiltonian, psi_init, time_range, [], [])

fig = plt.figure(num=0)
ax = fig.subplots(nrows=1, ncols=1)
ax.set_title(label="Rabi Oscillations")
ax.plot(time_range, qx.expect(dens_plusx, result.states))
ax.set(xlabel="time", ylabel="mean +X")
plt.show()