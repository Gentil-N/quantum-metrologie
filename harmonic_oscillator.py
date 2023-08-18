import qutip as qx
from matplotlib import pyplot as plt
import numpy as np
import math
import scipy.constants as cs

MASS = 1.0
FREQ = 1.0
HBAR = 1.0
KAPPA = 0.1
N_STATE = 10

op_a = qx.destroy(N_STATE)
op_ad = qx.create(N_STATE)

op_N = op_ad * op_a

#op_x = (op_ad + op_a) * math.sqrt(HBAR/ (2 * MASS * FREQ))
#op_p = (op_ad - op_a) * 1j * math.sqrt(HBAR * MASS * FREQ / 2)

#op_H = HBAR * FREQ * (op_N + 1/2)
#
#time_range = np.linspace(0, 10, 3)
#psi_init = qx.fock(N_STATE, 0)
#
#result = qx.mesolve(op_H, psi_init, time_range)
#print(result.states)
#print(qx.expect(op_N, result.states))
#exit()

energy_list = []
for i in range(0, N_STATE):
    state = qx.fock(N_STATE, i)
    energy_list.append((state.dag() * op_H * state).tr())

pos = [*range(1, len(energy_list) + 1)]

plt.bar(pos, energy_list)
plt.show()