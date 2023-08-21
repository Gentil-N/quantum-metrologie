import qutip as qx
from matplotlib import pyplot as plt
import numpy as np
import math

OMEGA_RABI = 1.0
OMEGA = 1.0
HBAR = 1.0
HILBERT_DIM_PHOTON = 4
KAPPA = 1.0

def tens_from_fock(f_op):
    return qx.tensor(qx.qeye(2), f_op)

def tens_from_spin(s_op):
    return qx.tensor(s_op, qx.qeye(HILBERT_DIM_PHOTON))

op_s_pos = qx.sigmax() + 1j * qx.sigmay()
op_s_neg = qx.sigmax() - 1j * qx.sigmay()
op_sz = qx.sigmaz()

op_a = qx.destroy(HILBERT_DIM_PHOTON)
op_ad = qx.create(HILBERT_DIM_PHOTON)
op_n = op_ad * op_a

#print("1 ", tens_from_spin(op_sz))
#print("2 ", tens_from_fock(op_ad * op_a))
#print("3 ", qx.tensor(op_s_pos, op_a))
#print("4 ", qx.tensor(op_s_neg, op_ad))
op_hamiltonian = HBAR * OMEGA / 2.0 * tens_from_spin(op_sz) + HBAR * OMEGA * tens_from_fock(op_ad * op_a) + HBAR * KAPPA / 2.0 * (qx.tensor(op_s_pos, op_a) + qx.tensor(op_s_neg, op_ad))

psi_init = qx.tensor(qx.fock(HILBERT_DIM_PHOTON, 1), qx.spin_state(1/2, 1/2))
print("Psi init: ", psi_init)
time_range = np.linspace(0, 10, 100)

result = qx.mesolve(op_hamiltonian, psi_init, time_range, [], [])
res_n = qx.expect(qx.tensor(op_n, qx.qeye(2)), result.states)
res_n_test = []
print(result.states[0])
state_mesure = qx.tensor(qx.fock(HILBERT_DIM_PHOTON, 0), qx.spin_state(1/2, 1/2))
for state in result.states:
    res_n_test.append((psi_init.dag() * state * state.dag() * psi_init).tr())

p_e = np.cos(OMEGA*time_range / 2)**2

fig = plt.figure(num=0)
ax = fig.subplots(nrows=1, ncols=1)
ax.set_title(label="LASER !!!")
ax.plot(time_range, res_n_test)
ax.plot(time_range, p_e, marker='o')
ax.set(xlabel="time", ylabel="mean N")
plt.show()