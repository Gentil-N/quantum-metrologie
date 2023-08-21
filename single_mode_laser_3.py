import qutip as qx
from matplotlib import pyplot as plt
import numpy as np
import math

OMEGA_RABI = 1.0
OMEGA = 1.0
HBAR = 1.0
HILBERT_DIM_PHOTON = 4
KAPPA = 10.0

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

spin_up = qx.spin_state(1/2, 1/2)
spin_down = qx.spin_state(1/2, -1/2)

#print("1 ", tens_from_spin(op_sz))
#print("2 ", tens_from_fock(op_ad * op_a))
#print("3 ", qx.tensor(op_s_pos, op_a))
#print("4 ", qx.tensor(op_s_neg, op_ad))
op_hamiltonian = HBAR * OMEGA / 2.0 * tens_from_spin(op_sz) + HBAR * OMEGA * tens_from_fock(op_ad * op_a) + HBAR * KAPPA / 2.0 * (qx.tensor(op_s_pos, op_a) + qx.tensor(op_s_neg, op_ad))

psi_init = qx.tensor(spin_up, qx.fock(HILBERT_DIM_PHOTON, 0)) # = |e,0>
print("Psi init: ", psi_init)
time_range = np.linspace(0, 1, 100)

result = qx.mesolve(op_hamiltonian, psi_init, time_range, [], [])

first_mesure = qx.tensor(spin_down, qx.fock(HILBERT_DIM_PHOTON, 1)) # |g,1>
second_mesure = qx.tensor(spin_up, qx.fock(HILBERT_DIM_PHOTON, 0)) # |e,0>

res_first = qx.expect(first_mesure * first_mesure.dag(), result.states)
res_second = qx.expect(second_mesure * second_mesure.dag(), result.states)

op_n_all = tens_from_spin(spin_up * spin_up.dag()) + tens_from_fock(op_n)

res_n = qx.expect(op_n_all, result.states)

fig0 = plt.figure(num=0)
ax0 = fig0.subplots(nrows=1, ncols=1)
ax0.set_title(label="LASER !!!")
ax0.plot(time_range, res_first)
ax0.plot(time_range, res_second)
ax0.set(xlabel="time", ylabel="prob")

fig1 = plt.figure(num=1)
ax1 = fig1.subplots(nrows=1, ncols=1)
ax1.set_title(label="LASER !!!")
ax1.plot(time_range, res_n)

plt.show()