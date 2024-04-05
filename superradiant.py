from matplotlib import pyplot as plt
import numpy as np
import qutip as qp
import math

FOCK_DIM = 100
HBAR = 1.0
DELTA_OMEGA = 0.0
G = 1.0
MEAN_N = 25

def normalize_state(state):
    state *= 1/np.sqrt((state.dag() * state).tr())
    return state

spin_up = qp.spin_state(1/2, 1/2)
spin_down = qp.spin_state(1/2, -1/2)

op_s_pos = 1/2 * (qp.sigmax() + 1j * qp.sigmay())
op_s_neg = 1/2 * (qp.sigmax() - 1j * qp.sigmay())
op_sz = qp.sigmaz()

op_a = qp.destroy(FOCK_DIM)
op_ad = qp.create(FOCK_DIM)
op_n = op_ad * op_a

spin_up = qp.spin_state(1/2, 1/2)
spin_down = qp.spin_state(1/2, -1/2)

op_hamiltonian = HBAR * DELTA_OMEGA * qp.tensor(qp.qeye(2), op_n) + HBAR * G * (qp.tensor(op_s_pos, op_a) + qp.tensor(op_s_neg, op_ad))

#psi_init = qp.tensor(spin_down, qp.coherent(FOCK_DIM, math.sqrt(MEAN_N)))
psi_init = qp.tensor(spin_up, qp.fock(FOCK_DIM, 0))

time_range = np.linspace(0, 100, 10000)

result = qp.mesolve(op_hamiltonian, psi_init, time_range, [], [])

ket_spin_up_fock_0 = qp.tensor(spin_up, qp.fock(FOCK_DIM, 0))
rho_spin_up_fock_0 = ket_spin_up_fock_0 * ket_spin_up_fock_0.dag()
ket_spin_down_fock_1 = qp.tensor(spin_down, qp.fock(FOCK_DIM, 1))
rho_spin_down_fock_1 = ket_spin_down_fock_1 * ket_spin_down_fock_1.dag()

ket_fock_1 = qp.tensor(spin_up, qp.qeye(FOCK_DIM))

fig1, ax1 = plt.subplots(1, 1)
fig1.set_size_inches(18.5, 10.5)
ax1.plot(time_range, np.cos(time_range)**2, label="model")
ax1.plot(time_range, qp.expect(rho_spin_up_fock_0, result.states), label="spin up")
ax1.plot(time_range, qp.expect(rho_spin_down_fock_1, result.states), label="spin down")
ax1.plot(time_range, qp.expect(ket_fock_1 * ket_fock_1.dag(), result.states), label="coherent revival")
ax1.legend()
ax1.grid()
plt.show()
