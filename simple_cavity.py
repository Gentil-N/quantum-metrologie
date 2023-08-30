import qutip as qx
from matplotlib import pyplot as plt
import numpy as np
import math

OMEGA_PHOTONS = 1.0
OMEGA_ATOM = 1.0
HBAR = 1.0
HILBERT_DIM_PHOTON = 2
G = 10.0
KAPPA = 3.0
REPUMP = 1.0
ENABLE_DECAY = True
ENABLE_REPUMP = True

def create_coherent_state(fock_space_dim, alpha):
    return np.exp(-np.abs(alpha)**2) *  sum((alpha**n)/np.sqrt(math.factorial(n))*qx.fock(fock_space_dim, n) for n in range(0, fock_space_dim))

def normalize_state(state):
    state *= 1/np.sqrt((state.dag() * state).tr())
    return state

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

op_collapsing = []
if ENABLE_DECAY:
    op_collapsing.append(KAPPA * tens_from_fock(op_a)) 
if ENABLE_REPUMP:
    op_collapsing.append(REPUMP * tens_from_spin(op_s_pos))

op_hamiltonian = HBAR * OMEGA_ATOM / 2.0 * tens_from_spin(op_sz) + HBAR * OMEGA_PHOTONS * tens_from_fock(op_ad * op_a) + HBAR * G / 2.0 * (qx.tensor(op_s_pos, op_a) + qx.tensor(op_s_neg, op_ad))

psi_init = qx.tensor(spin_up, qx.fock(HILBERT_DIM_PHOTON, 0)) # = |e,0>

time_range = np.linspace(0, 10, 1000)

result = qx.mesolve(op_hamiltonian, psi_init, time_range, op_collapsing, [])

#print(create_coherent_state(HILBERT_DIM_PHOTON, COHERENT_ALPHA))
#print(normalize_state(create_coherent_state(2, 1)))

first_mesure = qx.tensor(spin_down, qx.fock(HILBERT_DIM_PHOTON, 1)) # |g,1>
second_mesure = qx.tensor(spin_up, qx.fock(HILBERT_DIM_PHOTON, 0)) # |e,0>
#third_mesure = psi_init

res_first = qx.expect(first_mesure * first_mesure.dag(), result.states)

res_second = qx.expect(second_mesure * second_mesure.dag(), result.states)

op_n_all = tens_from_fock(op_n)
print(result.states[0])
print(result.states[0].dag() * op_n_all * result.states[0])
res_n_mesolve = qx.mesolve(op_hamiltonian, psi_init, time_range, op_collapsing, [op_n_all]).expect[0]
res_n_mesolve = np.array(res_n_mesolve)

res_a_test = []
i = 0
for t in time_range:
    res_a_test.append((result.states[i].dag() * tens_from_fock(op_ad) * np.exp(-1j * OMEGA_PHOTONS * t) * result.states[i]).tr())
    i += 1

G1 = qx.correlation_2op_2t(op_hamiltonian, psi_init, None, time_range, op_collapsing, tens_from_fock(op_ad), tens_from_fock(op_a))
print("**********************", G1)
if len(G1.shape) == 2:
    G1 = G1[0]
print("**********************", G1)
n = np.array(res_n_mesolve)
g1 = G1 / np.sqrt(n[0] * n)
print("**********************", g1)

#res_third = qx.expect(third_mesure * third_mesure.dag(), result.states)

res_n = qx.expect(op_n_all, result.states)

fig0 = plt.figure(num=0)
ax0 = fig0.subplots(nrows=1, ncols=1)
ax0.set_title(label="Simple Cavity (1 atom - 1 photon) start |e, 0>")
ax0.plot(time_range, res_first, label = "|<g,1|ψ(t)>|²")
ax0.plot(time_range, res_second, label= "|<e,0|ψ(t)>|²")
ax0.set(xlabel="time", ylabel="probabilities")
ax0.legend()

fig1 = plt.figure(num=1)
ax1 = fig1.subplots(nrows=1, ncols=1)
ax1.set_title(label="Simple Cavity (1 atom - 1 photon) start |e, 0>")
ax1.plot(time_range, res_n, label="|<ψ(t)|N|ψ(t)>|²")
ax1.plot(time_range, res_n_mesolve, label="direct mesolve: |<ψ(t)|N|ψ(t)>|²")
ax1.plot(time_range, np.imag(res_a_test).tolist(), label="test!!!")
ax1.set(xlabel="time", ylabel="probabilities")
ax1.legend()

fig2 = plt.figure(num=2)
ax2 = fig2.subplots(nrows=1, ncols=1)
ax2.set_title(label="Simple Cavity (1 atom - 1 photon) start |e, 0>")
ax2.plot(time_range, np.real(g1), label="g¹(t)")
ax2.set(xlabel="time", ylabel="probabilities")
ax2.legend()

plt.show()