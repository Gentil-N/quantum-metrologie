import qutip as qx
from matplotlib import pyplot as plt
import numpy as np
import math

OMEGA_PHOTONS = 1.0
OMEGA_ATOM = 1.0
HBAR = 1.0
HILBERT_DIM_PHOTON = 40
G = 10.0
KAPPA = 0.01 *G
REPUMP = 1.0
COHERENT_ALPHA = math.sqrt(20)
ENABLE_DECAY = False
ENABLE_REPUMP = False

def create_coherent_state(fock_space_dim, alpha):
    return np.exp(-np.abs(alpha)**2) *  sum((alpha**n)/np.sqrt(math.factorial(n))*qx.fock(fock_space_dim, n) for n in range(0, fock_space_dim))

def normalize_state(state):
    state *= 1/np.sqrt((state.dag() * state).tr())
    return state

def tens_from_fock(f_op):
    return qx.tensor(qx.qeye(2), f_op)

def tens_from_spin(s_op):
    return qx.tensor(s_op, qx.qeye(HILBERT_DIM_PHOTON))

op_s_pos = 1/2 * (qx.sigmax() + 1j * qx.sigmay())
op_s_neg = 1/2 * (qx.sigmax() - 1j * qx.sigmay())
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

op_hamiltonian = HBAR * OMEGA_ATOM * tens_from_spin(op_sz) + HBAR * OMEGA_ATOM * tens_from_fock(op_ad * op_a) + HBAR * G * (qx.tensor(op_s_pos, op_a) + qx.tensor(op_s_neg, op_ad))

psi_init = qx.tensor(spin_up, qx.fock(HILBERT_DIM_PHOTON, 0)) # = |e,0>

psi_init = qx.tensor(spin_up, qx.coherent(HILBERT_DIM_PHOTON, COHERENT_ALPHA))
#psi_init = qx.tensor(spin_down, create_coherent_state(HILBERT_DIM_PHOTON, COHERENT_ALPHA))

time_range = np.linspace(0, 10, 1000)

if ENABLE_DECAY:
    result = qx.mesolve(op_hamiltonian, psi_init, time_range, op_collapsing, [])
else:
    result = qx.mesolve(op_hamiltonian, psi_init, time_range, [], [])


#print(create_coherent_state(HILBERT_DIM_PHOTON, COHERENT_ALPHA))
#print(normalize_state(create_coherent_state(2, 1)))

print(result.states[0])
first_mesure = qx.tensor(spin_down, qx.qeye(HILBERT_DIM_PHOTON)) # |g>
second_mesure = qx.tensor(spin_up, qx.qeye(HILBERT_DIM_PHOTON)) # |e>
#third_mesure = psi_init

res_first = qx.expect(first_mesure * first_mesure.dag(), result.states)
res_second = qx.expect(second_mesure * second_mesure.dag(), result.states)
#res_third = qx.expect(third_mesure * third_mesure.dag(), result.states)

#op_n_all = tens_from_fock(op_n)
op_n_all = qx.tensor(spin_up * spin_up.dag(), op_n)
res_n = qx.expect(op_n_all, result.states)

res_pe_test = []
for t in time_range:
    value = 0.0
    for n in range(HILBERT_DIM_PHOTON):
        value += np.exp(-abs(COHERENT_ALPHA)**2) / math.factorial(n) * abs(COHERENT_ALPHA)**(2*n) * math.sin(math.sqrt(n) * G * t)**2
    res_pe_test.append(value)

fig0 = plt.figure(num=0)
ax0 = fig0.subplots(nrows=1, ncols=1)
ax0.set_title(label="Cavity start |α=" + str(COHERENT_ALPHA) + ">")
#ax0.plot(G * time_range, res_first, label = "|<g|ψ(t)>|²")
ax0.plot(time_range, res_second, label= "|<e|ψ(t)>|²")
ax0.set(xlabel="time", ylabel="probabilities")
ax0.legend()

fig1 = plt.figure(num=1)
ax1 = fig1.subplots(nrows=1, ncols=1)
ax1.set_title(label="Cavity start |α=" + str(COHERENT_ALPHA) + ">")
ax1.plot(time_range, res_n, label="|<ψ(t)|N|ψ(t)>|²")
ax1.set(xlabel="time", ylabel="probabilities")
ax1.legend()

fig2 = plt.figure(num=2)
ax2 = fig2.subplots(nrows=1, ncols=1)
ax2.set_title(label="Cavity start |α=" + str(COHERENT_ALPHA) + "> exact result")
ax2.plot(time_range, res_pe_test, label="|<e|ψ(t)>|²")
ax2.set(xlabel="time", ylabel="probabilities")
ax2.legend()

plt.show()
