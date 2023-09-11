import qutip as qx
from matplotlib import pyplot as plt
import numpy as np
import math

OMEGA_PHOTONS = 1.0
OMEGA_ATOM = 1.0
HBAR = 1.0
HILBERT_DIM_PHOTON = 3
G = 10.0
KAPPA = 3.0
REPUMP = 1.0
ENABLE_DECAY = True
ENABLE_REPUMP = True

def get_u_from_hamiltonian(hamitlonian, time):
    return (-1j * hamitlonian * time / HBAR).expm()

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

op_hamiltonian_base = HBAR * OMEGA_ATOM * tens_from_spin(op_sz) + HBAR * OMEGA_PHOTONS * tens_from_fock(op_ad * op_a)
op_hamiltonian_interaction = HBAR * G * (qx.tensor(op_s_pos, op_a) + qx.tensor(op_s_neg, op_ad))
op_hamiltonian =op_hamiltonian_base + op_hamiltonian_interaction

psi_init = qx.tensor(spin_up, qx.fock(HILBERT_DIM_PHOTON, 0)) # = |e,0>
print(psi_init.type == 'ket')
psi_init = psi_init * psi_init.dag()
print(psi_init.type == 'oper')

time_range = np.linspace(0, 5, 1000)

result = qx.mesolve(op_hamiltonian, psi_init, time_range, op_collapsing, [])

#print(create_coherent_state(HILBERT_DIM_PHOTON, COHERENT_ALPHA))
#print(normalize_state(create_coherent_state(2, 1)))

ket_e0 = qx.tensor(spin_up, qx.fock(HILBERT_DIM_PHOTON, 0)) # |e,0>
ket_g0 = qx.tensor(spin_down, qx.fock(HILBERT_DIM_PHOTON, 0)) # |g,0>
ket_e1 = qx.tensor(spin_up, qx.fock(HILBERT_DIM_PHOTON, 1)) # |e,1>
ket_g1 = qx.tensor(spin_down, qx.fock(HILBERT_DIM_PHOTON, 1)) # |g,1>
ket_e2 = qx.tensor(spin_up, qx.fock(HILBERT_DIM_PHOTON, 2)) # |e,2>
ket_g2 = qx.tensor(spin_down, qx.fock(HILBERT_DIM_PHOTON, 2)) # |g,2>

prob_e0 = qx.expect(ket_e0 * ket_e0.dag(), result.states)
prob_g0 = qx.expect(ket_g0 * ket_g0.dag(), result.states)
prob_e1 = qx.expect(ket_e1 * ket_e1.dag(), result.states)
prob_g1 = qx.expect(ket_g1 * ket_g1.dag(), result.states)
prob_e2 = qx.expect(ket_e2 * ket_e2.dag(), result.states)
prob_g2 = qx.expect(ket_g2 * ket_g2.dag(), result.states)

op_n_all = tens_from_fock(op_n)
res_n = qx.expect(op_n_all, result.states)
res_n_mesolve = qx.mesolve(op_hamiltonian, psi_init, time_range, op_collapsing, [op_n_all]).expect[0]
res_n_mesolve = np.array(res_n_mesolve)

#res_a_test = []
#i = 0
#for t in time_range:
#    res_a_test.append((result.states[i].dag() * tens_from_fock(op_ad) * np.exp(-1j * OMEGA_PHOTONS * t) * result.states[i]).tr())
#    i += 1

#G1 = qx.correlation_2op_2t(op_hamiltonian, psi_init, None, time_range, op_collapsing, tens_from_fock(op_ad), tens_from_fock(op_a))
#print("**********************", G1)
#if len(G1.shape) == 2:
#    G1 = G1[0]
#print("**********************", G1)
#n = np.array(res_n_mesolve)
#g1 = G1 / np.sqrt(n[0] * n) 
#print("**********************", g1)

g1 = []
i = 0
for t in time_range:

    if t == 0.0:
        g1.append(0.0)
        i += 1
        continue

    op_u = get_u_from_hamiltonian(op_hamiltonian, t)
    op_ad_evol_heis = op_u.dag() * tens_from_fock(op_ad) * op_u
    op_a_evol_heis = op_u.dag() * tens_from_fock(op_a) * op_u
    norm_factor = np.sqrt((result.states[i] * op_ad_evol_heis * op_a_evol_heis).tr() * (result.states[i] * tens_from_fock(op_ad) * tens_from_fock(op_a)).tr())
    g1.append((result.states[i] * op_ad_evol_heis * tens_from_fock(op_a)).tr() / norm_factor)
    i += 1

g2 = []
i = 0
for t in time_range:

    if t == 0.0:
        g2.append(0.0)
        i += 1
        continue

    op_u = get_u_from_hamiltonian(op_hamiltonian, t)
    op_ad_evol_heis = op_u.dag() * tens_from_fock(op_ad) * op_u
    op_a_evol_heis = op_u.dag() * tens_from_fock(op_a) * op_u
    norm_factor = (result.states[i] * tens_from_fock(op_ad) * tens_from_fock(op_a)).tr()**2
    g2.append((result.states[i] * tens_from_fock(op_ad) * op_ad_evol_heis * op_a_evol_heis * tens_from_fock(op_a)).tr() / norm_factor)
    i += 1

###############
### PLOTING ###
###############

common_title = "Simple Cavity (1 atom - 1 excitation) start |e, 0>"

fig0 = plt.figure(num=0)
ax0 = fig0.subplots(nrows=1, ncols=1)
ax0.set_title(label=common_title)
ax0.plot(time_range, prob_e0, label= "|<e,0|ψ(t)>|²")
ax0.plot(time_range, prob_g0, label= "|<g,0|ψ(t)>|²")
ax0.plot(time_range, prob_e1, label= "|<e,1|ψ(t)>|²")
ax0.plot(time_range, prob_g1, label= "|<g,1|ψ(t)>|²")
ax0.plot(time_range, prob_e2, label= "|<e,2|ψ(t)>|²")
ax0.plot(time_range, prob_g2, label= "|<g,2|ψ(t)>|²")
ax0.set(xlabel="time", ylabel="probabilities")
ax0.legend()

fig1 = plt.figure(num=1)
ax1 = fig1.subplots(nrows=1, ncols=1)
ax1.set_title(label=common_title)
ax1.plot(time_range, res_n, label="|<ψ(t)|N|ψ(t)>|²")
ax1.set(xlabel="time", ylabel="photon number")
ax1.legend()

fig2 = plt.figure(num=2)
ax2 = fig2.subplots(nrows=1, ncols=1)
ax2.set_title(label=common_title)
ax2.plot(time_range, np.real(g1), label="g¹(t) real part")
ax2.plot(time_range, np.imag(g1), label="g¹(t) imaginary part")
ax2.set(xlabel="time", ylabel="probabilities")
ax2.legend()

fig3 = plt.figure(num=3)
ax3 = fig3.subplots(nrows=1, ncols=1)
ax3.set_title(label=common_title)
ax3.plot(time_range, prob_e0 + prob_g0 + prob_e1 + prob_g1 + prob_e2 + prob_g2, label="| (<e,0| + <g,0| + <e,1| + <g,1| + <e,2| + <g,2|) * |ψ(t)> |²")
ax3.set(xlabel="time", ylabel="probabilities")
ax3.legend()

fig4 = plt.figure(num=4)
ax4 = fig4.subplots(nrows=1, ncols=1)
ax4.set_title(label=common_title)
ax4.plot(time_range, np.real(g2), label="g²(t) real part")
ax4.plot(time_range, np.imag(g2), label="g²(t) imaginary part")
ax4.set(xlabel="time", ylabel="probabilities")
ax4.legend()

plt.show()