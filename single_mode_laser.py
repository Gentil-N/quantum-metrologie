import qutip as qx
from matplotlib import pyplot as plt
import numpy as np
import math
import cmath
import scipy.constants as cs

E0 = 0.0
HBAR = 0.01
OMEGA_0 = 1.0
OMEGA = 0.5
T0 = 0.0
G = 0.8
N_ATOM_STATE = 2
N_PHOTON = 2

A_IDENTITY = qx.qeye(N_PHOTON)
B_IDENTITY = qx.qeye(N_ATOM_STATE)

op_a = qx.destroy(N_PHOTON)
op_ad = qx.create(N_PHOTON)
op_a_tens = qx.tensor(op_a, B_IDENTITY)
op_ad_tens = qx.tensor(op_ad, B_IDENTITY)
op_Na = op_ad * op_a
op_Na_tens = qx.tensor(op_Na, B_IDENTITY)

op_b = qx.destroy(N_ATOM_STATE)
op_bd = qx.create(N_ATOM_STATE)
op_b_tens = qx.tensor(A_IDENTITY, op_b)
op_bd_tens = qx.tensor(A_IDENTITY, op_bd)
op_Nb = op_bd * op_b
op_Nb_tens = qx.tensor(A_IDENTITY, op_Nb)

op_h_system = qx.tensor(A_IDENTITY, E0 + HBAR * OMEGA_0 / 2 * (op_bd * op_b - op_b * op_bd))

op_h_bath = qx.tensor(HBAR * OMEGA * (op_ad * op_a + 1/2), B_IDENTITY)

op_h_sb_independent = HBAR * G * (op_a_tens * op_bd_tens * np.exp(1j * (OMEGA_0 - OMEGA)) + op_a_tens * op_b_tens * np.exp(-1j * (OMEGA_0 + OMEGA)))

def h_sb_dependent(t, args):
    return np.exp(t-T0)

h_total = [op_h_system + op_h_bath, [op_h_sb_independent, h_sb_dependent]]

psi_init = qx.tensor(A_IDENTITY, qx.fock_dm(N_ATOM_STATE, 1))
time_range = np.linspace(T0, T0 + 10, 3)

result = qx.mesolve(h_total, psi_init, time_range)

print(result.states[0], " ", result.states[0].dag() * op_Na_tens * result.states[0])
print(" ")
print(result.states[1], " ", result.states[1].dag() * op_Na_tens * result.states[1])

fig = plt.figure(num=0)
ax = fig.subplots(nrows=1, ncols=1)
ax.set_title(label="Photons")
ax.plot(time_range, qx.expect(op_Na_tens, result.states))
ax.set(xlabel="time", ylabel="mean N")
plt.show()