from matplotlib import pyplot as plt
import numpy as np
import qutip as qp
import math
import cavity
from scipy.fft import fft, ifft, fftfreq

def rabi_test():
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
    ax1.plot(time_range, qp.expect(rho_spin_up_fock_0, result.states), label="spin up")
    ax1.plot(time_range, qp.expect(rho_spin_down_fock_1, result.states), label="spin down")
    ax1.plot(time_range, qp.expect(ket_fock_1 * ket_fock_1.dag(), result.states), label="coherent revival")
    ax1.legend()
    ax1.grid()
    plt.show()


N_fock = 10
N_atoms = 6
gamma = 1.0
delta = 2.0
g  = 10.0
mu = 9.0
kappa = 40.0

# Jaynes-Cummings Hamiltonian
a  = qp.tensor(qp.destroy(N_fock), qp.qeye(N_atoms + 1))
sm = qp.tensor(qp.qeye(N_fock), qp.spin_Jm(N_atoms/2.0))
sz = qp.tensor(qp.qeye(N_fock), qp.spin_Jz(N_atoms/2.0))
H = delta * a.dag()*a + g * (a.dag() * sm + a * sm.dag())

# collapse operators
c_ops = [
    np.sqrt(kappa) * a,
    np.sqrt(gamma) * sm,
    np.sqrt(mu) * sm.dag(),
]

# calculate the correlation function using the mesolve solver, and then fft to
# obtain the spectrum. Here we need to make sure to evaluate the correlation
# function for a sufficient long time and sufficiently high sampling rate so
# that the discrete Fourier transform (FFT) captures all the features in the
# resulting spectrum.
tlist = np.linspace(0, 5, 1000)
psi_init = qp.tensor(qp.fock(N_fock, 0), qp.spin_state(N_atoms/2.0, -N_atoms/2.0))
result = qp.mesolve(H, psi_init, tlist, c_ops, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar())
result_rhobar = qp.mesolve(H, a * result.states[-1], tlist, c_ops, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar())

common_title = "Cavity (Dicke Model - Master Equation)"



fig1 = plt.figure(num=1)
ax1 = fig1.subplots(nrows=1, ncols=1)
ax1.set_title(label=common_title)
ax1.plot(tlist, qp.expect(a.dag()*a, result.states), label="<n>")
ax1.plot(tlist, qp.expect(sz, result.states), label="<Sz>")
ax1.set(xlabel="time", ylabel="")
ax1.legend()

corr = qp.expect(a.dag(), result_rhobar.states)
corr /= abs(corr[0])

fig4 = plt.figure(num=4)
ax4 = fig4.subplots(nrows=1, ncols=1)
ax4.set_title(label=common_title)
ax4.plot(tlist, np.real(corr), label="g¹(t) real part")
ax4.plot(tlist, np.imag(corr), label="g¹(t) imaginary part")
ax4.set(xlabel="time", ylabel="")
ax4.legend()
plt.show()



freq_list = fftfreq(len(corr), tlist[1] - tlist[0])
fsignal = fft(corr)

# plot the spectra
fig, ax = plt.subplots(1, 1)
ax.plot(freq_list, fsignal)
ax.legend()
ax.set_xlabel('Frequency')
ax.set_ylabel('Power spectrum')

plt.show()
