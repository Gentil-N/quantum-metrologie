import numpy as np
from matplotlib import pyplot
import qutip
from matplotlib import pyplot as plt

N = 4                      # number of cavity fock states
delta = 1.0
g  = 1.5       # coupling strength
gamma = 0.25               # atom dissipation rate
kappa = 1.0               # cavity dissipation rate
mu = 4

# Jaynes-Cummings Hamiltonian
a  = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
sm = qutip.tensor(qutip.qeye(N), qutip.destroy(2))
H = delta*a.dag()*a + g*(a.dag()*sm + a*sm.dag())

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
tlist = np.linspace(0, 500, 5000)
corr = qutip.correlation_2op_1t(H, None, tlist, c_ops, a.dag(), a)

common_title = "Cavity (Dicke Model - Master Equation)"

#fig4 = plt.figure(num=4)
#ax4 = fig4.subplots(nrows=1, ncols=1)
#ax4.set_title(label=common_title)
#ax4.plot(tlist, np.real(corr), label="g¹(t) real part")
#ax4.plot(tlist, np.imag(corr), label="g¹(t) imaginary part")
#ax4.set(xlabel="time", ylabel="")
#ax4.legend()

wlist1, spec1 = qutip.spectrum_correlation_fft(tlist, corr)

# calculate the power spectrum using spectrum, which internally uses essolve
# to solve for the dynamics (by default)
wlist2 = np.linspace(-3, 3, 200)
spec2 = qutip.spectrum(H, wlist2, c_ops, a.dag(), a)

# plot the spectra
fig, ax = pyplot.subplots(1, 1)
ax.plot(wlist1, spec1, 'b', lw=2, label='eseries method')
ax.plot(wlist2, spec2, 'r--', lw=2, label='me+fft method')
ax.legend()
ax.set_xlabel('Frequency')
ax.set_ylabel('Power spectrum')
ax.set_title('Vacuum Rabi splitting')
ax.set_xlim(wlist2[0], wlist2[-1])

plt.show()