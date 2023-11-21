import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot
from matplotlib import pyplot as plt

delta = 0.0
g = 1.5
gamma = 0.25
kappa = 1.0
nu = 4.0

def correlation_system(t, x):
    return [-1j * g * x[2] + 1j * g * np.conj(x[2]) - kappa * x[0],
            nu + 1j * g * x[2] - 1j * g * np.conj(x[2]) - gamma * x[1] - nu * x[1],
            1j * g * x[1] - 1j * g * x[0] + 1j * delta * x[2] - gamma / 2 * x[2] - kappa / 2 * x[2] - nu / 2 * x[2] + 2j * g * x[1] * x[0]]

x_init = [0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j]

print("Solver started...", end='')
solution = solve_ivp(correlation_system, [0.0, 10.0], x_init, first_step=0.1, max_step=0.1)
print("Done")

fig, ax = pyplot.subplots(1, 1)
ax.plot(solution.t, np.real(solution.y[0]), label='photon number')
ax.plot(solution.t, np.real(solution.y[1]), label='excited state')
ax.plot(solution.t, np.real(solution.y[2]), label='???')
ax.legend()
ax.set_xlabel('time')
ax.set_ylabel('correlations')

plt.show()