import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot
from matplotlib import pyplot as plt

delta = 0.0
g = 1.5
gamma = 0.25
kappa = 1.0
nu = 4.0

def correlation_system(t, x_vec):
    return [-1j * x_vec[2] + 1j * g * np.conj(x_vec[2]) - kappa * x_vec[0],
            nu + 1j * g * x_vec[2] - 1j * g * np.conj(x_vec[2]) - gamma * x_vec[1] - nu * x_vec[1],
            1j * g * x_vec[1] - 1j * g * x_vec[0] + 1j * delta * x_vec[2] - gamma / 2 * x_vec[2] - kappa / 2 * x_vec[2] - nu / 2 * x_vec[2] + 2j * g * x_vec[1] * x_vec[0]]

x_vec_init = [0.0 + 0.0j, 0.0 + 0.0j, 0.0 + 0.0j]

solution = solve_ivp(correlation_system, [0.0, 10.0], x_vec_init, first_step=0.1, max_step=0.1)

fig, ax = pyplot.subplots(1, 1)
ax.plot(solution.t, np.real(solution.y[0]), label='photon number')
ax.plot(solution.t, np.real(solution.y[1]), label='excited state')
ax.plot(solution.t, np.real(solution.y[2]), label='???')
ax.legend()
ax.set_xlabel('time')
ax.set_ylabel('correlations')

plt.show()