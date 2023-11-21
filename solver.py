import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot
from matplotlib import pyplot as plt
import qutip as qp

ATOM_COUNT = 1
PHOTON_CAPACITY = 1
SPIN_NUM = ATOM_COUNT / 2.0
SPIN_DIM = ATOM_COUNT + 1 # 2 * N / 2 + 1
FOCK_DIM = PHOTON_CAPACITY + 1

HBAR = 1
delta = 0.0
G = 1.0 / np.sqrt(ATOM_COUNT)
gamma = 1.0
kappa = 1.0
nu = 1.0

op_jz = qp.tensor(qp.spin_Jz(SPIN_NUM), qp.qeye(FOCK_DIM))
op_jp = qp.tensor(qp.spin_Jp(SPIN_NUM), qp.qeye(FOCK_DIM))
op_jm = qp.tensor(qp.spin_Jm(SPIN_NUM), qp.qeye(FOCK_DIM))
op_ad = qp.tensor(qp.qeye(SPIN_DIM), qp.create(FOCK_DIM))
op_a = qp.tensor(qp.qeye(SPIN_DIM), qp.destroy(FOCK_DIM))
op_n = op_ad * op_a

def mean_value(init_state, op):
        if init_state.type == 'ket':
            return (init_state * init_state.dag() * op).tr()
        elif init_state.type == 'bra':
            return (init_state.dag() * init_state * op).tr()
        return (init_state * op).tr()

def correlation_system(t, x):

    return [ (-4j)*G*HBAR*x[6]*x[12] + (-4j)*G*HBAR*np.conj(x[4])*np.conj(x[3]) + (-4j)*G*HBAR*x[10]*x[1] + (8j)*G*HBAR*x[10]*np.conj(x[3])*x[12] + (2j)*G*HBAR*HBAR*x[1] + (-2.0)*gamma*HBAR*x[6]*np.conj(x[3]) + (-2.0)*gamma*HBAR*x[6]*np.conj(x[3]) + (-2.0)*gamma*HBAR*x[10]*x[0] + (4.0)*gamma*HBAR*x[10]*np.conj(x[3])*np.conj(x[3]) + (1.0)*gamma*HBAR*HBAR*x[0] + (2.0)*nu*HBAR*x[6]*np.conj(x[3]) + (2.0)*nu*HBAR*x[6]*np.conj(x[3]) + (2.0)*nu*HBAR*x[10]*x[0] + (-4.0)*nu*HBAR*x[10]*np.conj(x[3])*np.conj(x[3]) + (-3.0)*nu*HBAR*HBAR*x[0], 
(1j)*delta*x[1] + (1j)*G*x[0] + (-0.5)*kappa*x[1] + (-2j)*G*HBAR*np.conj(x[4])*x[12] + (-2j)*G*HBAR*np.conj(x[4])*x[12] + (-2j)*G*HBAR*x[10]*x[2] + (4j)*G*HBAR*x[10]*x[12]*x[12] + (-1.0)*gamma*HBAR*x[6]*x[12] + (-1.0)*gamma*HBAR*np.conj(x[4])*np.conj(x[3]) + (-1.0)*gamma*HBAR*x[10]*x[1] + (2.0)*gamma*HBAR*x[10]*np.conj(x[3])*x[12] + (1.0)*nu*HBAR*x[6]*x[12] + (1.0)*nu*HBAR*np.conj(x[4])*np.conj(x[3]) + (1.0)*nu*HBAR*x[10]*x[1] + (-2.0)*nu*HBAR*x[10]*np.conj(x[3])*x[12] + (-1.0)*nu*HBAR*HBAR*x[1], 
(2j)*delta*x[2] + (2j)*G*x[1] + (-1.0)*kappa*x[2], 
(2j)*G*HBAR*x[4] + (-1.0)*gamma*HBAR*x[5] + (-1.0)*gamma*HBAR*HBAR*x[3] + (1.0)*nu*HBAR*x[5], 
(-1j)*delta*x[4] + (-1j)*G*x[5] + (-0.5)*kappa*x[4] + (-1j)*G*HBAR*x[8]*np.conj(x[12]) + (-1j)*G*HBAR*x[8]*np.conj(x[12]) + (-1j)*G*HBAR*np.conj(x[3])*np.conj(x[2]) + (2j)*G*HBAR*np.conj(x[3])*np.conj(x[12])*np.conj(x[12]) + (1j)*G*HBAR*np.conj(x[8])*np.conj(x[12]) + (1j)*G*HBAR*np.conj(x[1])*x[12] + (1j)*G*HBAR*x[3]*x[7] + ((-0-2j))*G*HBAR*x[3]*x[12]*np.conj(x[12]) + (1.0)*gamma*HBAR*x[9]*np.conj(x[12]) + (1.0)*gamma*HBAR*np.conj(x[1])*np.conj(x[3]) + (1.0)*gamma*HBAR*x[3]*x[8] + (-2.0)*gamma*HBAR*x[3]*np.conj(x[3])*np.conj(x[12]) + (-1.0)*nu*HBAR*x[9]*np.conj(x[12]) + (-1.0)*nu*HBAR*np.conj(x[1])*np.conj(x[3]) + (-1.0)*nu*HBAR*x[3]*x[8] + (2.0)*nu*HBAR*x[3]*np.conj(x[3])*np.conj(x[12]) + (-2.0)*nu*HBAR*HBAR*x[4], 
(-1j)*G*HBAR*x[9]*np.conj(x[12]) + (-1j)*G*HBAR*np.conj(x[1])*np.conj(x[3]) + (-1j)*G*HBAR*x[3]*x[8] + (2j)*G*HBAR*x[3]*np.conj(x[3])*np.conj(x[12]) + (-2j)*G*HBAR*HBAR*x[4] + (2j)*G*HBAR*x[11]*np.conj(x[12]) + (2j)*G*HBAR*x[4]*x[10] + (2j)*G*HBAR*x[10]*x[4] + ((-0-4j))*G*HBAR*x[10]*x[10]*np.conj(x[12]) + (1j)*G*HBAR*np.conj(x[0])*x[12] + (1j)*G*HBAR*np.conj(x[8])*x[3] + (1j)*G*HBAR*x[3]*np.conj(x[8]) + ((-0-2j))*G*HBAR*x[3]*x[3]*x[12] + (1.0)*gamma*HBAR*np.conj(x[0])*np.conj(x[3]) + (1.0)*gamma*HBAR*x[9]*x[3] + (1.0)*gamma*HBAR*x[3]*x[9] + (-2.0)*gamma*HBAR*x[3]*x[3]*np.conj(x[3]) + (-1.0)*gamma*HBAR*HBAR*x[5] + (-1.0)*gamma*HBAR*x[11]*x[3] + (-1.0)*gamma*HBAR*x[5]*x[10] + (-1.0)*gamma*HBAR*x[10]*x[5] + (2.0)*gamma*HBAR*x[10]*x[10]*x[3] + (-1.0)*nu*HBAR*np.conj(x[0])*np.conj(x[3]) + (-1.0)*nu*HBAR*x[9]*x[3] + (-1.0)*nu*HBAR*x[3]*x[9] + (2.0)*nu*HBAR*x[3]*x[3]*np.conj(x[3]) + (-4.0)*nu*HBAR*HBAR*x[5] + (-2.0)*nu*HBAR*HBAR*HBAR*x[3] + (1.0)*nu*HBAR*x[11]*x[3] + (1.0)*nu*HBAR*x[5]*x[10] + (1.0)*nu*HBAR*x[10]*x[5] + (-2.0)*nu*HBAR*x[10]*x[10]*x[3], 
(-1j)*G*HBAR*x[0]*np.conj(x[12]) + (-1j)*G*HBAR*x[8]*np.conj(x[3]) + (-1j)*G*HBAR*np.conj(x[3])*x[8] + (2j)*G*HBAR*np.conj(x[3])*np.conj(x[3])*np.conj(x[12]) + (1j)*G*HBAR*x[9]*x[12] + (1j)*G*HBAR*np.conj(x[8])*np.conj(x[3]) + (1j)*G*HBAR*x[3]*x[1] + ((-0-2j))*G*HBAR*x[3]*np.conj(x[3])*x[12] + (-2j)*G*HBAR*x[11]*x[12] + (-2j)*G*HBAR*np.conj(x[4])*x[10] + (-2j)*G*HBAR*x[10]*np.conj(x[4]) + (4j)*G*HBAR*x[10]*x[10]*x[12] + (1.0)*gamma*HBAR*x[9]*np.conj(x[3]) + (1.0)*gamma*HBAR*x[9]*np.conj(x[3]) + (1.0)*gamma*HBAR*x[3]*x[0] + (-2.0)*gamma*HBAR*x[3]*np.conj(x[3])*np.conj(x[3]) + (-1.0)*gamma*HBAR*x[11]*np.conj(x[3]) + (-1.0)*gamma*HBAR*x[6]*x[10] + (-1.0)*gamma*HBAR*x[10]*x[6] + (2.0)*gamma*HBAR*x[10]*x[10]*np.conj(x[3]) + (-1.0)*nu*HBAR*x[9]*np.conj(x[3]) + (-1.0)*nu*HBAR*x[9]*np.conj(x[3]) + (-1.0)*nu*HBAR*x[3]*x[0] + (2.0)*nu*HBAR*x[3]*np.conj(x[3])*np.conj(x[3]) + (-5.0)*nu*HBAR*HBAR*x[6] + (2.0)*nu*HBAR*HBAR*HBAR*np.conj(x[3]) + (1.0)*nu*HBAR*x[11]*np.conj(x[3]) + (1.0)*nu*HBAR*x[6]*x[10] + (1.0)*nu*HBAR*x[10]*x[6] + (-2.0)*nu*HBAR*x[10]*x[10]*np.conj(x[3]), 
(1j)*G*x[8] + (-1j)*G*np.conj(x[8]) + (-1.0)*kappa*x[7], 
(-1j)*delta*x[8] + (-1j)*G*x[9] + (-0.5)*kappa*x[8] + (-2j)*G*HBAR*np.conj(x[4])*np.conj(x[12]) + (-2j)*G*HBAR*x[4]*x[12] + (-2j)*G*HBAR*x[10]*x[7] + (4j)*G*HBAR*x[10]*x[12]*np.conj(x[12]) + (-1.0)*gamma*HBAR*x[6]*np.conj(x[12]) + (-1.0)*gamma*HBAR*x[4]*np.conj(x[3]) + (-1.0)*gamma*HBAR*x[10]*x[8] + (2.0)*gamma*HBAR*x[10]*np.conj(x[3])*np.conj(x[12]) + (1.0)*nu*HBAR*x[6]*np.conj(x[12]) + (1.0)*nu*HBAR*x[4]*np.conj(x[3]) + (1.0)*nu*HBAR*x[10]*x[8] + (-2.0)*nu*HBAR*x[10]*np.conj(x[3])*np.conj(x[12]) + (-1.0)*nu*HBAR*HBAR*x[8] + (-2j)*G*HBAR*x[10], 
(2j)*G*HBAR*x[6]*np.conj(x[12]) + (2j)*G*HBAR*x[4]*np.conj(x[3]) + (2j)*G*HBAR*x[10]*x[8] + ((-0-4j))*G*HBAR*x[10]*np.conj(x[3])*np.conj(x[12]) + (-2j)*G*HBAR*x[5]*x[12] + (-2j)*G*HBAR*np.conj(x[4])*x[3] + (-2j)*G*HBAR*x[10]*np.conj(x[8]) + (4j)*G*HBAR*x[10]*x[3]*x[12] + (-2j)*G*HBAR*HBAR*np.conj(x[8]) + (-2.0)*gamma*HBAR*x[5]*np.conj(x[3]) + (-2.0)*gamma*HBAR*x[6]*x[3] + (-2.0)*gamma*HBAR*x[10]*x[9] + (4.0)*gamma*HBAR*x[10]*x[3]*np.conj(x[3]) + (-2.0)*gamma*HBAR*HBAR*x[9] + (2.0)*nu*HBAR*x[5]*np.conj(x[3]) + (2.0)*nu*HBAR*x[6]*x[3] + (2.0)*nu*HBAR*x[10]*x[9] + (-4.0)*nu*HBAR*x[10]*x[3]*np.conj(x[3]) + (4.0)*nu*HBAR*HBAR*x[11], 
((-0-1j))*G*HBAR*x[8] + (1j)*G*HBAR*np.conj(x[8]) + (1.0)*gamma*HBAR*x[9] + (-1.0)*nu*HBAR*x[9] + (-2.0)*nu*HBAR*HBAR*x[10], 
(-2j)*G*HBAR*x[6]*np.conj(x[12]) + (-2j)*G*HBAR*x[4]*np.conj(x[3]) + (-2j)*G*HBAR*x[10]*x[8] + (4j)*G*HBAR*x[10]*np.conj(x[3])*np.conj(x[12]) + (1j)*G*HBAR*HBAR*x[8] + (2j)*G*HBAR*x[5]*x[12] + (2j)*G*HBAR*np.conj(x[4])*x[3] + (2j)*G*HBAR*x[10]*np.conj(x[8]) + ((-0-4j))*G*HBAR*x[10]*x[3]*x[12] + (1j)*G*HBAR*HBAR*np.conj(x[8]) + (2.0)*gamma*HBAR*x[5]*np.conj(x[3]) + (2.0)*gamma*HBAR*x[6]*x[3] + (2.0)*gamma*HBAR*x[10]*x[9] + (-4.0)*gamma*HBAR*x[10]*x[3]*np.conj(x[3]) + (1.0)*gamma*HBAR*HBAR*x[9] + (-2.0)*nu*HBAR*x[5]*np.conj(x[3]) + (-2.0)*nu*HBAR*x[6]*x[3] + (-2.0)*nu*HBAR*x[10]*x[9] + (4.0)*nu*HBAR*x[10]*x[3]*np.conj(x[3]) + (1.0)*nu*HBAR*HBAR*x[9] + (2.0)*nu*HBAR*HBAR*HBAR*x[10] + (-4.0)*nu*HBAR*HBAR*x[11], 
(1j)*delta*x[12] + (1j)*G*np.conj(x[3]) + (-0.5)*kappa*x[12] ]


init_state = qp.tensor(qp.spin_state(ATOM_COUNT/2, ATOM_COUNT/2), qp.fock(PHOTON_CAPACITY + 1, 0))
x_init = [
    (1.0 + 0.0j) * mean_value(init_state, op_jp * op_jp),
    (1.0 + 0.0j) * mean_value(init_state, op_jp * op_ad),
    (1.0 + 0.0j) * mean_value(init_state, op_ad * op_ad),
    (1.0 + 0.0j) * mean_value(init_state, op_jm),
    (1.0 + 0.0j) * mean_value(init_state, op_jz * op_a),
    (1.0 + 0.0j) * mean_value(init_state, op_jz * op_jm),
    (1.0 + 0.0j) * mean_value(init_state, op_jz * op_jp),
    (1.0 + 0.0j) * mean_value(init_state, op_ad * op_a),
    (1.0 + 0.0j) * mean_value(init_state, op_jp * op_a),
    (1.0 + 0.0j) * mean_value(init_state, op_jm * op_jp),
    (1.0 + 0.0j) * mean_value(init_state, op_jz),
    (1.0 + 0.0j) * mean_value(init_state, op_jz * op_jz),
    (1.0 + 0.0j) * mean_value(init_state, op_ad)]

print(x_init)
#exit()

print("Solver started...", end='')
solution = solve_ivp(correlation_system, [0.0, 10.0], x_init, first_step=0.1, max_step=0.1)
print("Done")

fig, ax = pyplot.subplots(1, 1)
#for i in range(len(x_init)):
#    ax.plot(solution.t, np.real(solution.y[i]), label=str(i))
ax.plot(solution.t, np.real(np.conj(solution.y[3])), label="<S+>")
ax.plot(solution.t, np.real(solution.y[7]), label="<n>")
ax.plot(solution.t, np.real(solution.y[10]), label="<Sz>")
ax.legend()
ax.set_xlabel('time')
ax.set_ylabel('correlations')

plt.show()