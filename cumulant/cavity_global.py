import numpy as np
import qutip as qp

ATOM_COUNT = 2
PHOTON_CAPACITY = 2
SPIN_NUM = ATOM_COUNT / 2.0
SPIN_DIM = ATOM_COUNT + 1 # 2 * N / 2 + 1
FOCK_DIM = PHOTON_CAPACITY + 1

hbar = 1
delta = 0.0
true_g = 10.0
g = true_g / np.sqrt(ATOM_COUNT)
kappa = 40.0
gamma = 9.0
nu = 1.0

op_sz = qp.tensor(qp.spin_Jz(SPIN_NUM), qp.qeye(FOCK_DIM))
op_sp = qp.tensor(qp.spin_Jp(SPIN_NUM), qp.qeye(FOCK_DIM))
op_sm = qp.tensor(qp.spin_Jm(SPIN_NUM), qp.qeye(FOCK_DIM))
op_ad = qp.tensor(qp.qeye(SPIN_DIM), qp.create(FOCK_DIM))
op_a = qp.tensor(qp.qeye(SPIN_DIM), qp.destroy(FOCK_DIM))
op_n = op_ad * op_a
op_ac = op_a # Yeah... same operator but different name... it is just to be sure with what we work, when we call op_ac or op_a for example
op_adc = op_ad

init_state = qp.tensor(qp.spin_state(ATOM_COUNT/2, ATOM_COUNT/2), qp.fock(PHOTON_CAPACITY + 1, 0)) # '+1' because we count the 'zero' ladder
#init_state = qp.tensor(qp.spin_state(ATOM_COUNT/2, -ATOM_COUNT/2), qp.fock(PHOTON_CAPACITY + 1, 0)) # Superradiance from julia tutorial
#init_state = qp.tensor(qp.spin_state(ATOM_COUNT/2, -ATOM_COUNT/2), qp.coherent(PHOTON_CAPACITY + 1, 1))

def mean_value(state, op):
    if state.type == 'ket':
        return (state * state.dag() * op).tr()
    elif state.type == 'bra':
        return (state.dag() * state * op).tr()
    return (state * op).tr()

#mean_ac = mean_value(init_state, op_ac)
#mean_adc = mean_value(init_state, op_adc)