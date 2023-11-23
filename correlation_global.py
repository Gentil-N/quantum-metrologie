import numpy as np
import qutip as qp

ATOM_COUNT = 2
PHOTON_CAPACITY = 2
SPIN_NUM = ATOM_COUNT / 2.0
SPIN_DIM = ATOM_COUNT + 1 # 2 * N / 2 + 1
FOCK_DIM = PHOTON_CAPACITY + 1

hbar = 1
delta = 0.0
g = 10.0 / np.sqrt(ATOM_COUNT)
kappa = 40.0
gamma = 9.0
nu = 1.0

op_sz = qp.tensor(qp.spin_Jz(SPIN_NUM), qp.qeye(FOCK_DIM))
op_sp = qp.tensor(qp.spin_Jp(SPIN_NUM), qp.qeye(FOCK_DIM))
op_sm = qp.tensor(qp.spin_Jm(SPIN_NUM), qp.qeye(FOCK_DIM))
op_ad = qp.tensor(qp.qeye(SPIN_DIM), qp.create(FOCK_DIM))
op_a = qp.tensor(qp.qeye(SPIN_DIM), qp.destroy(FOCK_DIM))
op_n = op_ad * op_a

def mean_value(init_state, op):
    if init_state.type == 'ket':
        return (init_state * init_state.dag() * op).tr()
    elif init_state.type == 'bra':
        return (init_state.dag() * init_state * op).tr()
    return (init_state * op).tr()