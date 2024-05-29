import numpy as np
import qutip as qp
import scipy.constants as const
import math
from enum import Enum


#
# ATOM_COUNT = 10_000
# gamma = kappa / 103.95373
# <sz> steady = 0.24
# <sz> max = 5_000
# ratio = 0.24 / 5000 ~ 0.000048
#
# ATOM_COUNT = 100_000
# gamma = kappa / 103.95373
# <sz> steady = 412
# <sz> max = 50_000
# ratio = 412 / 50_000 ~ 0.00824
#
# ATOM_COUNT = 100_000
# gamma = kappa / 103.9537874
# <sz> steady = 2.01
# <sz> max = 50_000
# ratio = 2 / 50_000 = 0.00004
#
# ATOM_COUNT = 500_000
# gamma = kappa / 103.9537874
# <sz> steady = 70_000
# <sz> max = 250_000
# ratio = 70_000 / 250_000 = 0.28
# BUUUUUUUUG!!!
#
# ATOM_COUNT = 300_000
# gamma = kappa / 103.954647
# <sz> steady = 73
# <sz> max = 150_000
# ratio = 73 / 150_000 ~ 0.000487
#
# ATOM_COUNT = 300_000
# gamma = kappa / 103.954651
# <sz> steady = 5.6
# <sz> max = 150_000
# ratio = 5.6 / 150_000 ~ 0.000037333
#
#

ATOM_COUNT = 10000
PHOTON_CAPACITY = 10
SPIN_NUM = ATOM_COUNT / 2.0
SPIN_DIM = ATOM_COUNT + 1 # 2 * N / 2 + 1
FOCK_DIM = PHOTON_CAPACITY + 1

hbar = 1.0 #1.05457182e-34 # hbar canceled inside the matrix definition in qtip (i.e S-matix doesn't implement hbar)
delta = 0.0
true_g = 0.0 #10.0
g = 5117.82 #true_g / np.sqrt(ATOM_COUNT)
kappa = 2 * np.pi * 780*const.kilo #40.0
gamma = kappa / 103.95373 #9.0
nu = 2* np.pi * 7.5*const.kilo

op_sz = qp.tensor(qp.spin_Jz(SPIN_NUM), qp.qeye(FOCK_DIM))
op_sp = qp.tensor(qp.spin_Jp(SPIN_NUM), qp.qeye(FOCK_DIM))
op_sm = qp.tensor(qp.spin_Jm(SPIN_NUM), qp.qeye(FOCK_DIM))
op_ad = qp.tensor(qp.qeye(SPIN_DIM), qp.create(FOCK_DIM))
op_a = qp.tensor(qp.qeye(SPIN_DIM), qp.destroy(FOCK_DIM))
op_n = op_ad * op_a
op_ac = op_a # Yeah... same operator but different name... it is just to be sure with what we work, when we call op_ac or op_a for example
op_adc = op_ad

def mean_value(state, op):
    if state.type == 'ket':
        return (state * state.dag() * op).tr()
    elif state.type == 'bra':
        return (state.dag() * state * op).tr()
    return (state * op).tr()

class QOpInit(Enum):
    op_a = 1
    op_ad = 2
    op_sp = 3
    op_sm = 4
    op_sz = 5

class State:
    j = 0
    m = 0
    n = 0
    def __init__(self, j, m, n):
        self.j = j
        self.m = m
        self.n = n
    def cop(self):
        return State(self.j, self.m, self.n)

def proj_op_state(qop_init, state):
    if qop_init == QOpInit.op_a:
        if state.n < 0:
            return 0
        factor = math.sqrt(state.n)
        state.n -= 1
        return factor
    if qop_init == QOpInit.op_ad:
        state.n += 1
        return math.sqrt(state.n)
    if qop_init == QOpInit.op_sp:
        if state.m > state.j:
            return 0
        factor = hbar * math.sqrt(state.j * (state.j + 1) - state.m * (state.m + 1))
        state.m += 1
        return factor
    if qop_init == QOpInit.op_sm:
        if state.m < -state.j:
            return 0
        factor = hbar * math.sqrt(state.j * (state.j + 1) - state.m * (state.m - 1))
        state.m -= 1
        return factor
    if qop_init == QOpInit.op_sz:
        return hbar * state.m

def get_projection(state: State, list_qop_init):
    state_right = state.cop()
    state_left = state.cop()
    factor = 1.0
    list_qop_init.reverse()
    for qop_init in list_qop_init:
        temp = proj_op_state(qop_init, state_right)
        factor *= temp
    if state_left.n == state_right.n and state_left.j == state_right.j and state_left.m == state_right.m:
        return factor + 0j
    else:
        return 0+0j

#init_state = qp.tensor(qp.spin_state(ATOM_COUNT/2, ATOM_COUNT/2), qp.fock(PHOTON_CAPACITY + 1, 0)) # '+1' because we count the 'zero' ladder
#init_state = qp.tensor(qp.spin_state(ATOM_COUNT/2, -ATOM_COUNT/2), qp.fock(PHOTON_CAPACITY + 1, 0)) # Superradiance from julia tutorial
#init_state = qp.tensor(qp.spin_state(ATOM_COUNT/2, -ATOM_COUNT/2), qp.coherent(PHOTON_CAPACITY + 1, 1))

state = State(ATOM_COUNT/2, 0, 0)

#mean_ac = mean_value(init_state, op_ac)
#mean_adc = mean_value(init_state, op_adc)
