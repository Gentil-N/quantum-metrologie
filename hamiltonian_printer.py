import qutip as qx
from matplotlib import pyplot as plt
import numpy as np
import math

FOCK_DIM = 2
ATOM_NUMBER = 1
REVERSE = False

def spin_dim(atom_number):
    return atom_number + 1 # 2 * N/2 + 1

def tens_from_spin(spin_op):
    if REVERSE:
        return qx.tensor(qx.qeye(FOCK_DIM), spin_op)
    else:
        return qx.tensor(spin_op, qx.qeye(FOCK_DIM))

def tens_from_fock(fock_op):
    if REVERSE:
        return qx.tensor(fock_op, qx.qeye(spin_dim(ATOM_NUMBER)))
    else:
        return qx.tensor(qx.qeye(spin_dim(ATOM_NUMBER)), fock_op)

ad = qx.create(FOCK_DIM)
a = qx.destroy(FOCK_DIM)
#print(ad)
#print(a)

n = ad * a

n_all = tens_from_fock(n)
ad_all = tens_from_fock(ad)
a_all = tens_from_fock(a)

sz = qx.spin_Jz(ATOM_NUMBER / 2)
sp = qx.spin_Jp(ATOM_NUMBER / 2)
sm = qx.spin_Jm(ATOM_NUMBER / 2)

#print(sz)
#print(sp)
#print(sm)

sz_all = tens_from_spin(sz)
sp_all = tens_from_spin(sp)
sm_all = tens_from_spin(sm)

hamiltonian = n_all + sz_all + a_all * sp_all + ad_all * sm_all

#hun = qx.Qobj([[1], [1], [1], [1], [1]])
#op_hun = hun*hun.dag()
#
#print(op_hun)
#
#op_hun_tot = qx.tensor(op_hun, op_hun)
#print(op_hun_tot)
#print(hamiltonian)
#res = hamiltonian * op_hun_tot - op_hun_tot * hamiltonian
#
#h_full = res.full()

h_full = hamiltonian.full()

for row in h_full:
    for num in row:
        if num != 0:
            print(" X ", end ='')
        else:
            print(" - ", end='')
    print("")