import numpy as np
from cavity_global import *

def corr_system_ad_order_1(t, x):
	print(t)
	return [ (1j)*delta*x[0] + (1j)*g*x[1] + (-0.5)*kappa*x[0], 
(-2j)*g*hbar*x[2]*x[0] + (-1.0)*gamma*hbar*x[2]*x[1] + (1.0)*nu*hbar*x[2]*x[1] + (-1.0)*nu*hbar*hbar*x[1], 
(-1j)*g*hbar*x[1]*np.conj(x[0]) + (1j)*g*hbar*np.conj(x[1])*x[0] + (1.0)*gamma*hbar*np.conj(x[1])*x[1] + (-1.0)*nu*hbar*np.conj(x[1])*x[1] + (-2.0)*nu*hbar*hbar*x[2] ]

def corr_get_init_vec_ad_order_1(state):
	return [
(1.0 + 0.0j) * get_projection(state, [QOpInit.op_ad]), # 0
(1.0 + 0.0j) * get_projection(state, [QOpInit.op_sp]), # 1
(1.0 + 0.0j) * get_projection(state, [QOpInit.op_sz]), # 2
]