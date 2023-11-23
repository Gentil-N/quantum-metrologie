import numpy as np
from correlation_global import *

def correlation_system_ad_order_1(t, x):
	return [ (2j)*g*hbar*x[2]*np.conj(x[1]) + (-1.0)*gamma*hbar*x[2]*x[0] + (-1.0)*gamma*hbar*hbar*x[0] + (1.0)*nu*hbar*x[2]*x[0], 
(1j)*delta*x[1] + (1j)*g*np.conj(x[0]) + (-0.5)*kappa*x[1], 
(-1j)*g*hbar*np.conj(x[0])*np.conj(x[1]) + (1j)*g*hbar*x[0]*x[1] + (1.0)*gamma*hbar*x[0]*np.conj(x[0]) + (-1.0)*nu*hbar*x[0]*np.conj(x[0]) + (-2.0)*nu*hbar*hbar*x[2] ]

def get_init_vec_ad_order_1(init_state):
	return [
(1.0 + 0.0j) * mean_value(init_state, op_sm),
(1.0 + 0.0j) * mean_value(init_state, op_ad),
(1.0 + 0.0j) * mean_value(init_state, op_sz),
]