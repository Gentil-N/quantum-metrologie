import numpy as np
from cavity_global import *

def corr_system_ada_order_2(t, x):
	return [ (1j)*g*x[1] + (-1j)*g*np.conj(x[1]) + (-1.0)*kappa*x[0], 
(-1j)*delta*x[1] + (-1j)*g*x[2] + (-0.5)*kappa*x[1] + (-2j)*g*hbar*x[5]*x[8] + (-2j)*g*hbar*np.conj(x[5])*np.conj(x[8]) + (-2j)*g*hbar*x[3]*x[0] + (4j)*g*hbar*x[3]*np.conj(x[8])*x[8] + (-1.0)*gamma*hbar*x[6]*x[8] + (-1.0)*gamma*hbar*np.conj(x[5])*x[7] + (-1.0)*gamma*hbar*x[3]*x[1] + (2.0)*gamma*hbar*x[3]*x[7]*x[8] + (1.0)*nu*hbar*x[6]*x[8] + (1.0)*nu*hbar*np.conj(x[5])*x[7] + (1.0)*nu*hbar*x[3]*x[1] + (-2.0)*nu*hbar*x[3]*x[7]*x[8] + (-1.0)*nu*hbar*hbar*x[1] + (-2j)*g*hbar*x[3], 
(2j)*g*hbar*x[6]*x[8] + (2j)*g*hbar*np.conj(x[5])*x[7] + (2j)*g*hbar*x[3]*x[1] + ((-0-4j))*g*hbar*x[3]*x[7]*x[8] + (-2j)*g*hbar*x[9]*np.conj(x[8]) + (-2j)*g*hbar*x[5]*np.conj(x[7]) + (-2j)*g*hbar*x[3]*np.conj(x[1]) + (4j)*g*hbar*x[3]*np.conj(x[7])*np.conj(x[8]) + (-2j)*g*hbar*hbar*np.conj(x[1]) + (-2.0)*gamma*hbar*x[9]*x[7] + (-2.0)*gamma*hbar*x[6]*np.conj(x[7]) + (-2.0)*gamma*hbar*x[3]*x[2] + (4.0)*gamma*hbar*x[3]*np.conj(x[7])*x[7] + (-2.0)*gamma*hbar*hbar*x[2] + (2.0)*nu*hbar*x[9]*x[7] + (2.0)*nu*hbar*x[6]*np.conj(x[7]) + (2.0)*nu*hbar*x[3]*x[2] + (-4.0)*nu*hbar*x[3]*np.conj(x[7])*x[7] + (4.0)*nu*hbar*hbar*x[4], 
((-0-1j))*g*hbar*x[1] + (1j)*g*hbar*np.conj(x[1]) + (1.0)*gamma*hbar*x[2] + (-1.0)*nu*hbar*x[2] + (-2.0)*nu*hbar*hbar*x[3], 
(-2j)*g*hbar*x[6]*x[8] + (-2j)*g*hbar*np.conj(x[5])*x[7] + (-2j)*g*hbar*x[3]*x[1] + (4j)*g*hbar*x[3]*x[7]*x[8] + (1j)*g*hbar*hbar*x[1] + (2j)*g*hbar*x[9]*np.conj(x[8]) + (2j)*g*hbar*x[5]*np.conj(x[7]) + (2j)*g*hbar*x[3]*np.conj(x[1]) + ((-0-4j))*g*hbar*x[3]*np.conj(x[7])*np.conj(x[8]) + (1j)*g*hbar*hbar*np.conj(x[1]) + (2.0)*gamma*hbar*x[9]*x[7] + (2.0)*gamma*hbar*x[6]*np.conj(x[7]) + (2.0)*gamma*hbar*x[3]*x[2] + (-4.0)*gamma*hbar*x[3]*np.conj(x[7])*x[7] + (1.0)*gamma*hbar*hbar*x[2] + (-2.0)*nu*hbar*x[9]*x[7] + (-2.0)*nu*hbar*x[6]*np.conj(x[7]) + (-2.0)*nu*hbar*x[3]*x[2] + (4.0)*nu*hbar*x[3]*np.conj(x[7])*x[7] + (1.0)*nu*hbar*hbar*x[2] + (2.0)*nu*hbar*hbar*hbar*x[3] + (-4.0)*nu*hbar*hbar*x[4], 
(1j)*delta*x[5] + (1j)*g*x[6] + (-0.5)*kappa*x[5] + (-1j)*g*hbar*x[10]*x[8] + (-1j)*g*hbar*x[1]*np.conj(x[8]) + (-1j)*g*hbar*x[7]*x[0] + (2j)*g*hbar*x[7]*np.conj(x[8])*x[8] + (1j)*g*hbar*np.conj(x[1])*np.conj(x[8]) + (1j)*g*hbar*np.conj(x[1])*np.conj(x[8]) + (1j)*g*hbar*np.conj(x[7])*x[12] + ((-0-2j))*g*hbar*np.conj(x[7])*np.conj(x[8])*np.conj(x[8]) + (1.0)*gamma*hbar*x[2]*np.conj(x[8]) + (1.0)*gamma*hbar*np.conj(x[1])*x[7] + (1.0)*gamma*hbar*np.conj(x[7])*x[10] + (-2.0)*gamma*hbar*np.conj(x[7])*x[7]*np.conj(x[8]) + (-1.0)*nu*hbar*x[2]*np.conj(x[8]) + (-1.0)*nu*hbar*np.conj(x[1])*x[7] + (-1.0)*nu*hbar*np.conj(x[7])*x[10] + (2.0)*nu*hbar*np.conj(x[7])*x[7]*np.conj(x[8]) + (-2.0)*nu*hbar*hbar*x[5] + ((-0-1j))*g*hbar*x[7], 
(-1j)*g*hbar*x[11]*x[8] + (-1j)*g*hbar*x[1]*x[7] + (-1j)*g*hbar*x[7]*x[1] + (2j)*g*hbar*x[7]*x[7]*x[8] + (1j)*g*hbar*x[2]*np.conj(x[8]) + (1j)*g*hbar*np.conj(x[1])*x[7] + (1j)*g*hbar*np.conj(x[7])*x[10] + ((-0-2j))*g*hbar*np.conj(x[7])*x[7]*np.conj(x[8]) + (-2j)*g*hbar*x[4]*np.conj(x[8]) + (-2j)*g*hbar*x[5]*x[3] + (-2j)*g*hbar*x[3]*x[5] + (4j)*g*hbar*x[3]*x[3]*np.conj(x[8]) + (1.0)*gamma*hbar*x[2]*x[7] + (1.0)*gamma*hbar*x[2]*x[7] + (1.0)*gamma*hbar*np.conj(x[7])*x[11] + (-2.0)*gamma*hbar*np.conj(x[7])*x[7]*x[7] + (-1.0)*gamma*hbar*x[4]*x[7] + (-1.0)*gamma*hbar*x[6]*x[3] + (-1.0)*gamma*hbar*x[3]*x[6] + (2.0)*gamma*hbar*x[3]*x[3]*x[7] + (-1.0)*nu*hbar*x[2]*x[7] + (-1.0)*nu*hbar*x[2]*x[7] + (-1.0)*nu*hbar*np.conj(x[7])*x[11] + (2.0)*nu*hbar*np.conj(x[7])*x[7]*x[7] + (-5.0)*nu*hbar*hbar*x[6] + (2.0)*nu*hbar*hbar*hbar*x[7] + (1.0)*nu*hbar*x[4]*x[7] + (1.0)*nu*hbar*x[6]*x[3] + (1.0)*nu*hbar*x[3]*x[6] + (-2.0)*nu*hbar*x[3]*x[3]*x[7], 
(-2j)*g*hbar*x[5] + (-1.0)*gamma*hbar*x[6] + (1.0)*nu*hbar*x[6] + (-1.0)*nu*hbar*hbar*x[7], 
(-1j)*delta*x[8] + (-1j)*g*np.conj(x[7]) + (-0.5)*kappa*x[8], 
(-1j)*g*hbar*x[2]*x[8] + (-1j)*g*hbar*np.conj(x[10])*x[7] + (-1j)*g*hbar*np.conj(x[7])*x[1] + (2j)*g*hbar*np.conj(x[7])*x[7]*x[8] + (-2j)*g*hbar*hbar*np.conj(x[5]) + (2j)*g*hbar*x[4]*x[8] + (2j)*g*hbar*np.conj(x[5])*x[3] + (2j)*g*hbar*x[3]*np.conj(x[5]) + ((-0-4j))*g*hbar*x[3]*x[3]*x[8] + (1j)*g*hbar*np.conj(x[11])*np.conj(x[8]) + (1j)*g*hbar*np.conj(x[1])*np.conj(x[7]) + (1j)*g*hbar*np.conj(x[7])*np.conj(x[1]) + ((-0-2j))*g*hbar*np.conj(x[7])*np.conj(x[7])*np.conj(x[8]) + (1.0)*gamma*hbar*np.conj(x[11])*x[7] + (1.0)*gamma*hbar*x[2]*np.conj(x[7]) + (1.0)*gamma*hbar*np.conj(x[7])*x[2] + (-2.0)*gamma*hbar*np.conj(x[7])*np.conj(x[7])*x[7] + (-1.0)*gamma*hbar*hbar*x[9] + (-1.0)*gamma*hbar*x[4]*np.conj(x[7]) + (-1.0)*gamma*hbar*x[9]*x[3] + (-1.0)*gamma*hbar*x[3]*x[9] + (2.0)*gamma*hbar*x[3]*x[3]*np.conj(x[7]) + (-1.0)*nu*hbar*np.conj(x[11])*x[7] + (-1.0)*nu*hbar*x[2]*np.conj(x[7]) + (-1.0)*nu*hbar*np.conj(x[7])*x[2] + (2.0)*nu*hbar*np.conj(x[7])*np.conj(x[7])*x[7] + (-4.0)*nu*hbar*hbar*x[9] + (-2.0)*nu*hbar*hbar*hbar*np.conj(x[7]) + (1.0)*nu*hbar*x[4]*np.conj(x[7]) + (1.0)*nu*hbar*x[9]*x[3] + (1.0)*nu*hbar*x[3]*x[9] + (-2.0)*nu*hbar*x[3]*x[3]*np.conj(x[7]), 
(1j)*delta*x[10] + (1j)*g*x[11] + (-0.5)*kappa*x[10] + (-2j)*g*hbar*x[5]*np.conj(x[8]) + (-2j)*g*hbar*x[5]*np.conj(x[8]) + (-2j)*g*hbar*x[3]*x[12] + (4j)*g*hbar*x[3]*np.conj(x[8])*np.conj(x[8]) + (-1.0)*gamma*hbar*x[6]*np.conj(x[8]) + (-1.0)*gamma*hbar*x[5]*x[7] + (-1.0)*gamma*hbar*x[3]*x[10] + (2.0)*gamma*hbar*x[3]*x[7]*np.conj(x[8]) + (1.0)*nu*hbar*x[6]*np.conj(x[8]) + (1.0)*nu*hbar*x[5]*x[7] + (1.0)*nu*hbar*x[3]*x[10] + (-2.0)*nu*hbar*x[3]*x[7]*np.conj(x[8]) + (-1.0)*nu*hbar*hbar*x[10], 
(-4j)*g*hbar*x[6]*np.conj(x[8]) + (-4j)*g*hbar*x[5]*x[7] + (-4j)*g*hbar*x[3]*x[10] + (8j)*g*hbar*x[3]*x[7]*np.conj(x[8]) + (2j)*g*hbar*hbar*x[10] + (-2.0)*gamma*hbar*x[6]*x[7] + (-2.0)*gamma*hbar*x[6]*x[7] + (-2.0)*gamma*hbar*x[3]*x[11] + (4.0)*gamma*hbar*x[3]*x[7]*x[7] + (1.0)*gamma*hbar*hbar*x[11] + (2.0)*nu*hbar*x[6]*x[7] + (2.0)*nu*hbar*x[6]*x[7] + (2.0)*nu*hbar*x[3]*x[11] + (-4.0)*nu*hbar*x[3]*x[7]*x[7] + (-3.0)*nu*hbar*hbar*x[11], 
(2j)*delta*x[12] + (2j)*g*x[10] + (-1.0)*kappa*x[12] ]

def corr_get_init_vec_ada_order_2(init_state):
	return [
(1.0 + 0.0j) * mean_value(init_state, op_ad * op_a), # 0
(1.0 + 0.0j) * mean_value(init_state, op_sp * op_a), # 1
(1.0 + 0.0j) * mean_value(init_state, op_sm * op_sp), # 2
(1.0 + 0.0j) * mean_value(init_state, op_sz), # 3
(1.0 + 0.0j) * mean_value(init_state, op_sz * op_sz), # 4
(1.0 + 0.0j) * mean_value(init_state, op_sz * op_ad), # 5
(1.0 + 0.0j) * mean_value(init_state, op_sz * op_sp), # 6
(1.0 + 0.0j) * mean_value(init_state, op_sp), # 7
(1.0 + 0.0j) * mean_value(init_state, op_a), # 8
(1.0 + 0.0j) * mean_value(init_state, op_sz * op_sm), # 9
(1.0 + 0.0j) * mean_value(init_state, op_sp * op_ad), # 10
(1.0 + 0.0j) * mean_value(init_state, op_sp * op_sp), # 11
(1.0 + 0.0j) * mean_value(init_state, op_ad * op_ad), # 12
]