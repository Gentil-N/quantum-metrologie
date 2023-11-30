import numpy as np
from cavity_global import *

def g1_system_adac_order_2(t, x, mean_ac, mean_adc):
	return [ (1j)*delta*x[0] + (1j)*g*x[1] + (-0.5)*kappa*x[0], 
(-2j)*g*hbar*x[2]*mean_ac + (-2j)*g*hbar*x[5]*x[6] + (-2j)*g*hbar*x[7]*x[0] + (4j)*g*hbar*x[7]*x[6]*mean_ac + (-1.0)*gamma*hbar*x[3]*mean_ac + (-1.0)*gamma*hbar*x[5]*x[4] + (-1.0)*gamma*hbar*x[7]*x[1] + (2.0)*gamma*hbar*x[7]*x[4]*mean_ac + (1.0)*nu*hbar*x[3]*mean_ac + (1.0)*nu*hbar*x[5]*x[4] + (1.0)*nu*hbar*x[7]*x[1] + (-2.0)*nu*hbar*x[7]*x[4]*mean_ac + (-1.0)*nu*hbar*hbar*x[1], 
(1j)*delta*x[2] + (1j)*g*x[3] + (-0.5)*kappa*x[2] + (-1j)*g*hbar*x[11]*np.conj(x[6]) + (-1j)*g*hbar*x[8]*x[6] + (-1j)*g*hbar*x[4]*x[14] + (2j)*g*hbar*x[4]*x[6]*np.conj(x[6]) + (1j)*g*hbar*np.conj(x[8])*x[6] + (1j)*g*hbar*np.conj(x[8])*x[6] + (1j)*g*hbar*np.conj(x[4])*x[15] + ((-0-2j))*g*hbar*np.conj(x[4])*x[6]*x[6] + (1.0)*gamma*hbar*x[9]*x[6] + (1.0)*gamma*hbar*np.conj(x[8])*x[4] + (1.0)*gamma*hbar*np.conj(x[4])*x[11] + (-2.0)*gamma*hbar*np.conj(x[4])*x[4]*x[6] + (-1.0)*nu*hbar*x[9]*x[6] + (-1.0)*nu*hbar*np.conj(x[8])*x[4] + (-1.0)*nu*hbar*np.conj(x[4])*x[11] + (2.0)*nu*hbar*np.conj(x[4])*x[4]*x[6] + (-2.0)*nu*hbar*hbar*x[2] + ((-0-1j))*g*hbar*x[4], 
(-1j)*g*hbar*x[12]*np.conj(x[6]) + (-1j)*g*hbar*x[8]*x[4] + (-1j)*g*hbar*x[4]*x[8] + (2j)*g*hbar*x[4]*x[4]*np.conj(x[6]) + (1j)*g*hbar*x[9]*x[6] + (1j)*g*hbar*np.conj(x[8])*x[4] + (1j)*g*hbar*np.conj(x[4])*x[11] + ((-0-2j))*g*hbar*np.conj(x[4])*x[4]*x[6] + (-2j)*g*hbar*x[10]*x[6] + (-2j)*g*hbar*x[2]*x[7] + (-2j)*g*hbar*x[7]*x[2] + (4j)*g*hbar*x[7]*x[7]*x[6] + (1.0)*gamma*hbar*x[9]*x[4] + (1.0)*gamma*hbar*x[9]*x[4] + (1.0)*gamma*hbar*np.conj(x[4])*x[12] + (-2.0)*gamma*hbar*np.conj(x[4])*x[4]*x[4] + (-1.0)*gamma*hbar*x[10]*x[4] + (-1.0)*gamma*hbar*x[3]*x[7] + (-1.0)*gamma*hbar*x[7]*x[3] + (2.0)*gamma*hbar*x[7]*x[7]*x[4] + (-1.0)*nu*hbar*x[9]*x[4] + (-1.0)*nu*hbar*x[9]*x[4] + (-1.0)*nu*hbar*np.conj(x[4])*x[12] + (2.0)*nu*hbar*np.conj(x[4])*x[4]*x[4] + (-5.0)*nu*hbar*hbar*x[3] + (2.0)*nu*hbar*hbar*hbar*x[4] + (1.0)*nu*hbar*x[10]*x[4] + (1.0)*nu*hbar*x[3]*x[7] + (1.0)*nu*hbar*x[7]*x[3] + (-2.0)*nu*hbar*x[7]*x[7]*x[4], 
(-2j)*g*hbar*x[2] + (-1.0)*gamma*hbar*x[3] + (1.0)*nu*hbar*x[3] + (-1.0)*nu*hbar*hbar*x[4], 
(-1j)*g*hbar*x[8]*mean_ac + (-1j)*g*hbar*x[1]*np.conj(x[6]) + (-1j)*g*hbar*x[4]*x[16] + (2j)*g*hbar*x[4]*np.conj(x[6])*mean_ac + (1j)*g*hbar*np.conj(x[8])*mean_ac + (1j)*g*hbar*x[17]*x[6] + (1j)*g*hbar*np.conj(x[4])*x[0] + ((-0-2j))*g*hbar*np.conj(x[4])*x[6]*mean_ac + (1.0)*gamma*hbar*x[9]*mean_ac + (1.0)*gamma*hbar*x[17]*x[4] + (1.0)*gamma*hbar*np.conj(x[4])*x[1] + (-2.0)*gamma*hbar*np.conj(x[4])*x[4]*mean_ac + (-1.0)*nu*hbar*x[9]*mean_ac + (-1.0)*nu*hbar*x[17]*x[4] + (-1.0)*nu*hbar*np.conj(x[4])*x[1] + (2.0)*nu*hbar*np.conj(x[4])*x[4]*mean_ac + (-2.0)*nu*hbar*hbar*x[5], 
(1j)*delta*x[6] + (1j)*g*x[4] + (-0.5)*kappa*x[6], 
((-0-1j))*g*hbar*x[8] + (1j)*g*hbar*np.conj(x[8]) + (1.0)*gamma*hbar*x[9] + (-1.0)*nu*hbar*x[9] + (-2.0)*nu*hbar*hbar*x[7], 
(-1j)*delta*x[8] + (-1j)*g*x[9] + (-0.5)*kappa*x[8] + (-2j)*g*hbar*x[2]*np.conj(x[6]) + (-2j)*g*hbar*np.conj(x[2])*x[6] + (-2j)*g*hbar*x[7]*x[14] + (4j)*g*hbar*x[7]*x[6]*np.conj(x[6]) + (-1.0)*gamma*hbar*x[3]*np.conj(x[6]) + (-1.0)*gamma*hbar*np.conj(x[2])*x[4] + (-1.0)*gamma*hbar*x[7]*x[8] + (2.0)*gamma*hbar*x[7]*x[4]*np.conj(x[6]) + (1.0)*nu*hbar*x[3]*np.conj(x[6]) + (1.0)*nu*hbar*np.conj(x[2])*x[4] + (1.0)*nu*hbar*x[7]*x[8] + (-2.0)*nu*hbar*x[7]*x[4]*np.conj(x[6]) + (-1.0)*nu*hbar*hbar*x[8] + (-2j)*g*hbar*x[7], 
(2j)*g*hbar*x[3]*np.conj(x[6]) + (2j)*g*hbar*np.conj(x[2])*x[4] + (2j)*g*hbar*x[7]*x[8] + ((-0-4j))*g*hbar*x[7]*x[4]*np.conj(x[6]) + (-2j)*g*hbar*x[13]*x[6] + (-2j)*g*hbar*x[2]*np.conj(x[4]) + (-2j)*g*hbar*x[7]*np.conj(x[8]) + (4j)*g*hbar*x[7]*np.conj(x[4])*x[6] + (-2j)*g*hbar*hbar*np.conj(x[8]) + (-2.0)*gamma*hbar*x[13]*x[4] + (-2.0)*gamma*hbar*x[3]*np.conj(x[4]) + (-2.0)*gamma*hbar*x[7]*x[9] + (4.0)*gamma*hbar*x[7]*np.conj(x[4])*x[4] + (-2.0)*gamma*hbar*hbar*x[9] + (2.0)*nu*hbar*x[13]*x[4] + (2.0)*nu*hbar*x[3]*np.conj(x[4]) + (2.0)*nu*hbar*x[7]*x[9] + (-4.0)*nu*hbar*x[7]*np.conj(x[4])*x[4] + (4.0)*nu*hbar*hbar*x[10], 
(-2j)*g*hbar*x[3]*np.conj(x[6]) + (-2j)*g*hbar*np.conj(x[2])*x[4] + (-2j)*g*hbar*x[7]*x[8] + (4j)*g*hbar*x[7]*x[4]*np.conj(x[6]) + (1j)*g*hbar*hbar*x[8] + (2j)*g*hbar*x[13]*x[6] + (2j)*g*hbar*x[2]*np.conj(x[4]) + (2j)*g*hbar*x[7]*np.conj(x[8]) + ((-0-4j))*g*hbar*x[7]*np.conj(x[4])*x[6] + (1j)*g*hbar*hbar*np.conj(x[8]) + (2.0)*gamma*hbar*x[13]*x[4] + (2.0)*gamma*hbar*x[3]*np.conj(x[4]) + (2.0)*gamma*hbar*x[7]*x[9] + (-4.0)*gamma*hbar*x[7]*np.conj(x[4])*x[4] + (1.0)*gamma*hbar*hbar*x[9] + (-2.0)*nu*hbar*x[13]*x[4] + (-2.0)*nu*hbar*x[3]*np.conj(x[4]) + (-2.0)*nu*hbar*x[7]*x[9] + (4.0)*nu*hbar*x[7]*np.conj(x[4])*x[4] + (1.0)*nu*hbar*hbar*x[9] + (2.0)*nu*hbar*hbar*hbar*x[7] + (-4.0)*nu*hbar*hbar*x[10], 
(1j)*delta*x[11] + (1j)*g*x[12] + (-0.5)*kappa*x[11] + (-2j)*g*hbar*x[2]*x[6] + (-2j)*g*hbar*x[2]*x[6] + (-2j)*g*hbar*x[7]*x[15] + (4j)*g*hbar*x[7]*x[6]*x[6] + (-1.0)*gamma*hbar*x[3]*x[6] + (-1.0)*gamma*hbar*x[2]*x[4] + (-1.0)*gamma*hbar*x[7]*x[11] + (2.0)*gamma*hbar*x[7]*x[4]*x[6] + (1.0)*nu*hbar*x[3]*x[6] + (1.0)*nu*hbar*x[2]*x[4] + (1.0)*nu*hbar*x[7]*x[11] + (-2.0)*nu*hbar*x[7]*x[4]*x[6] + (-1.0)*nu*hbar*hbar*x[11], 
(-4j)*g*hbar*x[3]*x[6] + (-4j)*g*hbar*x[2]*x[4] + (-4j)*g*hbar*x[7]*x[11] + (8j)*g*hbar*x[7]*x[4]*x[6] + (2j)*g*hbar*hbar*x[11] + (-2.0)*gamma*hbar*x[3]*x[4] + (-2.0)*gamma*hbar*x[3]*x[4] + (-2.0)*gamma*hbar*x[7]*x[12] + (4.0)*gamma*hbar*x[7]*x[4]*x[4] + (1.0)*gamma*hbar*hbar*x[12] + (2.0)*nu*hbar*x[3]*x[4] + (2.0)*nu*hbar*x[3]*x[4] + (2.0)*nu*hbar*x[7]*x[12] + (-4.0)*nu*hbar*x[7]*x[4]*x[4] + (-3.0)*nu*hbar*hbar*x[12], 
(-1j)*g*hbar*x[9]*np.conj(x[6]) + (-1j)*g*hbar*np.conj(x[11])*x[4] + (-1j)*g*hbar*np.conj(x[4])*x[8] + (2j)*g*hbar*np.conj(x[4])*x[4]*np.conj(x[6]) + (-2j)*g*hbar*hbar*np.conj(x[2]) + (2j)*g*hbar*x[10]*np.conj(x[6]) + (2j)*g*hbar*np.conj(x[2])*x[7] + (2j)*g*hbar*x[7]*np.conj(x[2]) + ((-0-4j))*g*hbar*x[7]*x[7]*np.conj(x[6]) + (1j)*g*hbar*np.conj(x[12])*x[6] + (1j)*g*hbar*np.conj(x[8])*np.conj(x[4]) + (1j)*g*hbar*np.conj(x[4])*np.conj(x[8]) + ((-0-2j))*g*hbar*np.conj(x[4])*np.conj(x[4])*x[6] + (1.0)*gamma*hbar*np.conj(x[12])*x[4] + (1.0)*gamma*hbar*x[9]*np.conj(x[4]) + (1.0)*gamma*hbar*np.conj(x[4])*x[9] + (-2.0)*gamma*hbar*np.conj(x[4])*np.conj(x[4])*x[4] + (-1.0)*gamma*hbar*hbar*x[13] + (-1.0)*gamma*hbar*x[10]*np.conj(x[4]) + (-1.0)*gamma*hbar*x[13]*x[7] + (-1.0)*gamma*hbar*x[7]*x[13] + (2.0)*gamma*hbar*x[7]*x[7]*np.conj(x[4]) + (-1.0)*nu*hbar*np.conj(x[12])*x[4] + (-1.0)*nu*hbar*x[9]*np.conj(x[4]) + (-1.0)*nu*hbar*np.conj(x[4])*x[9] + (2.0)*nu*hbar*np.conj(x[4])*np.conj(x[4])*x[4] + (-4.0)*nu*hbar*hbar*x[13] + (-2.0)*nu*hbar*hbar*hbar*np.conj(x[4]) + (1.0)*nu*hbar*x[10]*np.conj(x[4]) + (1.0)*nu*hbar*x[13]*x[7] + (1.0)*nu*hbar*x[7]*x[13] + (-2.0)*nu*hbar*x[7]*x[7]*np.conj(x[4]), 
(1j)*g*x[8] + (-1j)*g*np.conj(x[8]) + (-1.0)*kappa*x[14], 
(2j)*delta*x[15] + (2j)*g*x[11] + (-1.0)*kappa*x[15], 
(-1j)*delta*x[16] + (-1j)*g*x[17] + (-0.5)*kappa*x[16], 
(2j)*g*hbar*np.conj(x[2])*mean_ac + (2j)*g*hbar*x[5]*np.conj(x[6]) + (2j)*g*hbar*x[7]*x[16] + ((-0-4j))*g*hbar*x[7]*np.conj(x[6])*mean_ac + (-1.0)*gamma*hbar*x[13]*mean_ac + (-1.0)*gamma*hbar*x[5]*np.conj(x[4]) + (-1.0)*gamma*hbar*x[7]*x[17] + (2.0)*gamma*hbar*x[7]*np.conj(x[4])*mean_ac + (-1.0)*gamma*hbar*hbar*x[17] + (1.0)*nu*hbar*x[13]*mean_ac + (1.0)*nu*hbar*x[5]*np.conj(x[4]) + (1.0)*nu*hbar*x[7]*x[17] + (-2.0)*nu*hbar*x[7]*np.conj(x[4])*mean_ac ]

def g1_get_init_vec_adac_order_2(x, t_index):
	return [
x[11][t_index],# adac # 0
x[7][t_index],# spac # 1
x[2][t_index],# szad # 2
x[3][t_index],# szsp # 3
x[1][t_index],# sp # 4
np.conj(x[2][t_index]),# szac # 5
x[0][t_index],# ad # 6
x[9][t_index],# sz # 7
x[7][t_index],# spa # 8
x[8][t_index],# smsp # 9
x[10][t_index],# szsz # 10
x[4][t_index],# spad # 11
x[5][t_index],# spsp # 12
x[6][t_index],# szsm # 13
x[11][t_index],# ada # 14
x[12][t_index],# adad # 15
np.conj(x[12][t_index]),# aac # 16
np.conj(x[4][t_index]),# smac # 17
]