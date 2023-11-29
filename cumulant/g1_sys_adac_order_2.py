import numpy as np
from cavity_global import *

def g1_system_adac_order_2(t, x, mean_ac, mean_adc):
	return [ (-2j)*delta*x[0] + (-2j)*g*x[1] + (-1.0)*kappa*x[0], 
(-1j)*delta*x[1] + (-1j)*g*x[2] + (-0.5)*kappa*x[1] + (2j)*g*hbar*x[6]*x[13] + (2j)*g*hbar*x[6]*x[13] + (2j)*g*hbar*x[11]*x[0] + ((-0-4j))*g*hbar*x[11]*x[13]*x[13] + (-1.0)*gamma*hbar*x[5]*x[13] + (-1.0)*gamma*hbar*x[6]*x[7] + (-1.0)*gamma*hbar*x[11]*x[1] + (2.0)*gamma*hbar*x[11]*x[7]*x[13] + (-1.0)*gamma*hbar*hbar*x[1] + (1.0)*nu*hbar*x[5]*x[13] + (1.0)*nu*hbar*x[6]*x[7] + (1.0)*nu*hbar*x[11]*x[1] + (-2.0)*nu*hbar*x[11]*x[7]*x[13], 
(4j)*g*hbar*x[5]*x[13] + (4j)*g*hbar*x[6]*x[7] + (4j)*g*hbar*x[11]*x[1] + ((-0-8j))*g*hbar*x[11]*x[7]*x[13] + (2j)*g*hbar*hbar*x[1] + (-2.0)*gamma*hbar*x[5]*x[7] + (-2.0)*gamma*hbar*x[5]*x[7] + (-2.0)*gamma*hbar*x[11]*x[2] + (4.0)*gamma*hbar*x[11]*x[7]*x[7] + (-3.0)*gamma*hbar*hbar*x[2] + (2.0)*nu*hbar*x[5]*x[7] + (2.0)*nu*hbar*x[5]*x[7] + (2.0)*nu*hbar*x[11]*x[2] + (-4.0)*nu*hbar*x[11]*x[7]*x[7] + (1.0)*nu*hbar*hbar*x[2], 
(2j)*g*hbar*x[6]*mean_ac + (2j)*g*hbar*x[15]*x[13] + (2j)*g*hbar*x[11]*x[4] + ((-0-4j))*g*hbar*x[11]*x[13]*mean_ac + (-1.0)*gamma*hbar*x[5]*mean_ac + (-1.0)*gamma*hbar*x[15]*x[7] + (-1.0)*gamma*hbar*x[11]*x[3] + (2.0)*gamma*hbar*x[11]*x[7]*mean_ac + (-1.0)*gamma*hbar*hbar*x[3] + (1.0)*nu*hbar*x[5]*mean_ac + (1.0)*nu*hbar*x[15]*x[7] + (1.0)*nu*hbar*x[11]*x[3] + (-2.0)*nu*hbar*x[11]*x[7]*mean_ac, 
(-1j)*delta*x[4] + (-1j)*g*x[3] + (-0.5)*kappa*x[4], 
(-1j)*g*hbar*x[10]*x[13] + (-1j)*g*hbar*x[1]*np.conj(x[7]) + (-1j)*g*hbar*x[7]*x[9] + (2j)*g*hbar*x[7]*np.conj(x[7])*x[13] + (-2j)*g*hbar*hbar*x[6] + (2j)*g*hbar*x[12]*x[13] + (2j)*g*hbar*x[6]*x[11] + (2j)*g*hbar*x[11]*x[6] + ((-0-4j))*g*hbar*x[11]*x[11]*x[13] + (1j)*g*hbar*x[2]*np.conj(x[13]) + (1j)*g*hbar*np.conj(x[9])*x[7] + (1j)*g*hbar*x[7]*np.conj(x[9]) + ((-0-2j))*g*hbar*x[7]*x[7]*np.conj(x[13]) + (1.0)*gamma*hbar*x[2]*np.conj(x[7]) + (1.0)*gamma*hbar*x[10]*x[7] + (1.0)*gamma*hbar*x[7]*x[10] + (-2.0)*gamma*hbar*x[7]*x[7]*np.conj(x[7]) + (-1.0)*gamma*hbar*hbar*x[5] + (-1.0)*gamma*hbar*x[12]*x[7] + (-1.0)*gamma*hbar*x[5]*x[11] + (-1.0)*gamma*hbar*x[11]*x[5] + (2.0)*gamma*hbar*x[11]*x[11]*x[7] + (-1.0)*nu*hbar*x[2]*np.conj(x[7]) + (-1.0)*nu*hbar*x[10]*x[7] + (-1.0)*nu*hbar*x[7]*x[10] + (2.0)*nu*hbar*x[7]*x[7]*np.conj(x[7]) + (-4.0)*nu*hbar*hbar*x[5] + (-2.0)*nu*hbar*hbar*hbar*x[7] + (1.0)*nu*hbar*x[12]*x[7] + (1.0)*nu*hbar*x[5]*x[11] + (1.0)*nu*hbar*x[11]*x[5] + (-2.0)*nu*hbar*x[11]*x[11]*x[7], 
(-1j)*delta*x[6] + (-1j)*g*x[5] + (-0.5)*kappa*x[6] + (-1j)*g*hbar*x[9]*x[13] + (-1j)*g*hbar*x[9]*x[13] + (-1j)*g*hbar*np.conj(x[7])*x[0] + (2j)*g*hbar*np.conj(x[7])*x[13]*x[13] + (1j)*g*hbar*np.conj(x[9])*x[13] + (1j)*g*hbar*x[1]*np.conj(x[13]) + (1j)*g*hbar*x[7]*x[8] + ((-0-2j))*g*hbar*x[7]*np.conj(x[13])*x[13] + (1.0)*gamma*hbar*x[10]*x[13] + (1.0)*gamma*hbar*x[1]*np.conj(x[7]) + (1.0)*gamma*hbar*x[7]*x[9] + (-2.0)*gamma*hbar*x[7]*np.conj(x[7])*x[13] + (-1.0)*nu*hbar*x[10]*x[13] + (-1.0)*nu*hbar*x[1]*np.conj(x[7]) + (-1.0)*nu*hbar*x[7]*x[9] + (2.0)*nu*hbar*x[7]*np.conj(x[7])*x[13] + (-2.0)*nu*hbar*hbar*x[6], 
(2j)*g*hbar*x[6] + (-1.0)*gamma*hbar*x[5] + (-1.0)*gamma*hbar*hbar*x[7] + (1.0)*nu*hbar*x[5], 
(1j)*g*x[9] + (-1j)*g*np.conj(x[9]) + (-1.0)*kappa*x[8], 
(-1j)*delta*x[9] + (-1j)*g*x[10] + (-0.5)*kappa*x[9] + (-2j)*g*hbar*np.conj(x[6])*x[13] + (-2j)*g*hbar*x[6]*np.conj(x[13]) + (-2j)*g*hbar*x[11]*x[8] + (4j)*g*hbar*x[11]*np.conj(x[13])*x[13] + (-1.0)*gamma*hbar*x[14]*x[13] + (-1.0)*gamma*hbar*x[6]*np.conj(x[7]) + (-1.0)*gamma*hbar*x[11]*x[9] + (2.0)*gamma*hbar*x[11]*np.conj(x[7])*x[13] + (1.0)*nu*hbar*x[14]*x[13] + (1.0)*nu*hbar*x[6]*np.conj(x[7]) + (1.0)*nu*hbar*x[11]*x[9] + (-2.0)*nu*hbar*x[11]*np.conj(x[7])*x[13] + (-1.0)*nu*hbar*hbar*x[9] + (-2j)*g*hbar*x[11], 
(2j)*g*hbar*x[14]*x[13] + (2j)*g*hbar*x[6]*np.conj(x[7]) + (2j)*g*hbar*x[11]*x[9] + ((-0-4j))*g*hbar*x[11]*np.conj(x[7])*x[13] + (-2j)*g*hbar*x[5]*np.conj(x[13]) + (-2j)*g*hbar*np.conj(x[6])*x[7] + (-2j)*g*hbar*x[11]*np.conj(x[9]) + (4j)*g*hbar*x[11]*x[7]*np.conj(x[13]) + (-2j)*g*hbar*hbar*np.conj(x[9]) + (-2.0)*gamma*hbar*x[5]*np.conj(x[7]) + (-2.0)*gamma*hbar*x[14]*x[7] + (-2.0)*gamma*hbar*x[11]*x[10] + (4.0)*gamma*hbar*x[11]*x[7]*np.conj(x[7]) + (-2.0)*gamma*hbar*hbar*x[10] + (2.0)*nu*hbar*x[5]*np.conj(x[7]) + (2.0)*nu*hbar*x[14]*x[7] + (2.0)*nu*hbar*x[11]*x[10] + (-4.0)*nu*hbar*x[11]*x[7]*np.conj(x[7]) + (4.0)*nu*hbar*hbar*x[12], 
((-0-1j))*g*hbar*x[9] + (1j)*g*hbar*np.conj(x[9]) + (1.0)*gamma*hbar*x[10] + (-1.0)*nu*hbar*x[10] + (-2.0)*nu*hbar*hbar*x[11], 
(-2j)*g*hbar*x[14]*x[13] + (-2j)*g*hbar*x[6]*np.conj(x[7]) + (-2j)*g*hbar*x[11]*x[9] + (4j)*g*hbar*x[11]*np.conj(x[7])*x[13] + (1j)*g*hbar*hbar*x[9] + (2j)*g*hbar*x[5]*np.conj(x[13]) + (2j)*g*hbar*np.conj(x[6])*x[7] + (2j)*g*hbar*x[11]*np.conj(x[9]) + ((-0-4j))*g*hbar*x[11]*x[7]*np.conj(x[13]) + (1j)*g*hbar*hbar*np.conj(x[9]) + (2.0)*gamma*hbar*x[5]*np.conj(x[7]) + (2.0)*gamma*hbar*x[14]*x[7] + (2.0)*gamma*hbar*x[11]*x[10] + (-4.0)*gamma*hbar*x[11]*x[7]*np.conj(x[7]) + (1.0)*gamma*hbar*hbar*x[10] + (-2.0)*nu*hbar*x[5]*np.conj(x[7]) + (-2.0)*nu*hbar*x[14]*x[7] + (-2.0)*nu*hbar*x[11]*x[10] + (4.0)*nu*hbar*x[11]*x[7]*np.conj(x[7]) + (1.0)*nu*hbar*hbar*x[10] + (2.0)*nu*hbar*hbar*hbar*x[11] + (-4.0)*nu*hbar*hbar*x[12], 
(-1j)*delta*x[13] + (-1j)*g*x[7] + (-0.5)*kappa*x[13], 
(-1j)*g*hbar*np.conj(x[2])*x[13] + (-1j)*g*hbar*x[9]*np.conj(x[7]) + (-1j)*g*hbar*np.conj(x[7])*x[9] + (2j)*g*hbar*np.conj(x[7])*np.conj(x[7])*x[13] + (1j)*g*hbar*x[10]*np.conj(x[13]) + (1j)*g*hbar*np.conj(x[9])*np.conj(x[7]) + (1j)*g*hbar*x[7]*np.conj(x[1]) + ((-0-2j))*g*hbar*x[7]*np.conj(x[7])*np.conj(x[13]) + (-2j)*g*hbar*x[12]*np.conj(x[13]) + (-2j)*g*hbar*np.conj(x[6])*x[11] + (-2j)*g*hbar*x[11]*np.conj(x[6]) + (4j)*g*hbar*x[11]*x[11]*np.conj(x[13]) + (1.0)*gamma*hbar*x[10]*np.conj(x[7]) + (1.0)*gamma*hbar*x[10]*np.conj(x[7]) + (1.0)*gamma*hbar*x[7]*np.conj(x[2]) + (-2.0)*gamma*hbar*x[7]*np.conj(x[7])*np.conj(x[7]) + (-1.0)*gamma*hbar*x[12]*np.conj(x[7]) + (-1.0)*gamma*hbar*x[14]*x[11] + (-1.0)*gamma*hbar*x[11]*x[14] + (2.0)*gamma*hbar*x[11]*x[11]*np.conj(x[7]) + (-1.0)*nu*hbar*x[10]*np.conj(x[7]) + (-1.0)*nu*hbar*x[10]*np.conj(x[7]) + (-1.0)*nu*hbar*x[7]*np.conj(x[2]) + (2.0)*nu*hbar*x[7]*np.conj(x[7])*np.conj(x[7]) + (-5.0)*nu*hbar*hbar*x[14] + (2.0)*nu*hbar*hbar*hbar*np.conj(x[7]) + (1.0)*nu*hbar*x[12]*np.conj(x[7]) + (1.0)*nu*hbar*x[14]*x[11] + (1.0)*nu*hbar*x[11]*x[14] + (-2.0)*nu*hbar*x[11]*x[11]*np.conj(x[7]), 
(-1j)*g*hbar*x[9]*mean_ac + (-1j)*g*hbar*x[17]*x[13] + (-1j)*g*hbar*np.conj(x[7])*x[4] + (2j)*g*hbar*np.conj(x[7])*x[13]*mean_ac + (1j)*g*hbar*np.conj(x[9])*mean_ac + (1j)*g*hbar*x[3]*np.conj(x[13]) + (1j)*g*hbar*x[7]*x[16] + ((-0-2j))*g*hbar*x[7]*np.conj(x[13])*mean_ac + (1.0)*gamma*hbar*x[10]*mean_ac + (1.0)*gamma*hbar*x[3]*np.conj(x[7]) + (1.0)*gamma*hbar*x[7]*x[17] + (-2.0)*gamma*hbar*x[7]*np.conj(x[7])*mean_ac + (-1.0)*nu*hbar*x[10]*mean_ac + (-1.0)*nu*hbar*x[3]*np.conj(x[7]) + (-1.0)*nu*hbar*x[7]*x[17] + (2.0)*nu*hbar*x[7]*np.conj(x[7])*mean_ac + (-2.0)*nu*hbar*hbar*x[15], 
(1j)*delta*x[16] + (1j)*g*x[17] + (-0.5)*kappa*x[16], 
(-2j)*g*hbar*np.conj(x[6])*mean_ac + (-2j)*g*hbar*x[15]*np.conj(x[13]) + (-2j)*g*hbar*x[11]*x[16] + (4j)*g*hbar*x[11]*np.conj(x[13])*mean_ac + (-1.0)*gamma*hbar*x[14]*mean_ac + (-1.0)*gamma*hbar*x[15]*np.conj(x[7]) + (-1.0)*gamma*hbar*x[11]*x[17] + (2.0)*gamma*hbar*x[11]*np.conj(x[7])*mean_ac + (1.0)*nu*hbar*x[14]*mean_ac + (1.0)*nu*hbar*x[15]*np.conj(x[7]) + (1.0)*nu*hbar*x[11]*x[17] + (-2.0)*nu*hbar*x[11]*np.conj(x[7])*mean_ac + (-1.0)*nu*hbar*hbar*x[17] ]

def g1_get_init_vec_adac_order_2(x, t_index):
	return [
x[2][t_index],# aa
x[1][t_index],# sma
x[0][t_index],# smsm
x[1][t_index],# smac
x[2][t_index],# aac
x[9][t_index],# szsm
x[8][t_index],# sza
x[7][t_index],# sm
x[10][t_index],# ada
x[4][t_index],# spa
x[5][t_index],# smsp
x[3][t_index],# sz
x[6][t_index],# szsz
np.conj(x[11][t_index]),# a
x[12][t_index],# szsp
x[8][t_index],# szac
x[10][t_index],# adac
x[4][t_index],# spac
]