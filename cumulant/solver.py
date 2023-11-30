import numpy as np
from scipy.integrate import solve_ivp
from matplotlib import pyplot
from matplotlib import pyplot as plt
import qutip as qp
from cavity_global import *
from corr_sys_ada_order_2 import *
from corr_sys_ada_order_3 import *
from corr_sys_ada_order_4 import *
from corr_sys_ad_order_2 import *
from corr_sys_ad_order_3 import *
from corr_sys_ad_order_4 import *
from g1_sys_adac_order_2 import *
from g1_sys_adac_order_3 import *
from g1_sys_adac_order_4 import *
from g2_sys_adcadaac_order_4 import *

#op_sz = qp.tensor(qp.spin_Jz(SPIN_NUM), qp.qeye(FOCK_DIM))
#op_sp = qp.tensor(qp.spin_Jp(SPIN_NUM), qp.qeye(FOCK_DIM))
#op_sm = qp.tensor(qp.spin_Jm(SPIN_NUM), qp.qeye(FOCK_DIM))
#op_ad = qp.tensor(qp.qeye(SPIN_DIM), qp.create(FOCK_DIM))
#op_a = qp.tensor(qp.qeye(SPIN_DIM), qp.destroy(FOCK_DIM))
#op_n = op_ad * op_a

#init_state = qp.tensor(qp.spin_state(ATOM_COUNT/2, ATOM_COUNT/2), qp.fock(PHOTON_CAPACITY + 1, 0))



### corr



#x_init_ada_order_2 = corr_get_init_vec_ada_order_2(init_state)
#x_init_ada_order_3 = corr_get_init_vec_ada_order_3(init_state)
#x_init_ada_order_4 = corr_get_init_vec_ada_order_4(init_state)

#print(x_init_order_3)
#exit()

#print("Solver started...", end='', flush=True)
#solution_order_2 = solve_ivp(corr_system_ada_order_2, [0.0, 10.0], x_init_ada_order_2, first_step=0.1, max_step=0.1)
#solution_order_3 = solve_ivp(corr_system_ada_order_3, [0.0, 10.0], x_init_ada_order_3, first_step=0.1, max_step=0.1)
#solution_order_4 = solve_ivp(corr_system_ada_order_4, [0.0, 10.0], x_init_ada_order_4, first_step=0.1, max_step=0.1)
#print("Done")

# /!\ /!\ /!\ If corr_sys file are regenrated, think to update the code to grab n and Sz: indices not updated for the last version of file generation!!!

#fig, ax = pyplot.subplots(1, 1)
#fig.set_size_inches(18.5, 10.5)
#for i in range(len(x_init)):
#    ax.plot(solution.t, np.real(solution.y[i]), label=str(i))
#ax.plot(solution.t, np.real(np.conj(solution.y[3])), label="<S+>")
#ax.plot(solution_order_2.t, np.real(solution_order_2.y[8]), label="<n> 2nd", color="green", linestyle="--")
#ax.plot(solution_order_2.t, np.real(solution_order_2.y[11]), label="<Sz> 2nd", color="green")
#ax.plot(solution_order_3.t, np.real(solution_order_3.y[26]), label="<n> 3rd", color="orange", linestyle="--")
#ax.plot(solution_order_3.t, np.real(solution_order_3.y[31]), label="<Sz> 3rd", color="orange")
#ax.plot(solution_order_4.t, np.real(solution_order_4.y[85]), label="<n> 4rd", color="red", linestyle="--")
#ax.plot(solution_order_4.t, np.real(solution_order_4.y[77]), label="<Sz> 4rd", color="red")

#ax.legend()
#ax.set_xlabel('time')
#ax.set_ylabel('correlations')

#plt.show()



### g1



#x_init_ad_order_2 = corr_get_init_vec_ad_order_2(init_state)
#x_init_ad_order_3 = corr_get_init_vec_ad_order_3(init_state)
#x_init_ad_order_4 = corr_get_init_vec_ad_order_4(init_state)
#
#print("Solver started...", end='', flush=True)
#solution_order_2 = solve_ivp(corr_system_ad_order_2, [0.0, 10.0], x_init_ad_order_2, first_step=0.1, max_step=0.1)
#solution_order_3 = solve_ivp(corr_system_ad_order_3, [0.0, 10.0], x_init_ad_order_3, first_step=0.1, max_step=0.1)
#solution_order_4 = solve_ivp(corr_system_ad_order_4, [0.0, 10.0], x_init_ad_order_4, first_step=0.1, max_step=0.1)
#print("Done")
#
## We compute g1 in the steady-state so we take the last value of a(t) and ad(t) for a0 (= ac) and ad0 (= adc)
## Note: for g1, we do not use ad0
#mean_ac_order_2 = np.conj(solution_order_2.y[0][-1])
#mean_adc_order_2 = solution_order_2.y[0][-1]
#mean_ac_order_3 = np.conj(solution_order_3.y[0][-1])
#mean_adc_order_3 = solution_order_3.y[0][-1]
#mean_ac_order_4 = np.conj(solution_order_4.y[0][-1])
#mean_adc_order_4 = solution_order_4.y[0][-1]
#
#def g1_sys_order_2(t, x):
#    return g1_system_adac_order_2(t, x, mean_ac_order_2, mean_adc_order_2)
#
#def g1_sys_order_3(t, x):
#    return g1_system_adac_order_3(t, x, mean_ac_order_3, mean_adc_order_3)
#
#def g1_sys_order_4(t, x):
#    return g1_system_adac_order_4(t, x, mean_ac_order_4, mean_adc_order_4)
#
#x_init_g1_order_2 = g1_get_init_vec_adac_order_2(solution_order_2.y, -1) # t_index = -1 : we take the last element of the list (steady-state)
#x_init_g1_order_3 = g1_get_init_vec_adac_order_3(solution_order_3.y, -1)
#x_init_g1_order_4 = g1_get_init_vec_adac_order_4(solution_order_4.y, -1)
#
#print("Solver started...", end='', flush=True)
#g1_solution_order_2 = solve_ivp(g1_sys_order_2, [0.0, 10.0], x_init_g1_order_2, first_step=0.1, max_step=0.1)
#g1_solution_order_3 = solve_ivp(g1_sys_order_3, [0.0, 10.0], x_init_g1_order_3, first_step=0.1, max_step=0.1)
#g1_solution_order_4 = solve_ivp(g1_sys_order_4, [0.0, 10.0], x_init_g1_order_4, first_step=0.1, max_step=0.1)
#print("Done")
#
#fig, ax = pyplot.subplots(1, 1)
#fig.set_size_inches(18.5, 10.5)
#ax.plot(g1_solution_order_2.t, np.real(g1_solution_order_2.y[0]), label="<ad(t)ac> 2nd steady-state", color="green", linestyle="--")
#ax.plot(g1_solution_order_3.t, np.real(g1_solution_order_3.y[0]), label="<ad(t)ac> 3nd steady-state", color="orange", linestyle="--")
#ax.plot(g1_solution_order_4.t, np.real(g1_solution_order_4.y[0]), label="<ad(t)ac> 4nd steady-state", color="red", linestyle="--")
#
#ax.legend()
#ax.set_xlabel('time')
#ax.set_ylabel('correlations')
#
#plt.show()



### g2



x_init_ada_order_4 = corr_get_init_vec_ada_order_4(init_state)

print("Solver started...", end='', flush=True)
solution_order_4 = solve_ivp(corr_system_ada_order_4, [0.0, 10.0], x_init_ada_order_4, first_step=0.1, max_step=0.1)
print("Done")

mean_ac_order_4 = solution_order_4.y[38][-1]
mean_adc_order_4 = np.conj(solution_order_4.y[0][-1])

def g2_sys_order_4(t, x):
    return g2_system_adcadaac_order_4(t, x, mean_ac_order_4, mean_adc_order_4)

x_init_g2_order_4 = g2_get_init_vec_adcadaac_order_4(solution_order_4.y, -1)

print("Solver started...", end='', flush=True)
g2_solution_order_4 = solve_ivp(g2_sys_order_4, [0.0, 10.0], x_init_g2_order_4, first_step=0.1, max_step=0.1)
print("Done")

fig, ax = pyplot.subplots(1, 1)
fig.set_size_inches(18.5, 10.5)
ax.plot(g2_solution_order_4.t, np.real(g2_solution_order_4.y[0]), label="<ad(t)ac> 4nd steady-state", color="red", linestyle="--")

ax.legend()
ax.set_xlabel('time')
ax.set_ylabel('correlations')

plt.show()