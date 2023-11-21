from equation_generator import *
from equation_converter import *

op = OP_SZ
#op = OP_Ad * OP_Ad * OP_A * OP_A
my_op_equ_list = develop_all_equations(op, CAVITY_H_OP_LIST, [(KAPPA, OP_A), (GAMMA, OP_SP), (NU, OP_SM)], 2)
comp_corr_equ_list = transform_equ_set_to_corr(my_op_equ_list)
print_equ_list(comp_corr_equ_list)
apply_cumulant_expansion(comp_corr_equ_list, 3)
#print("#######")
#print_equ_list(comp_corr_equ_list)
for i in range(10):
    print(i)
    new_list = complete_equations_one_pass(comp_corr_equ_list, 2)
    #print("#######")
    #print_equ_list(new_list)
    #exit()
    if len(new_list) == 0:
        break
    else:
        comp_corr_equ_list = new_list

print(len(comp_corr_equ_list))
comp_corr_equ_list = remove_dag_corr_equ(comp_corr_equ_list)
print_equ_list(comp_corr_equ_list)
print(len(comp_corr_equ_list))

print(convert_corr_equ_list_to_python_function(comp_corr_equ_list))