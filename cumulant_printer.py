from equation_generator import *
from equation_converter import *

op = OP_Ad * OP_A
order = 4

#op = OP_Ad * OP_Ad * OP_A * OP_A
print("Equation generation...", end="", flush=True)
my_op_equ_list = develop_all_equations(op, CAVITY_H_OP_LIST, [(KAPPA, OP_A), (GAMMA, OP_SP), (NU, OP_SM)], order)
comp_corr_equ_list = transform_equ_set_to_corr(my_op_equ_list)
#print_equ_list(comp_corr_equ_list)
apply_cumulant_expansion(comp_corr_equ_list, order + 1)
#print("#######")
#print_equ_list(comp_corr_equ_list)
for i in range(10):
    print(i, "...", end="", flush=True)
    new_list = complete_equations_one_pass(comp_corr_equ_list, order)
    #print("#######")
    #print_equ_list(new_list)
    #exit()
    if len(new_list) == 0:
        break
    else:
        comp_corr_equ_list = new_list

print("Done")

print("Simplification...", end="", flush=True)
print(len(comp_corr_equ_list), "...", end="", flush=True)
comp_corr_equ_list = remove_dag_corr_equ(comp_corr_equ_list)
#print_equ_list(comp_corr_equ_list)
print(len(comp_corr_equ_list), "...Done")

print("File creation...", end="", flush=True)
python_print = convert_corr_equ_list_to_python_function(comp_corr_equ_list, order, op)
file = open("corr_sys_" + op.get_simple_str() + "_order_" + str(order) + ".py", 'w')
file.write(python_print[0] + "\n\n" + python_print[1])
file.close()
print("Done")