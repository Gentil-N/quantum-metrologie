from equation_generator import *
from equation_converter import *



def get_equ_sys(op, order):
    print("Equation generation...", end="", flush=True)
    my_op_equ_list = develop_all_equations(op, CAVITY_H_OP_LIST, [(KAPPA, OP_A), (GAMMA, OP_SP), (NU, OP_SM)], order, [])
    #print_equ_list(my_op_equ_list)
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
    comp_corr_equ_list = remove_dag_corr_equ(comp_corr_equ_list, op_to_corr(op).cumulant_list[0])
    #print_equ_list(comp_corr_equ_list)
    print(len(comp_corr_equ_list), "...Done")
    return comp_corr_equ_list



def corr_cumulant_printer(op, order):

    comp_corr_equ_list = get_equ_sys(op, order)

    print("File creation...", end="", flush=True)
    python_func = convert_corr_equ_list_to_python_function(comp_corr_equ_list, order, op, "corr", False)
    python_vec = get_corr_python_vec_init(comp_corr_equ_list, order, op, "corr")
    file = open("corr_sys_" + op.get_simple_str() + "_order_" + str(order) + ".py", 'w')
    file.write(python_func + "\n\n" + python_vec)
    file.close()
    print("Done")
    #print_equ_list(comp_corr_equ_list)



def gn_cumulant_printer(op_init, op, order, g_name):
    print("Computing op_init equation system:")
    comp_init_corr_equ_list = get_equ_sys(op_init, order)
    init_dict = get_cumu_dict(comp_init_corr_equ_list)
    print("Computing op equation system:")
    comp_corr_equ_list = get_equ_sys(op, order)

    print("File creation...", end="", flush=True)
    python_func = convert_corr_equ_list_to_python_function(comp_corr_equ_list, order, op, g_name, True)
    python_vec_init = get_corr_python_vec_init_from_other_corr(comp_corr_equ_list, comp_init_corr_equ_list, init_dict, g_name, order, op)
    file = open(g_name + "_sys_" + op.get_simple_str() + "_order_" + str(order) + ".py", 'w')
    file.write(python_func + "\n\n" + python_vec_init)
    file.close()
    print("Done")

    print("File creation...", end="", flush=True)
    python_func = convert_corr_equ_list_to_python_function(comp_init_corr_equ_list, order, op_init, "corr", False)
    python_vec = get_corr_python_vec_init(comp_init_corr_equ_list, order, op_init, "corr")
    file = open("corr_sys_" + op_init.get_simple_str() + "_order_" + str(order) + ".py", 'w')
    file.write(python_func + "\n\n" + python_vec)
    file.close()
    print("Done")
    #print_equ_list(comp_corr_equ_list)



### PRINT
#corr_cumulant_printer(OP_Ad * OP_A, 2)
#corr_cumulant_printer(OP_Ad * OP_A, 3)
#corr_cumulant_printer(OP_Ad * OP_A, 4)
#gn_cumulant_printer(OP_Ad, OP_Ad * OP_Ac, 2, "g1")
#gn_cumulant_printer(OP_Ad, OP_Ad * OP_Ac, 3, "g1")
#gn_cumulant_printer(OP_Ad, OP_Ad * OP_Ac, 4, "g1")
gn_cumulant_printer(OP_Ad * OP_A, OP_Adc * OP_Ad * OP_A * OP_Ac, 5, "g2")