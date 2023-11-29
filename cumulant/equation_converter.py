from equation_generator import *

#def cumu_dag_same_qop(cumua: Cumulant, cumub_dag: Cumulant):
#    if len(cumua.qop_atomic_list) != len(cumub_dag.qop_atomic_list) or len(cumua.qop_cavity_list) != len(cumub_dag.qop_cavity_list):
#        return False
#    for i in range(len(cumua.qop_atomic_list)):
#        print(i, " ", len(cumua.qop_atomic_list), " ", len(cumub_dag.qop_atomic_list))
#        if qop_dag(cumua.qop_atomic_list[i]) != cumub_dag.qop_atomic_list[len(cumub_dag.qop_atomic_list) - i - 1]:
#            return False
#    for i in range(len(cumua.qop_cavity_list)):
#        if qop_dag(cumua.qop_cavity_list[i]) != cumub_dag.qop_cavity_list[len(cumub_dag.qop_cavity_list) - i - 1]:
#            return False
#    return True

class CumulantDict:
    def __init__(self):
        self.data = []
    def add_key_value(self, key: Cumulant, value: str):
        self.data.append((key, value))
    def get_value(self, key: Cumulant):
        for pair in self.data:
            if cumu_same_qop(key, pair[0]):
                return pair[1]
        for pair in self.data:
            if cumu_dag_qop(key, pair[0]):
                return "np.conj(" + pair[1] + ")"
        raise Exception("No key found for ", str(key))
    def __str__(self):
        res = "CumulantDict = { "
        for pair in self.data:
            res += "[" + str(pair[0]) + ", " + pair[1] + "] "
        return res + "}"
    
def get_qop_op_equ(op: QOp):
    if op == QOp.A:
        return "op_a"
    elif op == QOp.Ad:
        return "op_ad"
    elif op == QOp.SP:
        return "op_sp"
    elif op == QOp.SM:
        return "op_sm"
    elif op == QOp.SZ:
        return "op_sz"
    elif op == QOp.Ac:
        return "op_ac"
    elif op == QOp.Adc:
        return "op_adc"
    else:
        raise Exception("Operator not defined for op_equ operation")

def get_op_mul_from_cumu(cumu: Cumulant):
    res = ""
    qop_a_list = cumu.qop_atomic_list.copy()
    qop_c_list = cumu.qop_cavity_list.copy()
    if cumu.is_dag:
        qop_a_list.reverse()
        qop_c_list.reverse()
    for i in range(len(qop_a_list)):
        res += get_qop_op_equ(qop_dag(qop_a_list[i])) if cumu.is_dag else get_qop_op_equ(qop_a_list[i])
        if i != len(qop_a_list) - 1:
            res += " * "
    if len(qop_a_list) != 0 and len(qop_c_list) != 0:
        res += " * "
    for i in range(len(qop_c_list)):
        res += get_qop_op_equ(qop_dag(qop_c_list[i])) if cumu.is_dag else get_qop_op_equ(qop_c_list[i])
        if i != len(qop_c_list) - 1:
            res += " * "
    return res

def get_cumu_dict(corr_equ_list):
    c_vector_dict = CumulantDict()
    for i in range(len(corr_equ_list)):
        c_vector_dict.add_key_value(corr_equ_list[i][0].cumulant_list[0], ("x[" + str(i) + "]"))
    return c_vector_dict

def convert_corr_equ_list_to_python_function(corr_equ_list, order, op: Cumulant, sys_name, mean_included: bool):
    c_vector_dict = get_cumu_dict(corr_equ_list)
    if mean_included:
        c_vector_dict.add_key_value(Cumulant([], [QOp.Ac]), "mean_ac")
        c_vector_dict.add_key_value(Cumulant([], [QOp.Adc]), "mean_adc")
        c_vector_dict.add_key_value(Cumulant([], [QOp.Adc, QOp.Ac]), "mean_adc_ac")
        func_param = "(t, x, mean_ac, mean_adc)"
    else:
        func_param = "(t, x)"

    #print(c_vector_dict)

    py_func = "import numpy as np\nfrom cavity_global import *\n\ndef " + str(sys_name) + "_system_" + op.get_simple_str() + "_order_" + str(order) + func_param + ":\n\treturn [ "
    for i in range(len(corr_equ_list)):
        for j in range(len(corr_equ_list[i][1])):
            if j != 0:
                py_func += " + "
            curr_opcumu = corr_equ_list[i][1][j]
            py_func += "(" + str(curr_opcumu.c_factor) + ")"
            for factor in curr_opcumu.factor_list:
                py_func += "*" + factor.python_name
            for cumu in curr_opcumu.cumulant_list:
                py_func += "*" + c_vector_dict.get_value(cumu)
        if i == len(corr_equ_list) - 1:
            py_func += " ]"
            break
        py_func += ", \n"

    return py_func

def get_corr_python_vec_init(corr_equ_list, order, op, sys_name):
    py_vec_init = "def " + str(sys_name) + "_get_init_vec_" + op.get_simple_str() + "_order_" + str(order) + "(init_state):\n\treturn [\n"
    for i in range(len(corr_equ_list)):
        py_vec_init += "(1.0 + 0.0j) * mean_value(init_state, " + get_op_mul_from_cumu(corr_equ_list[i][0].cumulant_list[0]) + "),\n"
    py_vec_init += "]"
    return py_vec_init

def is_equivalent_cumu(cumua: Cumulant, cumub: Cumulant): #return [if equivalent, if dag]
    cumua_c = cumua.copy()
    cumub_c = cumub.copy()
    if cumu_same_qop(cumua_c, cumub_c):
        return [True, False]
    if cumu_dag_qop(cumua_c, cumub_c):
        return [True, True]
    for i in range(len(cumua_c.qop_cavity_list)):
        if cumua_c.qop_cavity_list[i] == QOp.Ac:
            cumua_c.qop_cavity_list[i] = QOp.A
        elif cumua_c.qop_cavity_list[i] == QOp.Adc:
            cumua_c.qop_cavity_list[i] = QOp.Ad
    for i in range(len(cumub_c.qop_cavity_list)):
        if cumub_c.qop_cavity_list[i] == QOp.Ac:
            cumub_c.qop_cavity_list[i] = QOp.A
        elif cumub_c.qop_cavity_list[i] == QOp.Adc:
            cumub_c.qop_cavity_list[i] = QOp.Ad
    if cumu_same_qop(cumua_c, cumub_c):
        return [True, False]
    if cumu_dag_qop(cumua_c, cumub_c):
        return [True, True]
    return [False, False]

def get_corr_python_vec_init_from_other_corr(corr_equ_list, other_corr_list, other_corr_dict: CumulantDict, sys_name, order, op):
    py_vec_init = "def " + str(sys_name) + "_get_init_vec_" + op.get_simple_str() + "_order_" + str(order) + "(x, t_index):\n\treturn [\n"
    for i in range(len(corr_equ_list)):
        found = False
        for j in range(len(other_corr_list)):
            res = is_equivalent_cumu(corr_equ_list[i][0].cumulant_list[0], other_corr_list[j][0].cumulant_list[0])
            if res[0]:
                found = True
                if res[1]:
                    py_vec_init += "np.conj(" + other_corr_dict.get_value(other_corr_list[j][0].cumulant_list[0]) + "[t_index]),"
                else:
                    py_vec_init += other_corr_dict.get_value(other_corr_list[j][0].cumulant_list[0]) + "[t_index],"
                py_vec_init += "# " + cumu_to_op(corr_equ_list[i][0].cumulant_list[0]).get_simple_str() + "\n"
                break
        if not(found):
            print(corr_equ_list[i][0].cumulant_list[0])
            raise Exception("Equivalent corr not found")
    return py_vec_init + "]"