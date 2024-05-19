from equation_generator import *
from cavity_global import *

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
    def get_value(self, key: Cumulant, is_rust: bool):
        for pair in self.data:
            if cumu_same_qop(key, pair[0]):
                return pair[1]
        for pair in self.data:
            if cumu_dag_qop(key, pair[0]):
                if is_rust:
                    return pair[1] + ".conj()"
                else:
                    return "np.conj(" + pair[1] + ")"
        raise Exception("No key found for ", str(key))
    def __str__(self):
        res = "CumulantDict = { "
        for pair in self.data:
            res += "[" + str(pair[0]) + ", " + pair[1] + "] "
        return res + "}"
    
def get_qop_op_equ_str(op: QOp):
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
        raise Exception("Operator not defined for op_equ_str operation")

def get_qop_op_equ(op: QOp):
    if op == QOp.A:
        return op_a
    elif op == QOp.Ad:
        return op_ad
    elif op == QOp.SP:
        return op_sp
    elif op == QOp.SM:
        return op_sm
    elif op == QOp.SZ:
        return op_sz
    elif op == QOp.Ac:
        return op_ac
    elif op == QOp.Adc:
        return op_adc
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
        res += get_qop_op_equ_str(qop_dag(qop_a_list[i])) if cumu.is_dag else get_qop_op_equ_str(qop_a_list[i])
        if i != len(qop_a_list) - 1:
            res += " * "
    if len(qop_a_list) != 0 and len(qop_c_list) != 0:
        res += " * "
    for i in range(len(qop_c_list)):
        res += get_qop_op_equ_str(qop_dag(qop_c_list[i])) if cumu.is_dag else get_qop_op_equ_str(qop_c_list[i])
        if i != len(qop_c_list) - 1:
            res += " * "
    return res

def get_op_init_from_cumu(cumu: Cumulant):
    res = "["
    qop_a_list = cumu.qop_atomic_list.copy()
    qop_c_list = cumu.qop_cavity_list.copy()
    if cumu.is_dag:
        qop_a_list.reverse()
        qop_c_list.reverse()
    for i in range(len(qop_a_list)):
        res += "QOpInit."
        res += get_qop_op_equ_str(qop_dag(qop_a_list[i])) if cumu.is_dag else get_qop_op_equ_str(qop_a_list[i])
        if i != len(qop_a_list) - 1:
            res += ", "
    if len(qop_a_list) != 0 and len(qop_c_list) != 0:
        res += ", "
    for i in range(len(qop_c_list)):
        res += "QOpInit."
        res += get_qop_op_equ_str(qop_dag(qop_c_list[i])) if cumu.is_dag else get_qop_op_equ_str(qop_c_list[i])
        if i != len(qop_c_list) - 1:
            res += ", "
    return res + "]"

def get_op_from_cumu(cumu: Cumulant):
    if len(cumu.qop_cavity_list) == 0 and len(cumu.qop_atomic_list) == 0:
        return 0
    if len(cumu.qop_cavity_list) == 0:
        op = get_qop_op_equ(cumu.qop_atomic_list[0])
        for qop in cumu.qop_atomic_list[1:]:
            op *= get_qop_op_equ(qop)
        return op
    elif len(cumu.qop_atomic_list) == 0:
        op = get_qop_op_equ(cumu.qop_cavity_list[0])
        for qop in cumu.qop_cavity_list[1:]:
            op *= get_qop_op_equ(qop)
        return op
    op = get_qop_op_equ(cumu.qop_atomic_list[0]) * get_qop_op_equ(cumu.qop_cavity_list[0])
    for qop in cumu.qop_atomic_list[1:]:
        op *= get_qop_op_equ(qop)
    for qop in cumu.qop_cavity_list[1:]:
        op *= get_qop_op_equ(qop)
    return op

def get_cumu_dict(corr_equ_list, is_rust: bool):
    c_vector_dict = CumulantDict()
    if is_rust:
        for i in range(len(corr_equ_list)):
            c_vector_dict.add_key_value(corr_equ_list[i][0].cumulant_list[0], ("yn[" + str(i) + "]"))
    else:
        for i in range(len(corr_equ_list)):
            c_vector_dict.add_key_value(corr_equ_list[i][0].cumulant_list[0], ("x[" + str(i) + "]"))
    return c_vector_dict

def convert_corr_equ_list_to_python_function(corr_equ_list, order, op: Cumulant, sys_name, mean_included: bool):
    c_vector_dict = get_cumu_dict(corr_equ_list, False)
    if mean_included:
        c_vector_dict.add_key_value(Cumulant([], [QOp.Ac]), "mean_ac")
        c_vector_dict.add_key_value(Cumulant([], [QOp.Adc]), "mean_adc")
        c_vector_dict.add_key_value(Cumulant([], [QOp.Adc, QOp.Ac]), "mean_adc_ac")
        func_param = "(t, x, mean_ac, mean_adc, mean_adc_ac)"
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
                py_func += "*" + c_vector_dict.get_value(cumu, False)
        if i == len(corr_equ_list) - 1:
            py_func += " ]"
            break
        py_func += ", \n"

    return py_func

def convert_corr_equ_list_to_rust_function(corr_equ_list, order, op: Cumulant, sys_name, mean_included: bool):
    c_vector_dict = get_cumu_dict(corr_equ_list, True)
    if mean_included:
        c_vector_dict.add_key_value(Cumulant([], [QOp.Ac]), "mean_ac")
        c_vector_dict.add_key_value(Cumulant([], [QOp.Adc]), "mean_adc")
        c_vector_dict.add_key_value(Cumulant([], [QOp.Adc, QOp.Ac]), "mean_adc_ac")
        func_param = "(t, x, mean_ac, mean_adc, mean_adc_ac)"
    else:
        func_param = "(time: f64, yn: &Cvecf, fty: &mut Cvecf)"

    #print(c_vector_dict)

    py_func = "fn " + str(sys_name) + "_system_" + op.get_simple_str() + "_order_" + str(order) + func_param + "{\n"
    for i in range(len(corr_equ_list)):
        py_func += "\t fty[" + str(i) + "] = "
        for j in range(len(corr_equ_list[i][1])):
            if j != 0:
                py_func += " + "
            curr_opcumu = corr_equ_list[i][1][j]
            py_func += "(Complex::new(" + str(curr_opcumu.c_factor.real) + ", " + str(curr_opcumu.c_factor.imag) + "))"
            for factor in curr_opcumu.factor_list:
                py_func += "*" + factor.python_name
            for cumu in curr_opcumu.cumulant_list:
                py_func += "*" + c_vector_dict.get_value(cumu, True)
        if i == len(corr_equ_list) - 1: 
            py_func += "; \n}"
            break
        py_func += "; \n"

    return py_func

def get_corr_python_vec_init(corr_equ_list, order, op, sys_name):
    py_vec_init = "def " + str(sys_name) + "_get_init_vec_" + op.get_simple_str() + "_order_" + str(order) + "(init_state):\n\treturn [\n"
    for i in range(len(corr_equ_list)):
        py_vec_init += "(1.0 + 0.0j) * mean_value(init_state, " + get_op_mul_from_cumu(corr_equ_list[i][0].cumulant_list[0]) + "), # " + str(i) + "\n"
    py_vec_init += "]"
    return py_vec_init

def get_corr_python_vec_init_no_matrix(corr_equ_list, order, op, sys_name):
    py_vec_init = "def " + str(sys_name) + "_get_init_vec_" + op.get_simple_str() + "_order_" + str(order) + "(state):\n\treturn [\n"
    for i in range(len(corr_equ_list)):
        py_vec_init += "(1.0 + 0.0j) * get_projection(state, " + get_op_init_from_cumu(corr_equ_list[i][0].cumulant_list[0]) + "), # " + str(i) + "\n"
    py_vec_init += "]"
    return py_vec_init

def get_corr_rust_vec_init(corr_equ_list, order, op, sys_name):
    filename = str(sys_name) + "_init_vec_" + op.get_simple_str() + "_order_" + str(order) + "__" + str(ATOM_COUNT) + "e_" + str(PHOTON_CAPACITY) # "e_" stands for "excited"
    py_vec_init = "fn " + filename + "(init_vec: &mut Cvecf) {\n"
    for i in range(len(corr_equ_list)):
        mval = (1.0 + 0.0j) * mean_value(init_state, get_op_from_cumu(corr_equ_list[i][0].cumulant_list[0]))
        py_vec_init += "\tinit_vec[" + str(i) + "] = " + "Complex::new(" + str(mval.real) + ", " + str(mval.imag) + "); // " + str(i) + " " + cumu_to_op(corr_equ_list[i][0].cumulant_list[0]).get_simple_str() + "\n"
    py_vec_init += "}"
    return [py_vec_init, filename]

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
                    py_vec_init += "np.conj(" + other_corr_dict.get_value(other_corr_list[j][0].cumulant_list[0], False) + "[t_index]),"
                else:
                    py_vec_init += other_corr_dict.get_value(other_corr_list[j][0].cumulant_list[0], False) + "[t_index],"
                py_vec_init += "# " + cumu_to_op(corr_equ_list[i][0].cumulant_list[0]).get_simple_str() + " # " + str(i) + "\n"
                break
        if not(found):
            print(corr_equ_list[i][0].cumulant_list[0])
            raise Exception("Equivalent corr not found")
    return py_vec_init + "]"