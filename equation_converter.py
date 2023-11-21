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

def convert_corr_equ_list_to_python_function(corr_equ_list):
    c_vector_dict = CumulantDict()
    for i in range(len(corr_equ_list)):
        c_vector_dict.add_key_value(corr_equ_list[i][0].cumulant_list[0], ("x[" + str(i) + "]"))
    print(c_vector_dict)
    py_func = "return [ "
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