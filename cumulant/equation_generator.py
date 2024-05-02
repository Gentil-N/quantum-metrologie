########################################
#           PARTITION GENERATOR        #
########################################



from collections import defaultdict

class Partition:

    def __init__(self, S):
        self.data = list(S)
        self.m = len(S)
        self.table = self.rgf_table()

    def __getitem__(self, i):
        #generates set partitions by index
        if i > len(self) - 1:
             raise IndexError
        L =  self.unrank_rgf(i)
        result = self.as_set_partition(L)
        return result
    
    def __len__(self):
        return self.table[self.m,0]

    def as_set_partition(self, L):
        # Transform a restricted growth function into a partition
        n = int(max(L[1:]+[1]))
        m = self.m
        data = self.data
        P = [[] for _ in range(n)]
        for i in range(m):
            P[int(L[i+1]-1)].append(data[i])
        return P

    def rgf_table(self):
        # Compute the table values 
        m = self.m
        D = defaultdict(lambda:1)
        for i in range(1,m+1):
            for j in range(0,m-i+1):
                D[i,j] = j * D[i-1,j] + D[i-1,j+1]
        return D

    def unrank_rgf(self, r):
        # Unrank a restricted growth function
        m = self.m
        L = [1 for _ in range(m+1)]
        j = 1
        D = self.table
        for i in range(2,m+1):
            v = D[m-i,j]
            cr = j*v
            if cr <= r:
                L[i] = j + 1
                r -= cr
                j += 1
            else:
                L[i] = r / v + 1
                r  %= v
        return L



########################################
#                 CODE                 #
########################################



from enum import Enum
import numpy as np
import math

class QOp(Enum):
    A = 1
    Ad = 2
    SP = 3
    SM = 4
    SZ = 5
    Ac = 6
    Adc = 7

def qop_dag(op: QOp):
    if op == QOp.A:
        return QOp.Ad
    elif op == QOp.Ad:
        return QOp.A
    elif op == QOp.SP:
        return QOp.SM
    elif op == QOp.SM:
        return QOp.SP
    elif op == QOp.SZ:
        return QOp.SZ
    elif op == QOp.Ac:
        return QOp.Adc
    elif op == QOp.Adc:
        return QOp.Ac
    else:
        raise Exception("Operator not defined for dagger operation")

def qop_str(op: QOp):
    if op == QOp.A:
        return "a"
    elif op == QOp.Ad:
        return "ad"
    elif op == QOp.SP:
        return "sp"
    elif op == QOp.SM:
        return "sm"
    elif op == QOp.SZ:
        return "sz"
    elif op == QOp.Ac:
        return "ac"
    elif op == QOp.Adc:
        return "adc"
    else:
        raise Exception("Operator not defined for str operation")

class Factor:
    def __init__(self, name="NoNameFactor", python_name="NoPythonName"):
        self.name = name
        self.python_name = python_name

class Operand:
    c_factor = 1.0 + 0.0j
    factor_list = []
    qop_atomic_list = []
    qop_cavity_list = []
    def __init__(self, c_factor=1.0+0.0j, factor_list=[], qop_a_list=[], qop_c_list=[]):
        self.c_factor = c_factor
        self.factor_list = factor_list
        self.qop_atomic_list = qop_a_list
        self.qop_cavity_list = qop_c_list
    def __mul__(self, other):
        new_c_factor = self.c_factor * other.c_factor
        new_factor_list = self.factor_list.copy()
        new_factor_list.extend(other.factor_list)
        new_qop_atomic_list = self.qop_atomic_list.copy()
        new_qop_atomic_list.extend(other.qop_atomic_list)
        new_qop_cavity_list = self.qop_cavity_list.copy()
        new_qop_cavity_list.extend(other.qop_cavity_list)
        return Operand(new_c_factor, new_factor_list, new_qop_atomic_list, new_qop_cavity_list)
    def __str__(self):
        res = str(self.c_factor)
        for fc in self.factor_list:
            res += "*" + fc.name
        res += " ("
        for qop_a in self.qop_atomic_list:
            res += " " + qop_str(qop_a)
        for qop_c in self.qop_cavity_list:
            res += " " + qop_str(qop_c)
        res += " )"
        return res
    def get_simple_str(self):
        res = ""
        for qop_a in self.qop_atomic_list:
            res += qop_str(qop_a)
        for qop_c in self.qop_cavity_list:
            res += qop_str(qop_c)
        return res
    def add_factor(self, factor: Factor):
        cop = self.copy()
        cop.factor_list.append(factor)
        return cop
    def empty_factors(self):
        cop = self.copy()
        cop.factor_list = []
        return cop
    def set_c_factor(self, c_factor):
        cop = self.copy()
        cop.c_factor = c_factor
        return cop
    def mul_c_factor(self, c_factor):
        cop = self.copy()
        cop.c_factor *= c_factor
        return cop
    def add_c_factor(self, other):
        cop = self.copy()
        cop.c_factor += other.c_factor
        return cop
    def is_null(self):
        return self.c_factor == (0.0 + 0.0j)
    def copy(self):
        return Operand(self.c_factor, self.factor_list.copy(), self.qop_atomic_list.copy(), self.qop_cavity_list.copy())
    def dag(self):
        cop = self.copy()
        cop.c_factor = np.conj(cop.c_factor)
        cop.qop_atomic_list.reverse()
        for i in range(len(cop.qop_atomic_list)):
            cop.qop_atomic_list[i] = qop_dag(cop.qop_atomic_list[i])
        cop.qop_cavity_list.reverse()
        for i in range(len(cop.qop_cavity_list)):
            cop.qop_cavity_list[i] = qop_dag(cop.qop_cavity_list[i])
        return cop
    def get_order(self):
        return len(self.qop_atomic_list) + len(self.qop_cavity_list)


def is_same_list(lista, listb):
    if len(lista) != len(listb):
        return False
    for i in range(len(lista)):
        if lista[i] != listb[i]:
            return False
    return True

def op_same_cavity_qop(opa: Operand, opb: Operand):
    return is_same_list(opa.qop_cavity_list, opb.qop_cavity_list)

def op_same_atomic_qop(opa: Operand, opb: Operand):
    return is_same_list(opa.qop_atomic_list, opb.qop_atomic_list)

def op_same_qop(opa: Operand, opb: Operand):
    return op_same_cavity_qop(opa, opb) and op_same_atomic_qop(opa, opb)

def op_list_contains_op_same_qop(op_test: Operand, op_list):
    for op in op_list:
        if op_same_qop(op_test, op):
            return True
    return False

def op_same_factors(opa: Operand, opb: Operand):
    return is_same_list(opa.factor_list, opb.factor_list)

def op_create_cumulant(op: Operand):
    return


G = Factor("G", "g")
OMEGA = Factor("ω", "omega")
DELTA = Factor("Δ", "delta")
KAPPA = Factor("K", "kappa")
GAMMA = Factor("γ", "gamma")
NU = Factor("ν", "nu")
HBAR = Factor("h", "hbar")

OP_A = Operand(1.0, [], [], [QOp.A])
OP_Ad = Operand(1.0, [], [], [QOp.Ad])
OP_N = OP_Ad * OP_A
OP_SP = Operand(1.0, [], [QOp.SP], [])
OP_SM = Operand(1.0, [], [QOp.SM], [])
OP_SZ = Operand(1.0, [], [QOp.SZ], [])
OP_Ac = Operand(1.0, [], [], [QOp.Ac])
OP_Adc = Operand(1.0, [], [], [QOp.Adc])
OP_GSPA = OP_SP * OP_A
OP_GSMAd = OP_SM * OP_Ad
CAVITY_H_OP_LIST = [OP_N.add_factor(DELTA), OP_GSPA.add_factor(G), OP_GSMAd.add_factor(G)]
HOSCILLATOR_H_OP_LIST = [OP_N.add_factor(OMEGA)]

# /!\ operators from the hamiltonian must not have hbar in it, as the function divide the whole by hbar!
def hamiltonian_commutator(hamiltonian_op_list: list, op: Operand):
    commutator_op_list = []
    for h_op in hamiltonian_op_list:
        commutator_op_list.append((h_op * op).mul_c_factor(1.0j))
        commutator_op_list.append((op * h_op).mul_c_factor(-1.0j))
    #op_n_delta_i = OP_N.add_factor(DELTA).mul_c_factor(1.0j)
    #return [op_n_delta_i * op, op * op_n_delta_i.mul_c_factor(-1.0), (OP_GSPA * op).mul_c_factor(1.0j), (OP_GSMAd * op).mul_c_factor(1.0j), (op * OP_GSPA).mul_c_factor(-1.0j), (op * OP_GSMAd).mul_c_factor(-1.0j)]
    return commutator_op_list
    
def compute_lb_terms(lb_factor: Factor, lb_op: Operand, op: Operand):
    lb_op_dag = lb_op.dag()
    return [(lb_op_dag * op * lb_op).add_factor(lb_factor), (lb_op_dag * lb_op * op).add_factor(lb_factor).mul_c_factor(-0.5), (op * lb_op_dag * lb_op).add_factor(lb_factor).mul_c_factor(-0.5)]

def add_equivalent(op_list):
    i = 0
    while i != len(op_list):
        j = i + 1
        while j != len(op_list):
            if op_same_qop(op_list[i], op_list[j]) and op_same_factors(op_list[i], op_list[j]):
                op_list[i] = op_list[i].add_c_factor(op_list[j])
                del op_list[j]
                continue
            j += 1
        i += 1

def remove_null(op_list):
    i = 0
    while i != len(op_list):
        if op_list[i].is_null():
            del op_list[i]
            continue
        i += 1

def order_cavity_qop(op: Operand):
    new_op_list = []
    i = 0
    while i < len(op.qop_cavity_list) - 1:
        if op.qop_cavity_list[i] == QOp.A and op.qop_cavity_list[i + 1] == QOp.Ad:
            op.qop_cavity_list[i] = QOp.Ad
            op.qop_cavity_list[i + 1] = QOp.A
            new_op = op.copy()
            del new_op.qop_cavity_list[i]
            del new_op.qop_cavity_list[i]
            new_op_list.append(new_op)
            new_op_list.extend(order_cavity_qop(new_op))
            i = 0
        # Normally we don't need to compare with Ac/Adc (do you have the ref ? Ok one "A" in addition...) because we never derive equations with them 
        #elif op.qop_cavity_list[i] == QOp.Ac and op.qop_cavity_list[i + 1] == QOp.Ad:
        #    op.qop_cavity_list[i] = QOp.Ad
        #    op.qop_cavity_list[i + 1] = QOp.Ac
        #    new_op = op.copy()
        #    del new_op.qop_cavity_list[i]
        #    del new_op.qop_cavity_list[i]
        #    new_op_list.append(new_op)
        #    new_op_list.extend(order_cavity_qop(new_op))
        #    i = 0
        #elif op.qop_cavity_list[i] == QOp.A and op.qop_cavity_list[i + 1] == QOp.Adc:
        #    op.qop_cavity_list[i] = QOp.Adc
        #    op.qop_cavity_list[i + 1] = QOp.A
        #    new_op = op.copy()
        #    del new_op.qop_cavity_list[i]
        #    del new_op.qop_cavity_list[i]
        #    new_op_list.append(new_op)
        #    new_op_list.extend(order_cavity_qop(new_op))
        #    i = 0
        #elif op.qop_cavity_list[i] == QOp.Ac and op.qop_cavity_list[i + 1] == QOp.Adc:
        #    op.qop_cavity_list[i] = QOp.Adc
        #    op.qop_cavity_list[i + 1] = QOp.Ac
        #    new_op = op.copy()
        #    del new_op.qop_cavity_list[i]
        #    del new_op.qop_cavity_list[i]
        #    new_op_list.append(new_op)
        #    new_op_list.extend(order_cavity_qop(new_op))
        #    i = 0
        else:
            i += 1
    return new_op_list

def order_atomic_qop(op: Operand):
    new_op_list = []
    i = 0
    while i < len(op.qop_atomic_list) - 1:
        if op.qop_atomic_list[i] == QOp.SP and op.qop_atomic_list[i + 1] == QOp.SM:
            op.qop_atomic_list[i] = QOp.SM
            op.qop_atomic_list[i + 1] = QOp.SP
            new_op = op.copy()
            del new_op.qop_atomic_list[i]
            new_op.qop_atomic_list[i] = QOp.SZ
            new_op = new_op.mul_c_factor(2.0).add_factor(HBAR)
            new_op_list.append(new_op)
            new_op_list.extend(order_atomic_qop(new_op))
            i = 0
            continue
        elif op.qop_atomic_list[i] == QOp.SM and op.qop_atomic_list[i + 1] == QOp.SZ:
            op.qop_atomic_list[i] = QOp.SZ
            op.qop_atomic_list[i + 1] = QOp.SM
            new_op = op.copy()
            del new_op.qop_atomic_list[i]
            new_op.qop_atomic_list[i] = QOp.SM
            new_op = new_op.add_factor(HBAR)
            new_op_list.append(new_op)
            new_op_list.extend(order_atomic_qop(new_op))
            i = 0
            continue
        elif op.qop_atomic_list[i] == QOp.SP and op.qop_atomic_list[i + 1] == QOp.SZ:
            op.qop_atomic_list[i] = QOp.SZ
            op.qop_atomic_list[i + 1] = QOp.SP
            new_op = op.copy()
            del new_op.qop_atomic_list[i]
            new_op.qop_atomic_list[i] = QOp.SP
            new_op = new_op.mul_c_factor(-1.0).add_factor(HBAR)
            new_op_list.append(new_op)
            new_op_list.extend(order_atomic_qop(new_op))
            i = 0
            continue
        else:
            i += 1
    return new_op_list

def order_cavity_qop_list(op_list):
    for op in op_list:
        new_op_list = order_cavity_qop(op)
        op_list.extend(new_op_list)

def order_atomic_qop_list(op_list):
    for op in op_list:
        new_op_list = order_atomic_qop(op)
        op_list.extend(new_op_list)

def develop_equation(initial_op: Operand, hamiltonian, lb_factor_tuples):
    all_op = hamiltonian_commutator(hamiltonian, initial_op)
    for lb_tuple in lb_factor_tuples:
        all_op.extend(compute_lb_terms(lb_tuple[0], lb_tuple[1], initial_op))
    add_equivalent(all_op)
    remove_null(all_op)
    order_cavity_qop_list(all_op)
    add_equivalent(all_op)
    remove_null(all_op)
    order_atomic_qop_list(all_op)
    add_equivalent(all_op)
    remove_null(all_op)
    return all_op

def left_equ_list_contains_op(equ_list, op): # Test ifleft side of equation list contain the op
    for equ in equ_list:
        if op_same_qop(equ[0], op):
            return True
    return False

def develop_all_equations(initial_op: Operand, hamiltonian, lb_factor_tuples, max_order, op_already_developed_list):
    if initial_op.get_order() > max_order:
        print(initial_op)
        raise Exception("Initial operator's order is higher than the maximal order: doesn't make sense!")
    op_to_develop_list = [initial_op.copy()]
    equ_list = []
    while len(op_to_develop_list) > 0:
        op_to_dev = op_to_develop_list[0]
        #print(op_to_dev)
        op_already_developed_list.append(op_to_dev.copy())
        del op_to_develop_list[0]

        count = 0
        right_ac_present = False
        left_adc_present = False
        if len(op_to_dev.qop_cavity_list) > 0: # We test if Adc and Ac are present only if we HAVE some cavity operators
            #print("###", cumu, " ", end='')
            if op_to_dev.qop_cavity_list[-1] == QOp.Ac:
                #print("ac ", end='')
                right_ac_present = True
                count += 1
            if op_to_dev.qop_cavity_list[0] == QOp.Adc:
                #print("adc ", end='')
                left_adc_present = True
                count += 1
            #print(len(cumu.qop_cavity_list) - count)
            if (len(op_to_dev.qop_cavity_list) - count == 0) and (len(op_to_dev.qop_atomic_list) == 0): # in that case, we don't need to test if we have the equation already developed, as the operator in question is only < ac >, < adc >, or < adc ac >: then we don't add it as an operator that we want to develop
                del op_to_develop_list[0]
                continue
        if left_adc_present: # We have a left Adc
            op_to_dev.qop_cavity_list = op_to_dev.qop_cavity_list[1:] # So we remove it!
        if right_ac_present: # This time we have a right Ac
            op_to_dev.qop_cavity_list = op_to_dev.qop_cavity_list[:-1] # We also remove it...

        #print("current = ", str(op_to_dev))
        equ_list.append([op_to_dev, develop_equation(op_to_dev, hamiltonian, lb_factor_tuples)])

        #print("####")
        #print_equ(equ_list[-1])
        if left_adc_present or right_ac_present: # if we had Adc left or Ac right
            multiply_equ_by_op(equ_list[-1], left_adc_present, right_ac_present) # we put it back!
        #print_equ(equ_list[-1])
        
        #print_op_list(equ_list[-1][1])
        for op in equ_list[-1][1]:
            if op.get_order() > max_order:
                continue
            if not(op_list_contains_op_same_qop(op, op_already_developed_list)) and not(op_list_contains_op_same_qop(op, op_to_develop_list)):
                #print_op_list(op_developed_list)
                #print("adding: ", op)
                op_to_develop_list.append(op.empty_factors().set_c_factor(1.0))
        #print(len(op_to_develop_list))
        #print("###")
    return equ_list

def multiply_equ_by_op(equ, left_adc: bool, right_ac: bool):
    if left_adc:
        equ[0] = OP_Adc * equ[0]
        for j in range(len(equ[1])):
            equ[1][j] = OP_Adc * equ[1][j]
    if right_ac:
        equ[0] = equ[0] * OP_Ac
        for j in range(len(equ[1])):
            equ[1][j] = equ[1][j] * OP_Ac
    


def print_op_list(op_list):
    for i in range(len(op_list) - 1):
        print(op_list[i], " + ", end='')
    print(op_list[-1])

def print_equ(equ):
    print(equ[0], " = ", end='')
    for i in range(len(equ[1]) - 1):
        print(equ[1][i], " + ", end='')
    print(equ[1][-1])



# <...>
class Cumulant:
    qop_atomic_list = []
    qop_cavity_list = []
    is_dag = False
    def __init__(self, qop_a_list=[], qop_c_list=[]):
        self.qop_atomic_list = qop_a_list
        self.qop_cavity_list = qop_c_list
        self.is_dag = False
    def get_order(self):
        return len(self.qop_atomic_list) + len(self.qop_cavity_list)
    def copy(self):
        return Cumulant(self.qop_atomic_list.copy(), self.qop_cavity_list.copy())
    def __str__(self):
        res = str(" <")
        a_list = reversed(self.qop_atomic_list) if self.is_dag else self.qop_atomic_list
        b_list = reversed(self.qop_cavity_list) if self.is_dag else self.qop_cavity_list
        for qop_a in a_list:
            res += " " + qop_str(qop_dag(qop_a) if self.is_dag else qop_a)
        for qop_c in b_list:
            res += " " + qop_str(qop_dag(qop_c) if self.is_dag else qop_c)
        res += " >"
        if self.is_dag:
            res += "†"
        return res
    def dag_constant(self):
        cop = self.copy()
        cop.is_dag = not(cop.is_dag)
        return cop



# <...> ... <...>
class OpCumulant:
    c_factor = 1.0 + 0.0j
    factor_list = []
    cumulant_list = []
    def __init__(self, c_factor=1.0+0.0j, factor_list=[], cumulant_list=[]):
        self.c_factor = c_factor
        self.factor_list = factor_list
        self.cumulant_list = cumulant_list
    def __str__(self):
        res = str(self.c_factor)
        for fc in self.factor_list:
            res += "*" + fc.name
        for cumu in self.cumulant_list:
            res += str(cumu)
        return res
    def copy(self):
        new_cumulant_list = []
        for cumu in self.cumulant_list:
            new_cumulant_list.append(cumu.copy())
        return OpCumulant(self.c_factor, self.factor_list.copy(), new_cumulant_list)



def cumulant_expansion(opcumu: OpCumulant, order):
    new_opcumu_list = [opcumu.copy()]
    k = -1
    while k < len(new_opcumu_list) - 1:
        k += 1
        #print(k, " ", len(new_opcumu_list))
        for l in range(len(new_opcumu_list[k].cumulant_list)):
            if new_opcumu_list[k].cumulant_list[l].get_order() >= order:
                S = set(range(0, new_opcumu_list[k].cumulant_list[l].get_order()))
                P = Partition(S)
                new_internal_opcumu_list = []
                for i in range(1, len(P)):
                    new_opcumu = OpCumulant(new_opcumu_list[k].c_factor * math.factorial(len(P[i]) - 1) * (-1)**len(P[i]), new_opcumu_list[k].factor_list.copy(), [])
                    for j in range(len(P[i])):
                        cumu = Cumulant([], [])
                        for m in range(len(P[i][j])):
                            len_atomic_list = len(new_opcumu_list[k].cumulant_list[l].qop_atomic_list)
                            #print(P[i][j][k])
                            if P[i][j][m] >= len_atomic_list:
                                cumu.qop_cavity_list.append(new_opcumu_list[k].cumulant_list[l].qop_cavity_list[P[i][j][m] - len_atomic_list])
                            else:
                                cumu.qop_atomic_list.append(new_opcumu_list[k].cumulant_list[l].qop_atomic_list[P[i][j][m]])
                        new_opcumu.cumulant_list.append(cumu.copy())
                    #print(new_opcumu)
                    for n in range(l):
                        new_opcumu.cumulant_list.append(new_opcumu_list[k].cumulant_list[n].copy())
                    for n in range(l+1, len(new_opcumu_list[k].cumulant_list)):
                        new_opcumu.cumulant_list.append(new_opcumu_list[k].cumulant_list[n].copy())
                    new_internal_opcumu_list.append(new_opcumu.copy())
                del new_opcumu_list[k]
                new_opcumu_list.extend(new_internal_opcumu_list)
                k = -1
                break
    return new_opcumu_list

    #if len(opcumu.cumulant_list) > 1 or opcumu.cumulant_list[0].get_order() != order:
    #    return [opcumu]
    #S = set(range(0, order))
    #P = Partition(S)
    ##for x in P:
    ##    print(x)
    #new_opcumu_list = []
    #for i in range(1, len(P)):
    #    new_opcumu = OpCumulant(opcumu.c_factor * math.factorial(len(P[i]) - 1) * (-1)**len(P[i]), opcumu.factor_list.copy(), [])
    #    for j in range(len(P[i])):
    #        cumu = Cumulant([], [])
    #        for k in range(len(P[i][j])):
    #            len_atomic_list = len(opcumu.cumulant_list[0].qop_atomic_list)
    #            #print(P[i][j][k])
    #            if P[i][j][k] >= len_atomic_list:
    #                cumu.qop_cavity_list.append(opcumu.cumulant_list[0].qop_cavity_list[P[i][j][k] - len_atomic_list])
    #            else:
    #                cumu.qop_atomic_list.append(opcumu.cumulant_list[0].qop_atomic_list[P[i][j][k]])
    #        new_opcumu.cumulant_list.append(cumu.copy())
    #    #print(new_opcumu)
    #    new_opcumu_list.append(new_opcumu.copy())
    #return new_opcumu_list



def op_to_corr(op: Operand):
    return OpCumulant(op.c_factor, op.factor_list.copy(), [Cumulant(op.qop_atomic_list.copy(), op.qop_cavity_list.copy())])

def cumu_to_op(cumu: Cumulant):
    return Operand(1.0, [], cumu.qop_atomic_list, cumu.qop_cavity_list)

def transform_to_correlation(equ):
    corr_equ = [op_to_corr(equ[0])]
    op_cumu_list = []
    for op in equ[1]:
        op_cumu_list.append(op_to_corr(op))
    corr_equ.append(op_cumu_list)
    return corr_equ

def transform_equ_set_to_corr(equ_list):
    corr_list = []
    for equ in equ_list:
        corr_list.append(transform_to_correlation(equ))
    return corr_list

def apply_cumulant_expansion(corr_equ_list, order):
    for i in range(len(corr_equ_list)):
        new_corr = []
        for opcumu in corr_equ_list[i][1]:
            new_corr.extend(cumulant_expansion(opcumu, order))
        corr_equ_list[i][1] = new_corr

def op_list_contains_cumu(op_list, cumu):
    for op in op_list:
        if is_same_list(op.qop_atomic_list, cumu.qop_atomic_list) and is_same_list(op.qop_cavity_list, cumu.qop_cavity_list):
            return True
    return False

def op_list_contains_cumu_special_adc_ac(op_list, cumu):
    for op in op_list:
        if is_same_list(op[0].qop_atomic_list, cumu.qop_atomic_list) and is_same_list(op[0].qop_cavity_list, cumu.qop_cavity_list):
            return True
    return False

# /!\ Must be applied after the cumulant expansion! Otherwise, the other functions will run for a long, long, very long time... (infinite) or your computer will crash before ;p
# /!\ ASSUMING THAT WE ONLY HAVE AC ON THE RIGHT OR ADC ON THE LEFT
def find_complete_op_corr_equ(corr_equ_list):
    op_list_to_develop = []
    for equ in corr_equ_list:
        for opcumu in equ[1]:
            for cumu in opcumu.cumulant_list:
                count = 0
                #right_ac_present = False
                #left_adc_present = False
                if len(cumu.qop_cavity_list) > 0: # We test if Adc and Ac are present only if we HAVE some cavity operators
                    #print("###", cumu, " ", end='')
                    if cumu.qop_cavity_list[-1] == QOp.Ac:
                        #print("ac ", end='')
                        #right_ac_present = True
                        count += 1
                    if cumu.qop_cavity_list[0] == QOp.Adc:
                        #print("adc ", end='')
                        #left_adc_present = True
                        count += 1
                    #print(len(cumu.qop_cavity_list) - count)
                    if (len(cumu.qop_cavity_list) - count == 0) and (len(cumu.qop_atomic_list) == 0): # in that case, we don't need to test if we have the equation already developed, as the operator in question is only < ac >, < adc >, or < adc ac >: then we don't add it as an operator that we wan't to develop
                        continue
                found = False
                for equ in corr_equ_list:
                    if cumu_same_qop(equ[0].cumulant_list[0], cumu): #is_same_list(equ[0].cumulant_list[0].qop_atomic_list, cumu.qop_atomic_list) and is_same_list(equ[0].cumulant_list[0].qop_cavity_list, cumu.qop_cavity_list):
                        found = True
                if not(found) and not(op_list_contains_cumu(op_list_to_develop, cumu)):
                    op_list_to_develop.append(Operand(1.0, [], cumu.qop_atomic_list, cumu.qop_cavity_list))
    #for op in op_list_to_develop:
    #    print(str(op))
    return op_list_to_develop

    #op_list_to_develop = []
    #for equ in corr_equ_list:
    #    for opcumu in equ[1]:
    #        for cumu in opcumu.cumulant_list:
    #            found = False
    #            for equ in corr_equ_list:
    #                if is_same_list(equ[0].cumulant_list[0].qop_atomic_list, cumu.qop_atomic_list) and is_same_list(equ[0].cumulant_list[0].qop_cavity_list, cumu.qop_cavity_list):
    #                    found = True
    #            if not(found) and not(op_list_contains_cumu(op_list_to_develop, cumu)):
    #                op_list_to_develop.append(Operand(1.0, [], cumu.qop_atomic_list, cumu.qop_cavity_list))
    ##for op in op_list_to_develop:
    ##    print(str(op))
    #return op_list_to_develop

def cumu_same_qop(cumua: Cumulant, cumub: Cumulant):
    if is_same_list(cumua.qop_atomic_list, cumub.qop_atomic_list) and is_same_list(cumua.qop_cavity_list, cumub.qop_cavity_list):
        return True
    return False

def remove_same_corr_equ(corr_equ_list_a, corr_equ_list_b):
    final_set_corr_equ = corr_equ_list_b.copy()
    for i in range(len(corr_equ_list_a)):
        found = False
        for j in range(len(corr_equ_list_b)):
            if cumu_same_qop(corr_equ_list_a[i][0].cumulant_list[0], corr_equ_list_b[j][0].cumulant_list[0]):
                found = True
        if not(found):
            final_set_corr_equ.append(corr_equ_list_a[i])
    return final_set_corr_equ

def cumu_dag_qop(cumua: Cumulant, cumub: Cumulant):
    if len(cumua.qop_atomic_list) != len(cumub.qop_atomic_list) or len(cumua.qop_cavity_list) != len(cumub.qop_cavity_list):
        return False
    for i in range(len(cumua.qop_cavity_list)):
        if cumua.qop_cavity_list[i] != qop_dag(cumub.qop_cavity_list[len(cumua.qop_cavity_list) - i - 1]):
            return False
    for i in range(len(cumua.qop_atomic_list)):
        if cumua.qop_atomic_list[i] != qop_dag(cumub.qop_atomic_list[len(cumua.qop_atomic_list) - i - 1]):
            return False
    return True

def remove_dag_corr_equ(corr_equ_list, original_op: Cumulant):
    new_corr_equ_list = []
    cumu_to_set_dag_list = []
    i = 0
    while i < len(corr_equ_list):
        new_corr_equ_list.append(corr_equ_list[i])
        for j in range(i + 1, len(corr_equ_list)):
            if cumu_dag_qop(corr_equ_list[i][0].cumulant_list[0], corr_equ_list[j][0].cumulant_list[0]):
                if cumu_same_qop(corr_equ_list[j][0].cumulant_list[0], original_op):
                    cumu_to_set_dag_list.append(corr_equ_list[i][0].cumulant_list[0].copy())
                    new_corr_equ_list[-1] = corr_equ_list[j]
                else:
                    cumu_to_set_dag_list.append(corr_equ_list[j][0].cumulant_list[0].copy())
                del corr_equ_list[j]
                break
        i += 1
    #replace targeted cumu by their dag
    for corr_equ in new_corr_equ_list:
        for opcumu in corr_equ[1]:
            for i in range(len(opcumu.cumulant_list)):
                for cumu_to_set in cumu_to_set_dag_list:
                    if cumu_same_qop(opcumu.cumulant_list[i], cumu_to_set):
                        opcumu.cumulant_list[i] = opcumu.cumulant_list[i].dag_constant()
                        break
    return new_corr_equ_list

#my_op = OP_N.copy()
#all_op = hamiltonian_commutator(CAVITY_H_OP_LIST, my_op)
#lb_a_terms = compute_lb_terms(KAPPA, OP_A, my_op)
#lb_sp_terms = compute_lb_terms(GAMMA, OP_SP, my_op)
#lb_sm_terms = compute_lb_terms(NU, OP_SM, my_op)
#all_op.extend(lb_a_terms)
#all_op.extend(lb_sp_terms)
#all_op.extend(lb_sm_terms)
#
#print_op_list(all_op)
#
#add_equivalent(all_op)
#
##print_op_list(all_op)
#
#remove_null(all_op)
#
##print_op_list(all_op)
#
#order_cavity_qop_list(all_op)
#
##print_op_list(all_op)
#
#add_equivalent(all_op)
#
##print_op_list(all_op)
#
#remove_null(all_op)
#
#print_op_list(all_op)
#
#test_list = [Operand(2, [], [QOp.SM, QOp.SZ, QOp.SP], []), Operand(-1, [], [QOp.SM, QOp.SP, QOp.SZ], []), Operand(-1, [], [QOp.SZ, QOp.SM, QOp.SP], [])]
#test_list = [Operand(2, [], [QOp.SM, QOp.SM, QOp.SP], []), Operand(-1, [], [QOp.SM, QOp.SP, QOp.SM], []), Operand(-1, [], [QOp.SM, QOp.SM, QOp.SP], [])]
#test_list = [Operand(2, [], [QOp.SP, QOp.SZ, QOp.SM], []), Operand(-1, [], [QOp.SP, QOp.SM, QOp.SZ], []), Operand(-1, [], [QOp.SZ, QOp.SP, QOp.SM], [])]
#
#print_op_list(test_list)
#
#order_atomic_qop_list(test_list)
#
##print_op_list(test_list)
#
#add_equivalent(test_list)
#
##print_op_list(test_list)
#
#remove_null(test_list)
#
#print_op_list(test_list)



#my_op = OP_Ad.copy()
#all_op = hamiltonian_commutator(HOSCILLATOR_H_OP_LIST, my_op)
#lb_a_terms = compute_lb_terms(KAPPA, OP_A, my_op)
#all_op.extend(lb_a_terms)
#
#add_equivalent(all_op)
#remove_null(all_op)
#order_cavity_qop_list(all_op)
#add_equivalent(all_op)
#remove_null(all_op)
#print_op_list(all_op)



#my_op = OP_A
#all_op = hamiltonian_commutator(CAVITY_H_OP_LIST, my_op)
#lb_a_terms = compute_lb_terms(KAPPA, OP_A, my_op)
#lb_sp_terms = compute_lb_terms(GAMMA, OP_SP, my_op)
#lb_sm_terms = compute_lb_terms(NU, OP_SM, my_op)
#all_op.extend(lb_a_terms)
#all_op.extend(lb_sp_terms)
#all_op.extend(lb_sm_terms)
#
#add_equivalent(all_op)
#remove_null(all_op)
#order_cavity_qop_list(all_op)
#add_equivalent(all_op)
#remove_null(all_op)
#print_op_list(all_op)



#my_op = OP_SM
#all_op = hamiltonian_commutator(CAVITY_H_OP_LIST, my_op)
#lb_a_terms = compute_lb_terms(KAPPA, OP_A, my_op)
#lb_sp_terms = compute_lb_terms(GAMMA, OP_SP, my_op)
#lb_sm_terms = compute_lb_terms(NU, OP_SM, my_op)
#all_op.extend(lb_a_terms)
#all_op.extend(lb_sp_terms)
#all_op.extend(lb_sm_terms)
#
#add_equivalent(all_op)
#remove_null(all_op)
#order_cavity_qop_list(all_op)
#add_equivalent(all_op)
#remove_null(all_op)
#order_atomic_qop_list(all_op)
#add_equivalent(all_op)
#remove_null(all_op)
#print_op_list(all_op)



#my_op = OP_SZ
#all_op = hamiltonian_commutator(CAVITY_H_OP_LIST, my_op)
#lb_a_terms = compute_lb_terms(KAPPA, OP_A, my_op)
#lb_sp_terms = compute_lb_terms(GAMMA, OP_SP, my_op)
#lb_sm_terms = compute_lb_terms(NU, OP_SM, my_op)
#all_op.extend(lb_a_terms)
#all_op.extend(lb_sp_terms)
#all_op.extend(lb_sm_terms)
#
#add_equivalent(all_op)
#remove_null(all_op)
#order_cavity_qop_list(all_op)
#add_equivalent(all_op)
#remove_null(all_op)
#order_atomic_qop_list(all_op)
#add_equivalent(all_op)
#remove_null(all_op)
#print_op_list(all_op)



#my_op = OP_SZ * OP_SZ
#all_op = hamiltonian_commutator(CAVITY_H_OP_LIST, my_op)
#lb_a_terms = compute_lb_terms(KAPPA, OP_A, my_op)
#lb_sp_terms = compute_lb_terms(GAMMA, OP_SP, my_op)
#lb_sm_terms = compute_lb_terms(NU, OP_SM, my_op)
#all_op.extend(lb_a_terms)
#all_op.extend(lb_sp_terms)
#all_op.extend(lb_sm_terms)
#
#add_equivalent(all_op)
#remove_null(all_op)
#order_cavity_qop_list(all_op)
#add_equivalent(all_op)
#remove_null(all_op)
#order_atomic_qop_list(all_op)
#add_equivalent(all_op)
#remove_null(all_op)
#print_op_list(all_op)



#opcumu = OpCumulant(1.0, [], [Cumulant([QOp.SZ, QOp.SM, QOp.SP], [])])
#opcumu_list = cumulant_expansion(opcumu, 3)
#print_op_list(opcumu_list)


def print_equ_list(equ_list):
    for equ in equ_list:
        print_equ(equ)
        print()

def get_op_developed_list_from_corr_equ_list(corr_equ_list):
    op_list = []
    for equ in corr_equ_list:
        op_list.append(cumu_to_op(equ[0].cumulant_list[0]))
    return op_list

def complete_equations_one_pass(corr_equ_list, order):
    op_comp_list = find_complete_op_corr_equ(corr_equ_list)
    #print_op_list(op_comp_list)
    if len(op_comp_list) == 0:
        #print(len(corr_equ_list))
        return []

    #for op in op_comp_list:
    #    print("op: ", op)
    #exit()

    op_already_developed_list = get_op_developed_list_from_corr_equ_list(corr_equ_list)

    new_corr_equ_list = corr_equ_list.copy()
    for op in op_comp_list:
        if op_list_contains_op_same_qop(op, op_already_developed_list):
            continue
        #print("\n\n\ndevelopment: ", str(op))
        equ_list = develop_all_equations(op, CAVITY_H_OP_LIST, [(KAPPA, OP_A), (GAMMA, OP_SP), (NU, OP_SM)], order, op_already_developed_list)
        #print("equ")
        #print_equ_list(equ_list)
        current_corr_equ_list = transform_equ_set_to_corr(equ_list)
        apply_cumulant_expansion(current_corr_equ_list, order + 1)
        #print("cumu")
        #print_equ_list(new_corr_equ_list)
        new_corr_equ_list.extend(current_corr_equ_list)
        #print_equ_list(new_corr_equ_list)

    #print_equ_list(corr_equ_list)
    #final_corr_equ_list = remove_same_corr_equ(corr_equ_list, corr_equ_list_set[0])
    ##print_equ_list(final_corr_equ_list)
    #for i in range(len(corr_equ_list_set)):
    #    #print("\n######################################\n")
    #    #print_equ_list(corr_equ_list_set[i])
    #    final_corr_equ_list = remove_same_corr_equ(final_corr_equ_list, corr_equ_list_set[i])

    #print(len(final_corr_equ_list))
    #print_equ_list(final_corr_equ_list)
    return new_corr_equ_list

#op = OP_Ad * OP_Ad * OP_A * OP_A
#my_op_equ_list = develop_all_equations(op, CAVITY_H_OP_LIST, [(KAPPA, OP_A), (GAMMA, OP_SP), (NU, OP_SM)], 4)
#comp_corr_equ_list = transform_equ_set_to_corr(my_op_equ_list)
#apply_cumulant_expansion(comp_corr_equ_list, 5)
#for i in range(10):
#    print(i)
#    comp_corr_equ_list = complete_equations_one_pass(comp_corr_equ_list, 4)
#
#print_equ_list(len(comp_corr_equ_list))

#op = OP_Ad * OP_A
#my_op_equ_list = develop_all_equations(op, CAVITY_H_OP_LIST, [(KAPPA, OP_A), (GAMMA, OP_SP), (NU, OP_SM)], 2)
#print("\n######################################\n")
#print_equ_list(my_op_equ_list)
#print("\n######################################\n")
#comp_corr_equ_list = transform_equ_set_to_corr(my_op_equ_list)
#print_equ_list(comp_corr_equ_list)
#print("\n######################################\n")
#apply_cumulant_expansion(comp_corr_equ_list, 3)
#print_equ_list(comp_corr_equ_list)
#print("\n######################################\n")
#comp_corr_equ_list = complete_equations_one_pass(comp_corr_equ_list, 2)
#print_equ_list(comp_corr_equ_list)
#print("\n######################################\n")
#comp_corr_equ_list = complete_equations_one_pass(comp_corr_equ_list, 2)
#print_equ_list(comp_corr_equ_list)
#print(len(comp_corr_equ_list))
#print("\n######################################\n")
#comp_corr_equ_list = remove_dag_corr_equ(comp_corr_equ_list)
#print_equ_list(comp_corr_equ_list)
#print(len(comp_corr_equ_list))
#print("\n######################################\n")

#op = OP_Ad * OP_Ac
#order = 4
#my_op_equ_list = develop_all_equations(op, CAVITY_H_OP_LIST, [(KAPPA, OP_A), (GAMMA, OP_SP), (NU, OP_SM)], order, [])
#print_equ_list(my_op_equ_list)
#print("\n######################################\n")
#comp_corr_equ_list = transform_equ_set_to_corr(my_op_equ_list)
#print_equ_list(comp_corr_equ_list)
#print("\n######################################\n")
#apply_cumulant_expansion(comp_corr_equ_list, order + 1)
#print_equ_list(comp_corr_equ_list)
#print("\n######################################\n")
#for i in range(10):
#    new_list = complete_equations_one_pass(comp_corr_equ_list, order) # not "order + ..." because we process the case of Ac/Adc in the function itself
#    print(len(new_list))
#    #print_equ_list(new_list)
#    print("\n######################################\n")
#    if len(new_list) == 0:
#        break
#    else:
#        comp_corr_equ_list = new_list
#
#for equ in comp_corr_equ_list:
#    print(equ[0])

#op = OP_Adc * OP_Ad * OP_A * OP_Ac
#order = 4
#my_op_equ_list = develop_all_equations(op, CAVITY_H_OP_LIST, [(KAPPA, OP_A), (GAMMA, OP_SP), (NU, OP_SM)], order, [])
##print_equ_list(my_op_equ_list)
#print("\n######################################\n")
#comp_corr_equ_list = transform_equ_set_to_corr(my_op_equ_list)
##print_equ_list(comp_corr_equ_list)
#print("\n######################################\n")
#apply_cumulant_expansion(comp_corr_equ_list, order + 1)
##print_equ_list(comp_corr_equ_list)
#print("\n######################################\n")
#for i in range(10):
#    new_list = complete_equations_one_pass(comp_corr_equ_list, order) # not "order + ..." because we process the case of Ac/Adc in the function itself
#    print(len(new_list))
#    #print_equ_list(new_list)
#    print("\n######################################\n")
#    if len(new_list) == 0:
#        break
#    else:
#        comp_corr_equ_list = new_list
#        
#for equ in comp_corr_equ_list:
#    print(equ[0])


op = OP_Ad
order = 1
my_op_equ_list = develop_all_equations(op, CAVITY_H_OP_LIST, [(KAPPA, OP_A), (GAMMA, OP_SP), (NU, OP_SM)], order, [])
#print_equ_list(my_op_equ_list)
#print("\n######################################\n")
comp_corr_equ_list = transform_equ_set_to_corr(my_op_equ_list)
print_equ_list(comp_corr_equ_list)
#print("\n######################################\n")
apply_cumulant_expansion(comp_corr_equ_list, order + 1)
#print_equ_list(comp_corr_equ_list)
#print("\n######################################\n")
for i in range(10):
    new_list = complete_equations_one_pass(comp_corr_equ_list, order) # not "order + ..." because we process the case of Ac/Adc in the function itself
    print(len(new_list))
    #print_equ_list(new_list)
    print("\n######################################\n")
    if len(new_list) == 0:
        break
    else:
        comp_corr_equ_list = new_list

comp_corr_equ_list = remove_dag_corr_equ(comp_corr_equ_list, op_to_corr(op).cumulant_list[0])
print_equ_list(comp_corr_equ_list)