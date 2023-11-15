from enum import Enum
import numpy as np

class QOp(Enum):
    A = 1
    Ad = 2
    SP = 3
    SM = 4
    SZ = 5

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
    else:
        raise Exception("Operator not defined for dagger operation")

def qop_str(op: QOp):
    if op == QOp.A:
        return "a"
    elif op == QOp.Ad:
        return "a'"
    elif op == QOp.SP:
        return "sp"
    elif op == QOp.SM:
        return "sm"
    elif op == QOp.SZ:
        return "sz"
    else:
        raise Exception("Operator not defined for str operation")

class Factor:
    name = "NoNameFactor"
    def __init__(self, name="NoNameFactor"):
        self.name = name

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
        res += "("
        for qop_a in self.qop_atomic_list:
            res += " " + qop_str(qop_a)
        for qop_c in self.qop_cavity_list:
            res += " " + qop_str(qop_c)
        res += " )"
        return res
    def add_factor(self, factor: Factor):
        cop = self.copy()
        cop.factor_list.append(factor)
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

def op_same_factors(opa: Operand, opb: Operand):
    return is_same_list(opa.factor_list, opb.factor_list)

G = Factor("G")
DELTA = Factor("Δ")
KAPPA = Factor("K")
GAMMA = Factor("γ")
NU = Factor("ν")
HBAR = Factor("h")

OP_A = Operand(1.0, [], [], [QOp.A])
OP_Ad = Operand(1.0, [], [], [QOp.Ad])
OP_N = OP_Ad * OP_A
OP_SP = Operand(1.0, [], [QOp.SP], [])
OP_SM = Operand(1.0, [], [QOp.SM], [])

OP_GSPA = Operand(1.0, [G], [QOp.SP], [QOp.A])
OP_GSMAd = Operand(1.0, [G], [QOp.SM], [QOp.Ad])

def hamiltonian_commutator(op: Operand):
    op_n_delta_i = OP_N.add_factor(DELTA).mul_c_factor(1.0j)
    return [op_n_delta_i * op, op * op_n_delta_i.mul_c_factor(-1.0), (OP_GSPA * op).mul_c_factor(1.0j), (OP_GSMAd * op).mul_c_factor(1.0j), (op * OP_GSPA).mul_c_factor(-1.0j), (op * OP_GSMAd).mul_c_factor(-1.0j)]
    
def lindblad_terms(lb_factor: Factor, lb_op: Operand, op: Operand):
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

def print_op_list(op_list):
    for i in range(len(op_list) - 1):
        print(op_list[i], " + ", end='')
    print(op_list[-1])

def remove_null(op_list):
    i = 0
    while i != len(op_list):
        if op_list[i].is_null():
            del op_list[i]
            continue
        i += 1

def order_cavity_qop(op: Operand):
    new_op_list = []
    for i in range(len(op.qop_cavity_list) - 1):
        if op.qop_cavity_list[i] == QOp.A and op.qop_cavity_list[i + 1] == QOp.Ad:
            op.qop_cavity_list[i] = QOp.Ad
            op.qop_cavity_list[i + 1] = QOp.A
            new_op = op.copy()
            del new_op.qop_cavity_list[i]
            del new_op.qop_cavity_list[i]
            new_op_list.append(new_op)
            new_op_list.extend(order_cavity_qop(new_op))
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

my_op = Operand(1, [], [], [QOp.Ad, QOp.A])
all_op = hamiltonian_commutator(Operand(1, [], [], [QOp.Ad, QOp.A]))
lb_a_terms = lindblad_terms(KAPPA, OP_A, my_op)
lb_sp_terms = lindblad_terms(GAMMA, OP_SP, my_op)
lb_sm_terms = lindblad_terms(NU, OP_SM, my_op)
all_op.extend(lb_a_terms)
all_op.extend(lb_sp_terms)
all_op.extend(lb_sm_terms)

print_op_list(all_op)

add_equivalent(all_op)

#print_op_list(all_op)

remove_null(all_op)

#print_op_list(all_op)

order_cavity_qop_list(all_op)

#print_op_list(all_op)

add_equivalent(all_op)

#print_op_list(all_op)

remove_null(all_op)

print_op_list(all_op)

test_list = [Operand(2, [], [QOp.SM, QOp.SZ, QOp.SP], []), Operand(-1, [], [QOp.SM, QOp.SP, QOp.SZ], []), Operand(-1, [], [QOp.SZ, QOp.SM, QOp.SP], [])]

print_op_list(test_list)

order_atomic_qop_list(test_list)

#print_op_list(test_list)

add_equivalent(test_list)

#print_op_list(test_list)

remove_null(test_list)

print_op_list(test_list)