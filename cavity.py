import qutip as qp
from matplotlib import pyplot as plt
import numpy as np
import math

HBAR = 1.0

class Cavity:

    def __init__(self, cavity_mode_frequency, atomic_transition_frequency, atom_count, photon_capacity, coupling, kappa = 0.0, repump = 0.0, enable_decay = False, enable_repump = False) -> None:
        self.init(cavity_mode_frequency, atomic_transition_frequency, atom_count, photon_capacity, coupling, kappa, repump, enable_decay, enable_repump)

    def init(self, cavity_mode_frequency, atomic_transition_frequency, atom_count, photon_capacity, coupling, kappa = 0.0, repump = 0.0, enable_decay = False, enable_repump = False):

        self.cf = cavity_mode_frequency
        self.af = atomic_transition_frequency
        self.a_count = atom_count
        self.spin_num = self.a_count / 2.0
        self.spin_dim = self.a_count + 1 # 2 * N / 2 + 1
        self.photon_capacity = photon_capacity
        self.fock_dim = self.photon_capacity + 1
        self.g = coupling
        self.kappa = kappa
        self.repump = repump
        self.enable_decay = enable_decay
        self.enable_repump = enable_repump

        self.op_jz = qp.spin_Jz(self.spin_num)
        self.op_jp = qp.spin_Jp(self.spin_num)
        self.op_jm = qp.spin_Jm(self.spin_num)

        self.op_ad = qp.create(self.fock_dim)
        self.op_a = qp.destroy(self.fock_dim)
        self.op_n = self.op_ad * self.op_a

        #self.op_hamiltonian_base = HBAR * self.mf * self.__tens_from_fock(self.op_n) + HBAR * self.af * self.__tens_from_spin(self.op_jz)
        #self.op_hamiltonian_intercation = HBAR * self.g / np.sqrt(self.a_count) * (qp.tensor(self.op_jp, self.op_a) + qp.tensor(self.op_jm, self.op_ad))
        #self.op_hamiltonian = self.op_hamiltonian_base + self.op_hamiltonian_intercation

        self.op_hamiltonian_base = HBAR * (self.cf - self.af) * self.__tens_from_fock(self.op_n)
        self.op_hamiltonian_intercation = HBAR * self.g / np.sqrt(self.a_count) * (qp.tensor(self.op_jp, self.op_a) + qp.tensor(self.op_jm, self.op_ad))
        self.op_hamiltonian = self.op_hamiltonian_base + self.op_hamiltonian_intercation

        self.op_collapsing = []
        if self.enable_decay:
            self.op_collapsing.append(self.kappa * self.__tens_from_fock(self.op_a))
        if self.enable_repump:
            self.op_collapsing.append(self.repump * self.__tens_from_spin(self.op_jp))

        self.result = []
        self.time_range = []

    def run_simulation(self, psi_init, time_start, time_stop, time_num):
        if psi_init.type == 'ket':
            psi_init = psi_init * psi_init.dag()
        elif psi_init.type == 'bra':
            psi_init = psi_init.dag() * psi_init
        self.time_range = np.linspace(time_start, time_stop, time_num)
        result = qp.mesolve(self.op_hamiltonian, psi_init, self.time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar())
        self.result = result.states
        return self.time_range

    def get_all_probabilities(self):
        probabilities = []
        m = -self.spin_num
        while m != self.spin_num + 1:
            for n in range(0, self.fock_dim, 1):
                ket = qp.tensor(qp.spin_state(self.spin_num, m), qp.fock(self.fock_dim, n))
                prob = [qp.expect(ket * ket.dag(), self.result), "j=" + str(self.spin_num) + ", m=" + str(m) + ", n=" + str(n)]
                probabilities.append(prob)
            m += 1
        return probabilities
    
    def compute_g2(self):
        g2 = []
        i = 0
        for t in self.time_range:
            op_u = self.__get_u_op_from_hamiltonian(t)
            op_ad_evol_heis = op_u.dag() * self.__tens_from_fock(self.op_ad) * op_u
            op_a_evol_heis = op_u.dag() * self.__tens_from_fock(self.op_a) * op_u
            norm_factor = (self.result[i] * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)).tr()**2
            if norm_factor == 0.0:
                g2.append(0.0)
                i += 1
                continue
            g2.append((self.result[i] * self.__tens_from_fock(self.op_ad) * op_ad_evol_heis * op_a_evol_heis * self.__tens_from_fock(self.op_a)).tr() / norm_factor)
            i += 1
        return g2
    
    def compute_g1(self):
        g1 = []
        i = 0
        for t in self.time_range:
            op_u = self.__get_u_op_from_hamiltonian(t)
            op_ad_evol_heis = op_u.dag() * self.__tens_from_fock(self.op_ad) * op_u
            op_a_evol_heis = op_u.dag() * self.__tens_from_fock(self.op_a) * op_u
            norm_factor = np.sqrt((self.result[i] * op_ad_evol_heis * op_a_evol_heis).tr() * (self.result[i] * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)).tr())
            if norm_factor == 0.0:
                g1.append(0.0)
                i += 1
                continue
            g1.append((self.result[i] * op_ad_evol_heis * self.__tens_from_fock(self.op_a)).tr() / norm_factor)
            i += 1
        return g1

    def get_state(self):
        return self.result

    def __get_u_op_from_hamiltonian(self, time):
        return (-1j * self.op_hamiltonian * time / HBAR).expm()

    def __tens_from_fock(self, op):
        return qp.tensor(qp.qeye(self.spin_dim), op)
    
    def __tens_from_spin(self, op):
        return qp.tensor(op, qp.qeye(self.fock_dim))
    
