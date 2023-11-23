import qutip as qp
from matplotlib import pyplot as plt
import numpy as np
import math

HBAR = 1.0



class Cavity:

    def __init__(self, cavity_mode_frequency, atomic_transition_frequency, atom_count, photon_capacity, coupling, kappa = 0.0, repump = 0.0, free_space_rate = 0.0, enable_decay = False, enable_repump = False, enable_free_space_decay = False) -> None:
        self.init(cavity_mode_frequency, atomic_transition_frequency, atom_count, photon_capacity, coupling, kappa, repump, free_space_rate, enable_decay, enable_repump, enable_free_space_decay)

    def init(self, cavity_mode_frequency, atomic_transition_frequency, atom_count, photon_capacity, coupling, kappa = 0.0, repump = 0.0, free_space_rate = 0.0, enable_decay = False, enable_repump = False, enable_free_space_decay = False):

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
        self.free_space_rate = free_space_rate
        self.enable_decay = enable_decay
        self.enable_repump = enable_repump
        self.enable_free_space_decay = enable_free_space_decay

        self.op_jz = qp.spin_Jz(self.spin_num)
        self.op_jp = qp.spin_Jp(self.spin_num)
        self.op_jm = qp.spin_Jm(self.spin_num)

        self.op_ad = qp.create(self.fock_dim)
        self.op_a = qp.destroy(self.fock_dim)
        self.op_n = self.op_ad * self.op_a

        # "Standard" hamiltonian
        #self.op_hamiltonian_base = HBAR * self.cf * self.__tens_from_fock(self.op_n) + HBAR * self.af * self.__tens_from_spin(self.op_jz)
        #self.op_hamiltonian_intercation = HBAR * self.g / np.sqrt(self.a_count) * (qp.tensor(self.op_jp, self.op_a) + qp.tensor(self.op_jm, self.op_ad)) # sqrt(N) !!! Ok for dicke model. But must be verified before use!
        #self.op_hamiltonian = self.op_hamiltonian_base + self.op_hamiltonian_intercation

        # RWA
        self.op_hamiltonian_base = HBAR * (self.cf - self.af) * self.__tens_from_fock(self.op_n)
        self.op_hamiltonian_intercation = HBAR * self.g / np.sqrt(self.a_count) * (qp.tensor(self.op_jp, self.op_a) + qp.tensor(self.op_jm, self.op_ad)) # sqrt(N) !!! Ok for dicke model. But must be verified before use!
        self.op_hamiltonian = self.op_hamiltonian_base + self.op_hamiltonian_intercation

        self.op_collapsing = []
        if self.enable_decay:
            self.op_collapsing.append(np.sqrt(self.kappa) * self.__tens_from_fock(self.op_a)) 
        if self.enable_repump:
            self.op_collapsing.append(np.sqrt(self.repump) * self.__tens_from_spin(self.op_jp))
        if self.enable_free_space_decay:
            self.op_collapsing.append(np.sqrt(self.free_space_rate) * self.__tens_from_spin(self.op_jm))

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

    def compute_steady_state(self):
        return qp.steadystate(self.op_hamiltonian, self.op_collapsing)

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

    def compute_expectation_values(self, intial_index_result):
        return qp.mesolve(self.op_hamiltonian, self.result[intial_index_result], self.time_range, self.op_collapsing, [self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a), self.__tens_from_spin(self.op_jz), self.__tens_from_spin(self.op_jp)], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar()).expect
    
    def compute_g2(self):
        g2 = []
        i = 0
        op_ad_evol_heis_list = qp.mesolve(-self.op_hamiltonian, self.__tens_from_fock(self.op_ad), self.time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar()).states
        op_a_evol_heis_list = qp.mesolve(-self.op_hamiltonian, self.__tens_from_fock(self.op_a), self.time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000),progress_bar=qp.ui.progressbar.TextProgressBar()).states
        #op_n_evol_heis_list = qp.mesolve(-self.op_hamiltonian, self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a), self.time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar()).states
        n = qp.mesolve(self.op_hamiltonian, self.result[0], self.time_range, self.op_collapsing, [self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar()).expect[0]

        propagator_list = qp.propagator(self.op_hamiltonian, self.time_range, self.op_collapsing)

        #Open Quantum System
        op_n_evol_list = []
        for prop in propagator_list:
            op_n_evol_list.append(qp.vector_to_operator(prop * qp.operator_to_vector(self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a))))
        op_ad_evol_list = []
        for prop in propagator_list:
            op_ad_evol_list.append(qp.vector_to_operator(prop * qp.operator_to_vector(self.__tens_from_fock(self.op_ad))))
        op_a_evol_list = []
        for prop in propagator_list:
            op_a_evol_list.append(qp.vector_to_operator(prop * qp.operator_to_vector(self.__tens_from_fock(self.op_a))))
        
        #Closed Quantum System
        #op_n_evol_list = []
        #for prop in propagator_list:
        #    op_n_evol_list.append(prop.dag() * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a) * prop)
        #op_ad_evol_list = []
        #for prop in propagator_list:
        #    op_ad_evol_list.append(prop.dag() * self.__tens_from_fock(self.op_ad) * prop)
        #op_a_evol_list = []
        #for prop in propagator_list:
        #    op_a_evol_list.append(prop.dag() * self.__tens_from_fock(self.op_a)  * prop)
        
        #propagator_list = qp.propagator(self.op_hamiltonian, self.time_range, self.op_collapsing)
        #for j in range(len(n)):
        #    #if not self.__op_is_approx_null(qp.vector_to_operator(propagator_list[j] * qp.operator_to_vector(self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a))) - #op_ad_evol_heis_list[j] * op_a_evol_heis_list[j], 0.01):
        #    #    print(j)
        #    if not self.__op_is_approx_null(qp.vector_to_operator(propagator_list[j] * qp.operator_to_vector(self.__tens_from_fock(self.op_ad))) * qp.vector_to_operator(propagator_list[j] * #qp.operator_to_vector(self.__tens_from_fock(self.op_a))) - op_n_evol_heis_list[j], 0.01):
        #        print(j)
        #exit()
        for t in self.time_range:
            #op_u = self.__get_u_op_from_hamiltonian(t)
            #op_ad_evol_heis = op_u.dag() * self.__tens_from_fock(self.op_ad) * op_u
            #op_a_evol_heis = op_u.dag() * self.__tens_from_fock(self.op_a) * op_u
            #op_ad_evol_heis = op_ad_evol_heis_list[i]
            #op_a_evol_heis = op_a_evol_heis_list[i]
            result = self.result[0]
            op_n_evol = op_n_evol_list[i]
            op_ad_evol = op_ad_evol_heis_list[i]
            op_a_evol = op_a_evol_heis_list[i]
            #norm_factor = (self.result[i] * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)).tr()**2
            #norm_factor = (result * op_ad_evol_heis * op_a_evol_heis).tr() * (result * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)).tr()
            norm_factor = 1.0#n[0] * n[i]
            if norm_factor == 0.0:
                g2.append(0.0)
                i += 1
                continue
            g2.append((result * self.__tens_from_fock(self.op_ad) * op_ad_evol * op_a_evol * self.__tens_from_fock(self.op_a)).tr() / norm_factor)
            i += 1
        return g2

    def compute_g1(self):
        g1 = []
        i = 0
        op_ad_evol_heis_list = qp.mesolve(-self.op_hamiltonian, self.__tens_from_fock(self.op_ad), self.time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar()).states
        #op_a_evol_heis_list = qp.mesolve(-self.op_hamiltonian, self.__tens_from_fock(self.op_a), self.time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000),progress_bar=qp.ui.progressbar.TextProgressBar()).states
        n = qp.mesolve(self.op_hamiltonian, self.result[0], self.time_range, self.op_collapsing, [self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar()).expect[0]

        propagator_list = qp.propagator(self.op_hamiltonian, self.time_range, self.op_collapsing)

        #Open Quantum System
        op_ad_evol_list = []
        for prop in propagator_list:
            op_ad_evol_list.append(qp.vector_to_operator(prop * qp.operator_to_vector(self.__tens_from_fock(self.op_ad))))
        op_a_evol_list = []
        for prop in propagator_list:
            op_a_evol_list.append(qp.vector_to_operator(prop * qp.operator_to_vector(self.__tens_from_fock(self.op_a))))

        #Closed Quantum System
        #op_ad_evol_list = []
        #for prop in propagator_list:
        #    op_ad_evol_list.append(prop.dag() * self.__tens_from_fock(self.op_ad) * prop)
        #op_a_evol_list = []
        #for prop in propagator_list:
        #    op_a_evol_list.append(prop.dag() * self.__tens_from_fock(self.op_a)  * prop)
        
        for t in self.time_range:
            op_u = self.__get_u_op_from_hamiltonian(t)
            op_ad_evol = op_u.dag() * self.__tens_from_fock(self.op_ad) * op_u
            op_a_evol = op_u.dag() * self.__tens_from_fock(self.op_a) * op_u
            #op_ad_evol = op_ad_evol_list[i]
            #op_a_evol = op_a_evol_list[i]
            result = self.result[0]
            norm_factor = np.sqrt((result * op_ad_evol * op_a_evol).tr() * (result * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)).tr())
            #norm_factor = np.sqrt((self.result[0] * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)).tr() * (self.result[i] * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)).tr())
            #norm_factor = np.sqrt(n[0] * n[i])
            #print(n[i], (result * op_ad_evol * op_a_evol).tr())
            if norm_factor == 0.0:
                #print(op_ad_evol_heis * op_a_evol_heis)
                #print(n[i], " ", (self.result[0] * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)).tr())
                g1.append(0.0)
                i += 1
                continue
            g1.append((result * op_ad_evol * self.__tens_from_fock(self.op_a)).tr() / 1.0)
            i += 1
        return g1

    def compute_g1_qutip(self, psi_init, time_range):
        g1 = qp.coherence_function_g1(self.op_hamiltonian, psi_init, time_range, self.op_collapsing, self.__tens_from_fock(self.op_a))[0]
        return g1

    def compute_g2_qutip(self, psi_init, time_range):
        #g2 = qp.correlation_3op_1t(self.op_hamiltonian, psi_init, self.time_range, self.op_collapsing, self.__tens_from_fock(self.op_ad), self.__tens_from_fock(self.op_ad) * self.#__tens_from_fock(self.op_a), self.__tens_from_fock(self.op_a))
        #for i in range(len(self.time_range)):
        #    norm_factor = (self.result[i] * self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)).tr()**2
        #    g2[i] /= norm_factor
        g2 = qp.coherence_function_g2(self.op_hamiltonian, psi_init, time_range, self.op_collapsing, self.__tens_from_fock(self.op_a))[0]
        return g2

    def compute_g2_array_qutip(self, psi_init, time_range, taulist):
        print("g2 array qutip...", end="")
        g2_array = qp.correlation_3op_2t(self.op_hamiltonian, psi_init, time_range, taulist, self.op_collapsing, self.__tens_from_fock(self.op_ad), self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a), self.__tens_from_fock(self.op_a))
        print("Done")
        return g2_array

    def compute_g2_array(self, psi_init, time_range, taulist):
        psi_init = self.__convert_psi_init(psi_init)
        g2_array = []
        i = 0
        for t in time_range:
            print("computing at time t = ", t)
            n = qp.mesolve(self.op_hamiltonian, self.result[i], t + taulist, self.op_collapsing, [self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)], options=qp.solver.Options(nsteps=1000)).expect[0]
            DpA = qp.mesolve(self.op_hamiltonian, self.__tens_from_fock(self.op_a) * self.result[i] * self.__tens_from_fock(self.op_ad), t + taulist, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000))
            BCDpA = self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a) * DpA.states
            g2 = []
            for BCDpA_unique in BCDpA:
                g2.append(BCDpA_unique.tr())
            g2 /= n[0] * n
            g2_array.append(np.real(g2))
            i += 1
        return g2_array

    def compute_g1_perso(self, psi_init, time_range):
        psi_init = self.__convert_psi_init(psi_init)
        #result = qp.mesolve(self.op_hamiltonian, psi_init, time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar())
        Bp = qp.mesolve(self.op_hamiltonian, self.__tens_from_fock(self.op_a) * psi_init, time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar())
        n = qp.mesolve(self.op_hamiltonian, psi_init, time_range, self.op_collapsing, [self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar()).expect[0]
        ABp = self.__tens_from_fock(self.op_ad) * Bp.states
        g1 = []
        for ABp_unique in ABp:
            g1.append(ABp_unique.tr())
        #g1 /= np.sqrt(n[0] * n)
        return g1

    def compute_g2_perso(self, psi_init, time_range):
        psi_init = self.__convert_psi_init(psi_init)
        #result = qp.mesolve(self.op_hamiltonian, psi_init, time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar())
        DpA = qp.mesolve(self.op_hamiltonian, self.__tens_from_fock(self.op_a) * psi_init * self.__tens_from_fock(self.op_ad), time_range, self.op_collapsing, [], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar())
        #print(DpA.states)
        n = qp.mesolve(self.op_hamiltonian, psi_init, time_range, self.op_collapsing, [self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a)], options=qp.solver.Options(nsteps=1000), progress_bar=qp.ui.progressbar.TextProgressBar()).expect[0]
        BCDpA = self.__tens_from_fock(self.op_ad) * self.__tens_from_fock(self.op_a) * DpA.states
        g2 = []
        for BCDpA_unique in BCDpA:
            g2.append(BCDpA_unique.tr())
        #g2 /= n[0] * n
        return g2

    def get_state(self):
        return self.result

    def __get_u_op_from_hamiltonian(self, time):
        return (-1j * self.op_hamiltonian * time / HBAR).expm()

    def __tens_from_fock(self, op):
        return qp.tensor(qp.qeye(self.spin_dim), op)
    
    def __tens_from_spin(self, op):
        return qp.tensor(op, qp.qeye(self.fock_dim))
    
    def __convert_psi_init(self, psi_init):
        if psi_init.type == 'ket':
            return psi_init * psi_init.dag()
        elif psi_init.type == 'bra':
            return psi_init.dag() * psi_init
        return psi_init

    def __op_is_approx_null(self, op, precision):
        for row in op.full():
            for num in row:
                if num <= -precision or num >= precision:
                    return False
        return True