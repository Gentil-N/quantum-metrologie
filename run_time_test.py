import qutip as qp
import cavity as cty
from matplotlib import pyplot as plt
import numpy as np
import time

def run_test(cavity: cty.Cavity, atom_count, photon_capacity):
    cavity.init(1.0, 1.0, atom_count, photon_capacity, 10.0, 1.0, 4.0, True, True)
    #psi_init = qp.tensor(qp.spin_state(atom_count/2, atom_count/2), qp.fock(photon_capacity + 1, 0))
    psi_init = qp.tensor(qp.spin_state(atom_count/2, -atom_count/2), qp.coherent(photon_capacity + 1, 1))
    print("STARTING TEST: " + str(atom_count) + " atoms & " + str(photon_capacity) + " photon capacity")
    print("Running simulation... ", end='')
    start = time.clock_gettime(0)
    time_range = cavity.run_simulation(psi_init, 0, 5, 100)
    end = time.clock_gettime(0)
    print("Done in " + str(end - start) + "s")
    total_time = end - start
    print("Computing g1... ", end='')
    start = time.clock_gettime(0)
    g1 = cavity.compute_g1()
    end = time.clock_gettime(0)
    print("Done in " + str(end - start) + "s")
    total_time += end - start
    print("Computing g2... ", end='')
    start = time.clock_gettime(0)
    g2 = cavity.compute_g2()
    end = time.clock_gettime(0)
    print("Done in " + str(end - start) + "s")
    total_time += end - start
    print("Total time: " + str(total_time) + "s\n")

cavity = cty.Cavity(0.0, 0.0, 0, 0, 0.0, 0.0, 0.0, False, False)

for i in range(1, 10):
    for j in range(1, 4):
        run_test(cavity, i, j * 10)