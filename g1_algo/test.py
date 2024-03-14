import algo
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *

START_TIME = 0.0
STOP_TIME = 4000.0
SAMPLE_COUNT = 40_000
SAMPLE_SPACING = (STOP_TIME - START_TIME) / SAMPLE_COUNT
MAX_TAU_SAMPLE_COUNT = int(SAMPLE_COUNT / 40)
FREQ = 2.0 * np.pi
COHERENCE_TIME = 10.0
SIGNAL_COUNT = 1

f_list = algo.generate_chaotic_light_list(START_TIME, STOP_TIME, SAMPLE_COUNT, SAMPLE_SPACING, COHERENCE_TIME, FREQ, algo.get_exp_ditributed_rand, SIGNAL_COUNT)
g_list = algo.create_hilbert_pair_for_list(f_list, START_TIME, STOP_TIME)

jorg_test = algo.jorg_stationary_test(f_list, START_TIME, STOP_TIME)
MAX_TAU_SAMPLE_COUNT = int(SAMPLE_COUNT / 10)
print("MAX_TAU_SAMPLE_COUNT ", MAX_TAU_SAMPLE_COUNT)
g_single = algo.compute_stationary_g1_norm(f_list[0], g_list[0], MAX_TAU_SAMPLE_COUNT, SAMPLE_COUNT)


x = np.linspace(0, MAX_TAU_SAMPLE_COUNT, MAX_TAU_SAMPLE_COUNT) * SAMPLE_SPACING
theo_lorentzian = np.exp(-np.abs(x)/COHERENCE_TIME)
fig1, ax1 = plt.subplots(1, 1)
fig1.set_size_inches(18.5, 10.5)
ax1.plot(x, g_single, label="simulation")
ax1.plot(x, theo_lorentzian, label="theory")
ax1.legend()
ax1.grid()
title = "Start stop times = " + str(START_TIME) + " to " + str(STOP_TIME) + " - Sample count = " + str(SAMPLE_COUNT) + " - Coh time = " + str(COHERENCE_TIME) + " - Frequency = " + str(FREQ)
ax1.set_title(title)
plt.savefig("./g1_algo/" + title + ".png")
plt.show()
