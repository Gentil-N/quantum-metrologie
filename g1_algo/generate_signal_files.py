import algo
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *

#START_TIME = 0.0
#STOP_TIME = 100.0
#SAMPLE_COUNT = 10000
#SAMPLE_SPACING = (STOP_TIME - START_TIME) / SAMPLE_COUNT
#FREQ = 2.0 * np.pi
#COHERENCE_TIME = 2.0
#SIGNAL_COUNT = 40

#f_list = algo.generate_chaotic_light_list(START_TIME, STOP_TIME, SAMPLE_COUNT, SAMPLE_SPACING, COHERENCE_TIME, FREQ, algo.get_exp_ditributed_rand, SIGNAL_COUNT)
#g_list = algo.create_hilbert_pair_for_list(f_list, START_TIME, STOP_TIME)

START_TIME = 0.0
STOP_TIME = 3.2E-05
SAMPLE_COUNT = 10000
DELTA_OMEGA = 43E6 # 43 MHz
DELTA_PHI = 0.0
ABS_ALPHA = 1.0
ABS_BETA = 1.0
SIGNAL_COUNT = 10

f_list = algo.generate_coherent_light_heterodyne(START_TIME, STOP_TIME, SAMPLE_COUNT, DELTA_OMEGA, DELTA_PHI, ABS_ALPHA, ABS_BETA, SIGNAL_COUNT)
g_list = algo.create_hilbert_pair_for_list(f_list, START_TIME, STOP_TIME)

fig0, ax0 = plt.subplots(1, 1)
fig0.set_size_inches(18.5, 10.5)
ax0.plot(np.linspace(START_TIME, STOP_TIME, SAMPLE_COUNT), f_list[0])
ax0.plot(np.linspace(START_TIME, STOP_TIME, SAMPLE_COUNT), f_list[1])
plt.show()

print("Generating f list")
for i in range(len(f_list)):
    algo.export_signal("./f" + str(i + 1), f_list[i])
print("Generating g list")
for i in range(len(g_list)):
    algo.export_signal("./g" + str(i + 1), g_list[i])
