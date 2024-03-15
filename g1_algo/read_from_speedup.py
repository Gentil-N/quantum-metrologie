import algo
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *

START_TIME = 0.0
STOP_TIME = 100.0
SAMPLE_COUNT = 1000
SAMPLE_SPACING = (STOP_TIME - START_TIME) / SAMPLE_COUNT
FREQ = 2.0 * np.pi
COHERENCE_TIME = 2.0
SIGNAL_COUNT = 40


g1_norm = algo.import_g1_norm_darray("./g1_norm_speedup", SAMPLE_COUNT)

t1_list = np.linspace(START_TIME, STOP_TIME, SAMPLE_COUNT, endpoint=False)
t2_list = np.linspace(START_TIME, STOP_TIME, SAMPLE_COUNT, endpoint=False)

fig0, ax0 = plt.subplots(1, 1)
fig0.set_size_inches(18.5, 10.5)
cs = ax0.contourf(t1_list, t2_list, g1_norm, cmap='inferno')
ax0.set_xlabel('time t1')
ax0.set_ylabel('time t2')
fig0.colorbar(cs)
title = "[SPEED_UP] Start stop times = " + str(START_TIME) + " to " + str(STOP_TIME) + " - Sample count = " + str(SAMPLE_COUNT) + " - Coh time = " + str(COHERENCE_TIME) + " - Frequency = " + str(FREQ) + " - Signal count = " + str(SIGNAL_COUNT)
ax0.set_title(title)
plt.savefig("./" + title + ".png")
plt.show()
