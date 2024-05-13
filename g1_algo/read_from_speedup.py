import algo
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *

COUNT = 10
START_TIME = -0.0005150016 + (COUNT) * 0.00032
STOP_TIME = START_TIME + 2 * 0.00032
SAMPLE_COUNT = 200000
SAMPLE_SPACING = (STOP_TIME - START_TIME) / SAMPLE_COUNT
FREQ = 43
COHERENCE_TIME = 0
SIGNAL_COUNT = 50
DOWNSCALE_FACTOR = 200
TRUE_SAMPLE_COUNT = int(SAMPLE_COUNT / DOWNSCALE_FACTOR)


g1_norm = algo.import_g1_norm_darray("./g1_norm_speedup" + str(COUNT), TRUE_SAMPLE_COUNT)

t1_list = np.linspace(START_TIME, STOP_TIME, TRUE_SAMPLE_COUNT, endpoint=False)
t2_list = np.linspace(START_TIME, STOP_TIME, TRUE_SAMPLE_COUNT, endpoint=False)

fig0, ax0 = plt.subplots(1, 1)
fig0.set_size_inches(18.5, 10.5)
cs = ax0.contourf(t1_list * 1000, t2_list * 1000, g1_norm, cmap='inferno')
ax0.set_xlabel('time t1 (ms)')
ax0.set_ylabel('time t2 (ms)')
fig0.colorbar(cs)
title = "Start stop times = " + str(START_TIME) + " to " + str(STOP_TIME) + " with " + str(DOWNSCALE_FACTOR) + " - Sample count = " + str(SAMPLE_COUNT) + " - Frequency = " + str(round(FREQ, 2)) + " - Signal count = " + str(SIGNAL_COUNT) + " - BPF 5MHz"
ax0.set_title(title)
plt.savefig("./g1-2-analysis/" + title + ".png") ####### THINK TO ENABLE IT WHEN NECESSARY!!!
plt.show()

#from PIL import Image
#np_g1_norm = np.array(g1_norm) * 255
#np_g1_norm = np_g1_norm.astype(np.uint8)
#print(np_g1_norm)
#im = Image.fromarray(np_g1_norm)
#im.save("./img_" + title + ".jpg")
