import algo
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *
import math

OFFSET_TIME = 0.00064
SAMPLE_COUNT = 200000
DOWNSCALE_FACTOR = 200
TRUE_SAMPLE_COUNT = int(SAMPLE_COUNT / DOWNSCALE_FACTOR) #nb of row/col pixels in the image

def compute_stationary_g1_chunk(g1_norm):
    N = len(g1_norm)
    assert(N == len(g1_norm[0]))
    g1_stationary_chunk = [0.0] * int(N/2 + 1) # +1 to be sure to include every pixel
    for i in range(N):
        for j in range(N):
            pos = abs(i - j)
            if (i >= N/2 and (j < (i - N/2) or j > (-i + (3*N/2)))) or (i < N/2 and (j < (N/2 - i) or j > (N/2 + i))):
                continue
            g1_stationary_chunk[pos] += g1_norm[i][j] / int(N + 1)
    g1_stationary_chunk[0] *= 2 #because by abs(), index 0 is run only once
    return g1_stationary_chunk


index_start = 2
index_count = 8
g1_norm = algo.import_g1_norm_darray("./g1_norm_speedup" + str(index_start), TRUE_SAMPLE_COUNT)
g1_stationary = np.array(compute_stationary_g1_chunk(g1_norm)) / (index_count + 1)
for i in range(index_count):
    print("Processing:", i)
    g1_norm = algo.import_g1_norm_darray("./g1_norm_speedup" + str(index_start + i + 1), TRUE_SAMPLE_COUNT)
    g1_stationary += np.array(compute_stationary_g1_chunk(g1_norm)) / (index_count + 1)


#g1_norm = algo.import_g1_norm_darray("./g1_norm_speedup" + str(INDEX), TRUE_SAMPLE_COUNT)
#fig1, ax1 = plt.subplots(1, 1)
#fig1.set_size_inches(18.5, 10.5)
#ax1.plot(np.linspace(0, 1, len(g1_norm[0])), g1_norm[0])
#ax1.legend()
#plt.show()
#g1_stationary = compute_stationary_g1_chunk(g1_norm)

#tau_max = len(g1_stationary) * HALF_OFFSET_TIME * 2 / SAMPLE_COUNT
print(len(g1_stationary))
tau_list = np.linspace(0, 0.00016, len(g1_stationary))

fig1, ax1 = plt.subplots(1, 1)
fig1.set_size_inches(18.5, 10.5)
ax1.plot(tau_list * 1000, g1_stationary**2)
ax1.set_xlabel("τ (ms)")
ax1.set_ylabel("|g¹(τ)|²")
ax1.legend()
plt.show()
