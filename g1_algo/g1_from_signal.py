from matplotlib import pyplot as plt
from scipy.fft import fft, ifft, fftfreq
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *
import random
from scipy import integrate


def get_hilbert_pair(signal, start_time, stop_time):
    sample_count = len(signal)
    freq_list = fftfreq(sample_count, (stop_time - start_time) / sample_count)
    return np.imag(ifft(np.sign(freq_list) * fft(signal)))

def create_hilbert_pair_for_list(f_list, start_time, stop_time):
    g_list = []
    for f in f_list:
        g_list.append(get_hilbert_pair(f, start_time, stop_time))
    return g_list

def compute_single_numerator_part_one(f, g, t1, t2):
    return f[t1] * f[t2] + g[t1] * g[t2]

def compute_single_numerator_part_two(f, g, t1, t2):
    return g[t1] * f[t2] - f[t1] * g[t2]

def compute_single_denominator_single_part(f, g, t):
    return f[t]**2 + g[t]**2

def compute_g1_norm(f_list, g_list, sample_count):
    len_f_list = len(f_list)
    assert(len_f_list == len(g_list))
    g1_map = []
    for t1 in range(sample_count):
        g1_submap = []
        for t2 in range(sample_count):
            numerator_part_one = 0.0
            numerator_part_two = 0.0
            denominator_t1 = 0.0
            denominator_t2 = 0.0
            for i in range(len_f_list):
                current_f = f_list[i]
                current_g = g_list[i]
                numerator_part_one += compute_single_numerator_part_one(current_f, current_g, t1, t2)
                numerator_part_two += compute_single_numerator_part_two(current_f, current_g, t1, t2)
                denominator_t1 += compute_single_denominator_single_part(current_f, current_g, t1)
                denominator_t2 += compute_single_denominator_single_part(current_f, current_g, t2)
            numerator_part_one /= len_f_list
            numerator_part_two /= len_f_list
            denominator_t1 /= len_f_list
            denominator_t2 /= len_f_list
            denominator = np.sqrt(denominator_t1 * denominator_t2)
            g1_t1_t2 = np.sqrt(numerator_part_one**2 + numerator_part_two**2) / denominator
            if 0 > g1_t1_t2 :
                exit()
            g1_submap.append(g1_t1_t2)
        print("t1=", t1, " done!")
        g1_map.append(g1_submap)
    return g1_map


### GENERATING CHAOTIC LIGHT
start_time = 0.0
stop_time = 200.0
sample_count = 2000
sample_spacing = (stop_time - start_time) / sample_count
freq = 2.0 * np.pi
coherence_time = 10.0 # factor multiplied to the random function
amount_of_signals = 20
f_list = []
x = np.linspace(start_time, stop_time, sample_count, endpoint=False)
#f_list = [np.cos(freq*x)]
#z_signal = np.exp(1j * freq * x)
for i in range(amount_of_signals):
    current_signal = []
    current_sample = 0
    stop = False
    while current_sample != sample_count:
        shift_sample = int(random.random() * coherence_time * 2 / sample_spacing)
        phase = 2 * np.pi * random.random() - np.pi
        if shift_sample + current_sample > sample_count:
            shift_sample = sample_count - current_sample
        #print("\n", current_sample, " ", shift_sample)
        current_time = current_sample * sample_spacing + start_time
        x = np.linspace(current_time, current_time + shift_sample * sample_spacing, shift_sample, endpoint=stop)
        #print(x)
        #y = np.cos(freq * x + freq * current_time + phase)
        y = np.cos(freq * x + freq * current_time + phase)
        current_sample += shift_sample
        current_signal.extend(y)
    f_list.append(current_signal)
g_list = create_hilbert_pair_for_list(f_list, start_time, stop_time)
### END OF GENERATION

for f in f_list[:3]:
    x = np.linspace(start_time, stop_time, sample_count, endpoint=False)
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(18.5, 10.5)
    ax1.plot(x, f)
    ax1.grid()
plt.show()

g_list = create_hilbert_pair_for_list(f_list, start_time, stop_time)
g_map = compute_g1_norm(f_list, g_list, sample_count)

t1_list = np.linspace(start_time, stop_time, sample_count, endpoint=False)
t2_list = np.linspace(start_time, stop_time, sample_count, endpoint=False)

fig0, ax0 = plt.subplots(1, 1)
fig0.set_size_inches(18.5, 10.5)
cs = ax0.contourf(t1_list, t2_list, g_map, cmap='inferno')
ax0.set_xlabel('time t1')
ax0.set_ylabel('time t2')
fig0.colorbar(cs)
#fig0.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='inferno'), ax=ax0)

plt.show()