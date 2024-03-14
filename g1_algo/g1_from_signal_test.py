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
            g1_submap.append(g1_t1_t2)
        print("t1=", t1, " done!")
        g1_map.append(g1_submap)
    return g1_map


def compute_single_numerator_by_index(f, g, t, tau):
    return np.sqrt((f[t]*f[t+tau] + g[t]*g[t+tau])**2 + (g[t]*f[t+tau] - f[t]*g[t+tau])**2)

def compute_single_denominator_by_index(f, g, t):
    return f[t]**2 + g[t]**2

def compute_stationary_g1_norm(f, g, max_tau_sample_count, sample_count, sample_spacing):
    g1 = []
    x = np.linspace(0, sample_count - max_tau_sample_count, sample_count - max_tau_sample_count, endpoint=False)
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = compute_single_denominator_by_index(f, g, x[i])
    denominator = integrate.simps(y, x)
    #print("denominator : ", denominator * sample_spacing, " ", denominator)
    for tau in reversed(range(max_tau_sample_count)):
        #stop_sample = sample_count - tau
        x = np.linspace(0, sample_count - max_tau_sample_count, sample_count - max_tau_sample_count, endpoint=False)
        y = np.zeros(len(x))
        for i in range(len(x)):
            y[i] = compute_single_numerator_by_index(f, g, int(x[i]), int(tau))
        numerator = integrate.simps(y, x)
        #print("numerator : ", numerator * sample_spacing, " ", numerator)
        g1.append(numerator / denominator)
        print("tau=", tau, " done!", g1[-1])
    return g1


### GENERATING CHAOTIC LIGHT
#start_time = 0.0
#stop_time = 50.0
#sample_count = 200
#sample_spacing = (stop_time - start_time) / sample_count
#freq = 2.0 * np.pi
#coherence_time = 2.0 # factor multiplied to the random function
#amount_of_signals = 10
#phase_list = [+np.pi, -np.pi]
#f_list = []
#average_coherence = 0
#num_gen = 0
#for i in range(amount_of_signals):
#    current_signal = []
#    current_sample = 0
#    stop = False
#    while current_sample != sample_count:
#        shift_sample = int(random.random() * coherence_time * 2 / sample_spacing)
#        average_coherence += shift_sample * sample_spacing
#        num_gen += 1
#        phase = 2 * np.pi * random.random() - np.pi
#        if shift_sample + current_sample > sample_count:
#            shift_sample = sample_count - current_sample
#        #print("\n", current_sample, " ", shift_sample)
#        current_time = current_sample * sample_spacing + start_time
#        x = np.linspace(current_time, current_time + shift_sample * sample_spacing, shift_sample, endpoint=stop)
#        #print(x)
#        y = np.cos(freq * x + freq * current_time + phase)
#        current_sample += shift_sample
#        current_signal.extend(y)
#    f_list.append(current_signal)
#### END OF GENERATION
#
#print(average_coherence/num_gen)
#for f in f_list[:3]:
#    x = np.linspace(start_time, stop_time, sample_count, endpoint=False)
#    fig1, ax1 = plt.subplots(1, 1)
#    fig1.set_size_inches(18.5, 10.5)
#    ax1.plot(x, f)
#    ax1.grid()
##exit()
#
#
##assuming that all f signals have equals sample_count, start and end time
#sample_count = len(f_list[0]) # also enumerates the possibilities for t1 and t2
#g_list = create_hilbert_pair_for_list(f_list, start_time, stop_time)
#g_map = compute_g1_norm(f_list, g_list, sample_count)
#
#t1_list = np.linspace(start_time, stop_time, sample_count, endpoint=False)
#t2_list = np.linspace(start_time, stop_time, sample_count, endpoint=False)
#
#g_map_np_array = np.array(g_map)
#print(g_map_np_array.max())
#print(g_map_np_array.min())
#
#fig1, ax1 = plt.subplots(1, 1)
#fig1.set_size_inches(18.5, 10.5)
#ax1.plot(t2_list, g_map[0])
#ax1.set_xlabel('time t2 for t1 = 0')
#ax1.legend()
#ax1.grid()
#
#fig2, ax2 = plt.subplots(1, 1)
#fig2.set_size_inches(18.5, 10.5)
#ax2.plot(t2_list, g_map[50])
#ax2.set_xlabel('time t2 for t1 = 500')
#ax2.legend()
#ax2.grid()
#
#fig3, ax3 = plt.subplots(1, 1)
#fig3.set_size_inches(18.5, 10.5)
#ax3.plot(t2_list, g_map[73])
#ax3.set_xlabel('time t2 for t1 = 730')
#ax3.legend()
#ax3.grid()
#
#fig0, ax0 = plt.subplots(1, 1)
#fig0.set_size_inches(18.5, 10.5)
#cs = ax0.contourf(t1_list, t2_list, g_map, cmap='inferno')
#ax0.set_xlabel('time t1')
#ax0.set_ylabel('time t2')
#fig0.colorbar(cs)
##fig0.colorbar(mappable=ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap='inferno'), ax=ax0)
#
#plt.show()


def get_exp_ditributed_rand(mean):
    return -1.0 * mean * np.log(random.random())

def get_unif_distributed_rand(mean):
    return random.random() * mean * 2

### GENERATING SINGLE CHAOTIC LIGHT
start_time = 0.0
stop_time = 4000.0
sample_count = 20000
sample_spacing = (stop_time - start_time) / sample_count
freq = 2.0 * np.pi
coherence_time = 4.0 # factor multiplied to the random function
amount_of_signals = 1
f_list = []
x = np.linspace(start_time, stop_time, sample_count, endpoint=False)
#f_list = [np.cos(freq*x)]
#z_signal = np.exp(1j * freq * x)
mean_check = 0.0
num_loop = 0
for i in range(amount_of_signals):
    current_signal = []
    current_sample = 0
    stop = False
    while current_sample != sample_count:
        mean_check += get_exp_ditributed_rand(coherence_time)
        num_loop += 1
        shift_sample = int(get_exp_ditributed_rand(coherence_time) / sample_spacing) # or get_unif_distributed_rand(coherence_time)
        phase = 2 * np.pi * random.random() - np.pi
        if shift_sample + current_sample > sample_count:
            shift_sample = sample_count - current_sample
        #print("\n", current_sample, " ", shift_sample)
        current_time = current_sample * sample_spacing + start_time
        x = np.linspace(current_time, current_time + shift_sample * sample_spacing, shift_sample, endpoint=stop)
        #print(x)
        #y = np.cos(freq * x + freq * current_time + phase)
        y = np.cos(freq * x + freq * current_time + phase) + 1j * np.sin(freq * x + freq * current_time + phase)
        current_sample += shift_sample
        current_signal.extend(y)
    f_list.append(current_signal)
g = get_hilbert_pair(np.real(f_list[0]), start_time, stop_time)
print("mean check = ", mean_check / num_loop)
### END OF GENERATION

def test_compute_g1(re_signal, im_signal, max_tau_sample_count, sample_count):
    g1 = []
    x = np.linspace(0, sample_count - max_tau_sample_count, sample_count - max_tau_sample_count, endpoint=False)
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = re_signal[int(x[i])]**2 + im_signal[int(x[i])]**2
    denominator = integrate.simps(y, x)
    for tau in range(max_tau_sample_count):
        #stop_sample = sample_count - tau
        yr = np.zeros(len(x))
        yi = np.zeros(len(x))
        #y = np.zeros(len(x))
        for i in range(len(x)):
            #y[i] = np.sqrt((re_signal[int(x[i])]*re_signal[int(x[i]+tau)] + im_signal[int(x[i])]*im_signal[int(x[i]+tau)])**2 + (im_signal[int(x[i])]*re_signal[int(x[i]+tau)] - re_signal[int(x[i])]*im_signal[int(x[i]+tau)])**2)
            #cx_num = (re_signal[int(x[i])] - 1j * im_signal[int(x[i])]) * (re_signal[int(x[i] + tau)] + 1j * im_signal[int(x[i] + tau)])
            #yr[i] = np.real(cx_num)
            #yi[i] = np.imag(cx_num)
            yr[i] = re_signal[int(x[i])]*re_signal[int(x[i]+tau)] + im_signal[int(x[i])]*im_signal[int(x[i]+tau)]
            yi[i] = im_signal[int(x[i])]*re_signal[int(x[i]+tau)] - re_signal[int(x[i])]*im_signal[int(x[i]+tau)]
        #numerator = integrate.simps(y, x)
        numerator = np.sqrt(integrate.simps(yr, x)**2 + integrate.simps(yi, x)**2)
        #print("numerator : ", numerator)
        g1.append(numerator / denominator)
        print("tau=", tau, " done!", g1[-1])
    return g1

x = np.linspace(start_time, stop_time, sample_count, endpoint=False)
fig1, ax1 = plt.subplots(1, 1)
fig1.set_size_inches(18.5, 10.5)
#ax1.plot(x, f_list[0], "+")
#ax1.plot(x, f_list[0])
ax1.plot(x, np.imag(f_list[0]))
ax1.plot(x, g)
ax1.grid()
plt.show()
#exit()

def jorg_stationary_test(f_list, start_time, stop_time):
    sample_count = len(f_list[0])
    freq_list = fftfreq(sample_count, (stop_time - start_time) / sample_count)
    average_squared_spectrum = np.zeros(len(freq_list))
    signal_count = len(f_list)
    for i in range(signal_count):
        average_squared_spectrum += np.abs(fft(f_list[i])) / signal_count
    return np.abs(ifft(np.sqrt(average_squared_spectrum)))

#assuming that all f signals have equals sample_count, start and end time
#sample_count = len(f_list[0]) # also enumerates the possibilities for t1 and t2
#print(sample_count)
jorg_test = jorg_stationary_test(f_list, start_time, stop_time)
max_tau_sample_count = int(sample_count / 10)
print("max_tau_sample_count", max_tau_sample_count)
g_single = test_compute_g1(np.real(f_list[0]), g, max_tau_sample_count, sample_count)
#g_list = create_hilbert_pair_for_list(f_list, start_time, stop_time)
#g_single = compute_stationary_g1_norm(f_list[0], g_list[0], max_tau_sample_count, sample_count, sample_spacing)


x = np.linspace(0, max_tau_sample_count, max_tau_sample_count) * sample_spacing
theo_lorentzian = np.exp(-np.abs(x)/coherence_time)
fig1, ax1 = plt.subplots(1, 1)
fig1.set_size_inches(18.5, 10.5)
ax1.plot(x, g_single, label="simulation")
ax1.plot(x, theo_lorentzian, label="theory")
#ax1.plot(x, jorg_test[:max_tau_sample_count], label="Jörg's stationary test")
ax1.grid()
ax1.legend()
ax1.set_title("Start stop times = " + str(start_time) + " to " + str(stop_time) + " - Sample count = " + str(sample_count) + " - Coh time = " + str(coherence_time) + " - Frequency = " + str(freq))
plt.savefig("./g1 from signal plots/Start stop times = " + str(start_time) + " to " + str(stop_time) + " - Sample count = " + str(sample_count) + " - Coh time = " + str(coherence_time) + " - Frequency = " + str(freq) +".png") 
plt.show()
