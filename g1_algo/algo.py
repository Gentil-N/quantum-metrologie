from scipy.fft import fft, ifft, fftfreq
import random
import numpy as np
from scipy import integrate

def import_signal(filename: str):
    file = open(filename)
    lines = file.readlines()
    y = []
    for line in lines:
        y.append(float(line))
    file.close()
    return y

def import_g1_norm_darray(filename: str, row_size):
    print("Start importing g1_norm_array...")
    file = open(filename)
    y = []
    count = 0
    total_processed = 0
    temp = []
    while True:
        line = file.readline()
        if not line:
            break
        temp.append(float(line))
        count += 1
        if count == row_size:
            y.append(temp)
            temp = []
            count = 0
            print("Imported: ", total_processed)
            total_processed += 1
    file.close()
    print("Import done")
    return y

def export_signal(filename: str, signal):
    file = open(filename, "w")
    for sample in signal:
        file.write(str(sample) + "\n")
    file.close()

def get_hilbert_pair(signal, start_time, stop_time):
    sample_count = len(signal)
    freq_list = fftfreq(sample_count, (stop_time - start_time) / sample_count)
    return np.imag(ifft(np.sign(freq_list) * fft(signal)))

def create_hilbert_pair_for_list(f_list, start_time, stop_time):
    g_list = []
    for f in f_list:
        g_list.append(get_hilbert_pair(f, start_time, stop_time))
    return g_list

def _compute_single_numerator_part_one(f, g, t1, t2):
    return f[t1] * f[t2] + g[t1] * g[t2]

def _compute_single_numerator_part_two(f, g, t1, t2):
    return g[t1] * f[t2] - f[t1] * g[t2]

def compute_g1_norm(f_list, g_list):
    sample_count = len(f_list[0])
    len_f_list = len(f_list)
    assert(len_f_list == len(g_list))
    denominator_list = []
    for t in range(sample_count):
        denominator = 0.0
        for i in range(len_f_list):
            denominator += (f_list[i][t]**2 + g_list[i][t]**2)
        denominator_list.append(denominator / len_f_list)
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
                numerator_part_one += _compute_single_numerator_part_one(current_f, current_g, t1, t2)
                numerator_part_two += _compute_single_numerator_part_two(current_f, current_g, t1, t2)
            numerator_part_one /= len_f_list
            numerator_part_two /= len_f_list
            denominator = np.sqrt(denominator_list[t1] * denominator_list[t2])
            g1_t1_t2 = np.sqrt(numerator_part_one**2 + numerator_part_two**2) / denominator
            g1_submap.append(g1_t1_t2)
        print("t1 = ", t1, " over ", sample_count-1 ," done!")
        g1_map.append(g1_submap)
    return g1_map

def compute_stationary_g1_norm(re_signal, im_signal, max_tau_sample_count, sample_count):
    g1 = []
    x = np.linspace(0, sample_count - max_tau_sample_count, sample_count - max_tau_sample_count, endpoint=False)
    y = np.zeros(len(x))
    for i in range(len(x)):
        y[i] = re_signal[int(x[i])]**2 + im_signal[int(x[i])]**2
    denominator = integrate.simps(y, x)
    for tau in range(max_tau_sample_count):
        yr = np.zeros(len(x))
        yi = np.zeros(len(x))
        for i in range(len(x)):
            yr[i] = re_signal[int(x[i])]*re_signal[int(x[i]+tau)] + im_signal[int(x[i])]*im_signal[int(x[i]+tau)]
            yi[i] = im_signal[int(x[i])]*re_signal[int(x[i]+tau)] - re_signal[int(x[i])]*im_signal[int(x[i]+tau)]
        numerator = np.sqrt(integrate.simps(yr, x)**2 + integrate.simps(yi, x)**2)
        g1.append(numerator / denominator)
        print("tau = ", tau, " over ", max_tau_sample_count, " done!")
    return g1

def get_exp_ditributed_rand(mean):
    return -1.0 * mean * np.log(random.random())

def get_unif_distributed_rand(mean):
    return random.random() * mean * 2

def generate_chaotic_light_list(start_time, stop_time, sample_count, sample_spacing, coherence_time, freq, distribution_func, signal_count):
    f_list = []
    x = np.linspace(start_time, stop_time, sample_count, endpoint=False)
    mean_check = 0.0
    num_loop = 0
    for i in range(signal_count):
        current_signal = []
        current_sample = 0
        while current_sample != sample_count:
            rand_val = distribution_func(coherence_time)
            shift_sample = int(rand_val / sample_spacing)
            mean_check += rand_val
            num_loop += 1
            phase = 2 * np.pi * random.random() - np.pi
            if shift_sample + current_sample > sample_count:
                shift_sample = sample_count - current_sample
            current_time = current_sample * sample_spacing + start_time
            x = np.linspace(current_time, current_time + shift_sample * sample_spacing, shift_sample, endpoint=False)
            y = np.cos(freq * x + freq * current_time + phase) # + 1j * np.sin(freq * x + freq * current_time + phase)
            current_sample += shift_sample
            current_signal.extend(y)
        f_list.append(current_signal)
    print("mean check = ", mean_check / num_loop)
    return f_list

def jorg_stationary_test(f_list, start_time, stop_time):
    sample_count = len(f_list[0])
    freq_list = fftfreq(sample_count, (stop_time - start_time) / sample_count)
    average_squared_spectrum = np.zeros(len(freq_list))
    signal_count = len(f_list)
    for i in range(signal_count):
        average_squared_spectrum += (np.abs(fft(f_list[i])) / signal_count)**2
    return np.abs(ifft(np.sqrt(average_squared_spectrum)))
