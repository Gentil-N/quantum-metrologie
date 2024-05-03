from scipy.fft import fft, ifft, fftfreq
import random
import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *

def import_signal(filename: str):
    file = open(filename)
    lines = file.readlines()
    y = []
    for line in lines:
        y.append(float(line))
    file.close()
    return y

def import_signal_csv_oscilloscope(filename: str):
    file = open(filename)
    lines = file.readlines()
    y = []
    for line in lines[1:]:
        y.append(float(line.split(",")[1]))
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

def bandpass_filter(signal, start_time, stop_time, low_freq_limit, high_freq_limit):
    sample_count = len(signal)
    freq_list = fftfreq(sample_count, (stop_time - start_time) / sample_count)
    fsignal = fft(signal)
    indices = np.where(((freq_list <= low_freq_limit) & (freq_list >= -low_freq_limit)) | (freq_list >= high_freq_limit) | (freq_list <= -high_freq_limit))
    fsignal[indices] = 0
    return np.real(ifft(fsignal))

def get_fft(signal, start_time, stop_time):
    sample_count = len(signal)
    return [fftfreq(sample_count, (stop_time - start_time) / sample_count), fft(signal)]

def cut_signal(signal, sample_start, sample_stop):
    return signal[sample_start:sample_stop].copy()

def get_spectrogram(signal, start_time, stop_time, amount, low_freq_limit, high_freq_limit):
    assert(len(signal) % amount == 0)
    sample_count_per_spec = int(len(signal) / amount)
    freqs = fftfreq(sample_count_per_spec, (stop_time - start_time) / amount / sample_count_per_spec)
    indices = np.where(((freqs >= low_freq_limit) & (freqs <= high_freq_limit)))
    freqs = freqs[indices]
    spec_list = []
    for i in range(amount):
        spec_list.append(np.abs(fft(signal[int(i*sample_count_per_spec):int((i+1)*sample_count_per_spec)]))[indices]**2)
        #fig1, ax1 = plt.subplots(1, 1)
        #fig1.set_size_inches(18.5, 10.5)
        #x1.plot(freqs, spec_list[-1], label="freq")
        #ax1.set_yscale('log')
        #ax1.legend()
        #plt.show()
    return [freqs, spec_list]

def plot_spectrogram(spectrogram, start_time, stop_time, name):
    fig0, ax0 = plt.subplots(1, 1)
    fig0.set_size_inches(18.5, 10.5)
    time_list = np.linspace(start_time, stop_time, len(spectrogram[1]))
    print(len(time_list), " ", len(spectrogram[0]), " ", len(spectrogram[1]), " ", len(spectrogram[1][0]))
    cs = ax0.contourf(spectrogram[0] / 1E6, time_list * 1000, spectrogram[1], cmap='inferno')
    ax0.set_xlabel('frequencies (MHZ)')
    ax0.set_ylabel('time (ms)')
    fig0.colorbar(cs)
    plt.savefig("./g1-2-analysis/spectrograms/spec-"+name+".png")
    plt.show()

def plot_power_spectrum(signal, time, cut_start, cut_stop, freq_cut_start, freq_cut_stop, name):
    time_indices = np.where((np.array(time) >= cut_start) & (np.array(time) <= cut_stop))
    time = np.array(time)[time_indices]
    signal = signal[time_indices]
    print(len(time), " ", len(signal))
    freqs = fftfreq(len(time), (cut_stop - cut_start) / len(time))
    fft_signal = fft(signal)
    freq_indices = np.where((freqs > freq_cut_start) & (freqs < freq_cut_stop))
    freqs = freqs[freq_indices]
    fft_signal = fft_signal[freq_indices]
    power_spectrum = np.abs(fft_signal)**2
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(18.5, 10.5)
    ax1.plot(freqs / 1E6, power_spectrum, label="freq")
    ax1.set_xlabel("frequency (MHz)")
    ax1.legend()
    plt.savefig("./g1-2-analysis/power-spectrums-continous/ps-cut(" + str(cut_start) + "-" + str(cut_stop) + ")-" + name + ".png")
    plt.show()


def plot_signal(signal, start_time, stop_time):
    sample_count = len(signal)
    time = np.linspace(start_time, stop_time, sample_count)
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(18.5, 10.5)
    ax1.plot(time, signal, label="freq")
    ax1.legend()
    plt.show()

def plot_freq(signal, start_time, stop_time):
    sample_count = len(signal)
    freq_list = fftfreq(sample_count, (stop_time - start_time) / sample_count)
    fsignal = fft(signal)
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(18.5, 10.5)
    index = int(len(freq_list)/2)
    ax1.plot(freq_list[:index], np.real(fsignal[:index]), label="freq")
    ax1.legend()
    plt.show()

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
    #x = np.linspace(start_time, stop_time, sample_count, endpoint=False)
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

def generate_coherent_light_heterodyne(start_time, stop_time, sample_count, delta_omega, delta_phi, abs_alpha, abs_beta, signal_count):
    f_list = []
    x = np.linspace(start_time, stop_time, sample_count, endpoint=False)
    for i in range(signal_count):
        current_delta_phi = delta_phi
        if signal_count > 1:
            current_delta_phi = random.random() * 2.0 * np.pi
        y = 2.0 * np.real(abs_alpha * abs_beta * np.exp(1j*(delta_omega * x + current_delta_phi)))
        f_list.append(y)
    return f_list

def jorg_stationary_test(f_list, start_time, stop_time):
    sample_count = len(f_list[0])
    freq_list = fftfreq(sample_count, (stop_time - start_time) / sample_count)
    average_squared_spectrum = np.zeros(len(freq_list))
    signal_count = len(f_list)
    for i in range(signal_count):
        average_squared_spectrum += (np.abs(fft(f_list[i])) / signal_count)**2
    return np.abs(ifft(np.sqrt(average_squared_spectrum)))
