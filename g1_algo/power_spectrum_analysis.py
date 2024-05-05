import algo
import h5py
import numpy as np
from matplotlib import pyplot as plt


power_spectrum_list = []
freqs = []
CUT_START = 400000
CUT_STOP = 1300000

LOW_FREQ_BAND_PASS = 41E6
HIGH_FREQ_BAND_PASS = 46E6
sample_delay = 0
start_time = 0
stop_time = 0

for i in range(50):
    print("Processing ", i)
    fileh5 = h5py.File("./../g1-2/G1Measurement_bare_" + str(i) + ".h5", "r")
    key_list = list(fileh5.keys())
    key_time = key_list[0]
    key_channels = key_list[1]
    time = fileh5[key_time][()].tolist()
    start_time = time[0]
    stop_time = time[-1]
    f = fileh5[key_channels][()][0]
    f = algo.bandpass_filter(f, start_time, stop_time, LOW_FREQ_BAND_PASS, HIGH_FREQ_BAND_PASS)
    #algo.plot_spectrogram(algo.get_spectrogram(f, start_time, stop_time, 100, 4.325E7, 4.35E7), start_time, stop_time, str(i))
    cut_f = algo.cut_signal(f, CUT_START, CUT_STOP)
    algo.plot_power_spectrum(f, time, 0.0009, 0.0025, 4.33E7, 4.345E7, str(i))
    #algo.plot_freq(f, start_time, stop_time)
    #sample_delay = (stop_time - start_time) / len(time)
    #algo.plot_freq(cut_f, start_time + CUT_START * sample_delay, start_time + CUT_STOP * sample_delay)
    #exit()
    sample_delay = (stop_time - start_time) / len(time)
    freq_data = algo.get_fft(cut_f, start_time + CUT_START * sample_delay, start_time + CUT_STOP * sample_delay)
    power_spectrum_list.append(np.abs(freq_data[1])**2)
    freqs = freq_data[0]
    print("done")

average_power_spectrum = power_spectrum_list[0]
for i in range(1, 50):
    average_power_spectrum += power_spectrum_list[i] / 50


fig1, ax1 = plt.subplots(1, 1)
fig1.set_size_inches(18.5, 10.5)
ax1.plot(freqs[:int(len(freqs)/2)], average_power_spectrum[:int(len(freqs)/2)]) #plot half + title with sample/time cut
ax1.set_title("Average power spectrum over 50 signals - 5MHz band pass filter - cut " + str(CUT_START) + " (" + str(start_time + CUT_START * sample_delay) + "):" + str(CUT_STOP) + " (" + str(start_time + CUT_STOP * sample_delay) + ")")
ax1.legend()
plt.show()
