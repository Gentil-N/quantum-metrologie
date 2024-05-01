import algo
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *
from scipy.fft import fft, ifft, fftfreq

START_SIGNAL = 1
STOP_SIGNAL = 10
PATH = "./../fastdiffamp/coh2/"
START_TIME = -8.000000E-04
STOP_TIME = 7.999968E-04
LOW_FREQ_BAND_PASS = 43E6#44.15E6
HIGH_FREQ_BAND_PASS = 45E6#44.25E6

#f_list = []
for i in range(START_SIGNAL, STOP_SIGNAL + 1):
    print("Processing ", i)
    print("Importing original signal...")
    f = algo.import_signal_csv_oscilloscope(PATH + "SC2CAL" + str(f'{i:02d}') + ".CSV")
    algo.plot_freq(f, START_TIME, STOP_TIME)
    exit()
    f = algo.bandpass_filter(f, START_TIME, STOP_TIME, LOW_FREQ_BAND_PASS, HIGH_FREQ_BAND_PASS)
    print("Creating hilbert pair...")
    g = algo.get_hilbert_pair(f, START_TIME, STOP_TIME)
    #f_list.append(f)
    print("Exporting f...")
    algo.export_signal("./f" + str(i), f)
    print("Exporting g...")
    algo.export_signal("./g" + str(i), g)

#for f in f_list[:-4]:
#    fig0, ax0 = plt.subplots(1, 1)
#    fig0.set_size_inches(18.5, 10.5)
#    ax0.plot(np.linspace(START_TIME, STOP_TIME, len(f)), f)
#plt.show()
