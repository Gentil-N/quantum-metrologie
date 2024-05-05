import algo
import h5py

PATH_0 = "./../g1-1/"
PATH_1 = "./../g1-2/"
NAME = "G1Measurement_"
ADD_0 = ""
ADD_1 = "burst_"
ADD_2 = "oscillatory-burst_"
ADD_3 = "short-burst_"
ADD_1_0 = "bare_"
ADD_1_1 = "filt-demo_"

LOW_FREQ_BAND_PASS_0 = 1.75E6
HIGH_FREQ_BAND_PASS_0 = 2.05E6

LOW_FREQ_BAND_PASS_1 = 43.32E6
HIGH_FREQ_BAND_PASS_1 = 43.45E6

LOW_FREQ_BAND_PASS_2 = 41E6
HIGH_FREQ_BAND_PASS_2 = 46E6

def get_complete_path(root_path, name, add, num):
    return [root_path + name + add + str(num) + ".h5", root_path + name + add + str(num) + ".csv"]

def process(root_path, name, add, low_pass, high_pass):
    for i in range(50):
        print("Processing ", i)
        names = get_complete_path(root_path, name, add, i)
        fileh5 = h5py.File(names[0], "r")
        key_list = list(fileh5.keys())
        key_time = key_list[0]
        key_channels = key_list[1]
        time = fileh5[key_time][()].tolist()
        start_time = time[0]
        stop_time = time[-1]
        print(time[0] - time[100000])
        f = fileh5[key_channels][()][0].tolist()
        algo.plot_signal(f, start_time, stop_time)
        algo.plot_freq(f, start_time, stop_time)
        #exit()
        f = algo.bandpass_filter(f, start_time, stop_time, low_pass, high_pass)
        algo.plot_freq(f, start_time, stop_time)
        algo.plot_signal(f, start_time, stop_time)
        exit()
        print("Creating hilbert pair...")
        g = algo.get_hilbert_pair(f, start_time, stop_time)
        #f_list.append(f)
        print("Exporting f...")
        algo.export_signal("./temp_fg/f" + str(i), f)
        print("Exporting g...")
        algo.export_signal("./temp_fg/g" + str(i), g)
        print("done")

process(PATH_1, NAME, ADD_1_0, LOW_FREQ_BAND_PASS_2, HIGH_FREQ_BAND_PASS_2)
