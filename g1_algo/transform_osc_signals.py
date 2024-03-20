import algo
import numpy as np

START_SIGNAL = 1
STOP_SIGNAL = 10
PATH = "./../fastdiffamp/coh2/"
START_TIME = -8.000000E-04
STOP_TIME = 7.999968E-04

for i in range(START_SIGNAL, STOP_SIGNAL + 1):
    print("Processing ", i)
    print("Importing original signal...")
    f = algo.import_signal_csv_oscilloscope(PATH + "SC2CAL" + str(f'{i:02d}') + ".CSV")
    print("Creating hilbert pair...")
    g = algo.get_hilbert_pair(f, START_TIME, STOP_TIME)
    print("Exporting f...")
    algo.export_signal("./f" + str(i), f)
    print("Exporting g...")
    algo.export_signal("./g" + str(i), g)
