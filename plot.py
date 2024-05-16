from matplotlib import pyplot as plt
import numpy as np
import math
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *

def plot_map_t1_greater_than_t_2():
    time = np.linspace(0, 1, 1000)
    plot_map = np.zeros((1000, 1000))
    for i in range(len(plot_map)):
        for j in range(len(plot_map[i])):
            if i >= j:
                plot_map[i][j] = 1.0

    fig0, ax0 = plt.subplots(1, 1)
    fig0.set_size_inches(18.5, 10.5)
    cs = ax0.contourf(time, time, plot_map, cmap='binary')
    ax0.set_xlabel('t1 [s]')
    ax0.set_ylabel('t2 [s]')
    fig0.colorbar(cs)
    plt.show()

def plot_notch_response():
    x = np.linspace(0, 100, 10000)
    y = 1/(np.sqrt(1 + (10 / (x - 1/x)**2)))
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(18.5, 10.5)
    ax1.plot(x, y)
    ax1.set_xlabel("ω [rad/s]")
    ax1.set_ylabel("Gain [no unit]")
    ax1.set_xscale('log')
    plt.show()

def plot_tank_response():
    x = np.linspace(0, 100, 10000)
    y = 1/(np.sqrt(1 + 10**2 * (1/x - x)**2))
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(18.5, 10.5)
    ax1.plot(x, y)
    ax1.set_xlabel("ω [rad/s]")
    ax1.set_ylabel("Gain [no unit]")
    ax1.set_xscale('log')
    plt.show()

def import_spec(filename):
    file = open(filename)
    lines = file.readlines()
    y = []
    x = []
    for line in lines:
        splitted = line.split(",")
        x.append(float(splitted[0]))
        y.append(float(splitted[1]))
    file.close()
    return [x, y]

def plot_specs():
    print("Importing..")
    odeep = import_spec("./fastdiffamp/original_deep")
    print("0")
    opick = import_spec("./fastdiffamp/original_pick")
    print("1")
    mdeep = import_spec("./fastdiffamp/modified_deep")
    print("2")
    mpick = import_spec("./fastdiffamp/modified_pick")
    print("done")
    fig1, ax1 = plt.subplots(1, 1)
    fig1.set_size_inches(18.5, 10.5)
    ax1.plot(np.array(odeep[0]) / 1000000, np.array(odeep[1]))
    ax1.set_xlabel("frequency [MHz]")
    ax1.set_ylabel("Power [dBm]")

    fig2, ax2 = plt.subplots(1, 1)
    fig2.set_size_inches(18.5, 10.5)
    ax2.plot(np.array(opick[0]) / 1000, opick[1])
    ax2.set_xlabel("frequency [kHz]")
    ax2.set_ylabel("Power [dBm]")

    fig3, ax3 = plt.subplots(1, 1)
    fig3.set_size_inches(18.5, 10.5)
    ax3.plot(np.array(mdeep[0]) / 1000000, mdeep[1])
    ax3.set_xlabel("frequency [MHz]")
    ax3.set_ylabel("Power [dBm]")

    fig4, ax4 = plt.subplots(1, 1)
    fig4.set_size_inches(18.5, 10.5)
    ax4.plot(np.array(mpick[0]) / 1000, mpick[1])
    ax4.set_xlabel("frequency [kHz]")
    ax4.set_ylabel("Power [dBm]")

    plt.show()

#plot_map_t1_greater_than_t_2()
#plot_notch_response()
#plot_tank_response()
plot_specs()
