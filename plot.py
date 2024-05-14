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

#plot_map_t1_greater_than_t_2()
plot_notch_response()
