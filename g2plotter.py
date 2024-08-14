import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import *
import math

def get_raw_array(filename):
    f = open(filename, 'r')
    content = f.readlines()
    list = []
    for x in content:
        row = x.split()
        list.append(int(row[0]))
    return np.array(list)

def get_bin_array(raw_array, time_range, bin_count, bin_size):
    array_bin = np.zeros(bin_count)
    for time in raw_array:
        array_bin[math.floor((time - raw_array[0]) / bin_size)] += 1
    time_list = np.linspace(0, time_range, bin_count)
    return [time_list, array_bin]

def get_times_array_bin(fl_list, bin_size):
    raw_array_list = [[], []]
    for filename in fl_list[0]:
        raw_array_list[0].append(get_raw_array(filename))
    for filename in fl_list[1]:
        raw_array_list[1].append(get_raw_array(filename))

    time_range_list = []
    for raw_array in raw_array_list[0]:
        time_range_list.append(raw_array[-1] - raw_array[0])
    for raw_array in raw_array_list[1]:
        time_range_list.append(raw_array[-1] - raw_array[0])
    time_range = max(time_range_list)
    bin_count = math.ceil(time_range / bin_size)

    print("time range:", time_range, "/ bin count:", bin_count)

    times_list = [[], []]
    array_bin_list = [[], []]
    for raw_array in raw_array_list[0]:
        res = get_bin_array(raw_array, time_range, bin_count, bin_size)
        times_list[0].append(res[0])
        array_bin_list[0].append(res[1])
    for raw_array in raw_array_list[1]:
        res = get_bin_array(raw_array, time_range, bin_count, bin_size)
        times_list[1].append(res[0])
        array_bin_list[1].append(res[1])
    return [times_list, array_bin_list]

def plot_times_array_bin(times_array_bin_list):
    times_list = times_array_bin_list[0]
    array_bin_list = times_array_bin_list[1]
    for i in range(len(times_list[0])):
        #print(times_list[i][int(len(times_list) / 2)])
        plt.plot(times_list[0][i] / 1E9, array_bin_list[0][i])
        plt.plot(times_list[1][i] / 1E9, array_bin_list[1][i])
        plt.ylabel("photon count")
        plt.xlabel("time (ms)")
        plt.title("Photon count for each data array (2024/07/03 - 21h19)")
    plt.show()

def get_g2(times_array_bin_list):
    times_list = times_array_bin_list[0]
    array_bin_list = times_array_bin_list[1]
    N = len(times_list[0])
    # Assumption: all "time_list" are equivalents
    half_denominator_first = []
    half_denominator_second = []
    for i in range(len(times_list[0][0])):
        temp = 0.0
        for array_bin in array_bin_list[0]:
            temp += array_bin[i]
        half_denominator_first.append(temp)
    for i in range(len(times_list[1][0])):
        temp = 0.0
        for array_bin in array_bin_list[1]:
            temp += array_bin[i]
        half_denominator_second.append(temp)
    map_t1 = []
    for i in range(len(times_list[0][0])):
        print(int(times_list[0][0][i] / times_list[0][0][-1] * 100), "%")
        map_t2 = []
        for j in range(len(times_list[0][0])):
            temp = 0.0
            for k in range(len(array_bin_list[0])):
                temp += array_bin_list[0][k][i] * array_bin_list[1][k][j]
            map_t2.append(N * temp / (half_denominator_first[i] * half_denominator_second[j]))
        map_t1.append(map_t2)
    return map_t1


def get_stationary_g2(g2, bound):
    N = len(g2)
    g2_stationary = [0.0] * int(bound * 2 + 1)
    for i in range(N):
        for j in range(N):
            pos = abs(i - j)
            if (((i < 2*bound) and ((j < 2*bound - i) or (j > 2*bound + i))) or ((i >= 2*bound) and ((j < i - 2*bound) or (j > 6*bound - i)))):
                continue
            g2_stationary[pos] += g2[i][j] / int(bound * 4 + 1)
    g2_stationary[0] *= 2
    return g2_stationary

def get_stationary_g2_lub(g2, l_b, u_b):
    N = len(g2)
    u_b_flip = (u_b - N / 2) * -1.0 + N / 2
    u_b_is_min = (u_b_flip == min(u_b_flip, l_b))
    len_stat = 2 * u_b_flip if u_b_is_min else 2 * l_b
    g2_stationary = [0.0] * int(len_stat + 1)
    print(u_b, l_b)
    for i in range(N):
        for j in range(N):
            pos = abs(i - j)
            if (j > (2 * u_b - i)) or (j < (2 * l_b - i)):
                continue
            if u_b_is_min:
                if (j < (i - 2 * u_b_flip)) or (j > (i + 2 * u_b_flip)):
                    continue
            else:
                if (j < (i - 2 * l_b)) or (j > (i + 2 * l_b)):
                    continue
            g2_stationary[pos] += g2[i][j] / int((u_b - l_b) * 2 + 1)
    g2_stationary[0] *= 2
    return g2_stationary


### SCRIPT

fl_19h14 = [['./results/comp/data_19h14_ch2_trig2',],
            ['./results/comp/data_19h14_ch1_trig2']]

fl_19h17 = [['./results/comp/data_19h17_ch1_trig0',
            './results/comp/data_19h17_ch1_trig1',
            './results/comp/data_19h17_ch1_trig2',
            './results/comp/data_19h17_ch1_trig3'],
            ['./results/comp/data_19h17_ch2_trig0',
            './results/comp/data_19h17_ch2_trig1',
            './results/comp/data_19h17_ch2_trig2',
            './results/comp/data_19h17_ch2_trig3']]

fl_15h23 = [['./results/comp/data_15h23_ch1_trig0',
            './results/comp/data_15h23_ch1_trig1',
            './results/comp/data_15h23_ch1_trig2',
            './results/comp/data_15h23_ch1_trig3',],
            ['./results/comp/data_15h23_ch2_trig0',
            './results/comp/data_15h23_ch2_trig1',
            './results/comp/data_15h23_ch2_trig2',
            './results/comp/data_15h23_ch2_trig3']]

fl_21h19 = [['./results/20240703/data_21h19_ch1_trig0',
            './results/20240703/data_21h19_ch1_trig2',
            './results/20240703/data_21h19_ch1_trig4',
            './results/20240703/data_21h19_ch1_trig10',
            './results/20240703/data_21h19_ch1_trig11',
            './results/20240703/data_21h19_ch1_trig12',
            './results/20240703/data_21h19_ch1_trig16',
            './results/20240703/data_21h19_ch1_trig18',],
            ['./results/20240703/data_21h19_ch2_trig0',
            './results/20240703/data_21h19_ch2_trig2',
            './results/20240703/data_21h19_ch2_trig4',
            './results/20240703/data_21h19_ch2_trig10',
            './results/20240703/data_21h19_ch2_trig11',
            './results/20240703/data_21h19_ch2_trig12',
            './results/20240703/data_21h19_ch2_trig16',
            './results/20240703/data_21h19_ch2_trig18']]

fl_19h20 = [['./results/20240702/data_19h20_ch1_trig3',
            './results/20240702/data_19h20_ch1_trig8',],
            ['./results/20240702/data_19h20_ch2_trig3',
            './results/20240702/data_19h20_ch2_trig8']]

fl_19h40 = [['./results/20240702/data_19h40_ch1_trig3',
            './results/20240702/data_19h40_ch1_trig8',],
            ['./results/20240702/data_19h40_ch2_trig3',
            './results/20240702/data_19h40_ch2_trig8']]

fl_19h42 = [['./results/20240702/data_19h42_ch1_trig2',
            './results/20240702/data_19h42_ch1_trig5',
            './results/20240702/data_19h42_ch1_trig10',
            './results/20240702/data_19h42_ch1_trig13',
            './results/20240702/data_19h42_ch1_trig15',
            './results/20240702/data_19h42_ch1_trig18',],
            ['./results/20240702/data_19h42_ch2_trig2',
            './results/20240702/data_19h42_ch2_trig5',
            './results/20240702/data_19h42_ch2_trig10',
            './results/20240702/data_19h42_ch2_trig13',
            './results/20240702/data_19h42_ch2_trig15',
            './results/20240702/data_19h42_ch2_trig18',]]

fl_19h48 = [['./results/20240702/data_19h48_ch1_trig2',
            './results/20240702/data_19h48_ch1_trig9',
            './results/20240702/data_19h48_ch1_trig16',],
            ['./results/20240702/data_19h48_ch2_trig2',
            './results/20240702/data_19h48_ch2_trig9',
            './results/20240702/data_19h48_ch1_trig16',]]

fl_19h50 = [['./results/20240702/data_19h50_ch1_trig1',
            './results/20240702/data_19h50_ch1_trig5',],
            ['./results/20240702/data_19h50_ch1_trig1',
            './results/20240702/data_19h50_ch1_trig5']]

BIN_SIZE = 50_000_000 # in picoseconds!
res = get_times_array_bin(fl_21h19, BIN_SIZE)
plot_times_array_bin(res)

g2 = get_g2(res)

fig0, ax0 = plt.subplots(1, 1)
fig0.set_size_inches(18.5, 10.5)
cs = ax0.contourf(res[0][0][0] / 1E9, res[0][0][0] / 1E9, g2, cmap='inferno')
ax0.set_ylabel('t2 (ms)')
ax0.set_xlabel('t1 (ms)')
ax0.set_title("|g2(t1,t2)| (2024/07/03 - 21h19)")
fig0.colorbar(cs)
plt.show()


#bound = int(len(g2) / 4)
#g2_stationary = get_stationary_g2(g2, bound)
range_time = (res[0][0][0][-1] - res[0][0][0][0]) / 1E9 # ms
l_time = 1.25 # ms
u_time = 1.45 # ms
l_b = int(l_time * len(g2) / range_time)
u_b = int(u_time * len(g2) / range_time)

g2_stationary = get_stationary_g2_lub(g2, l_b, u_b)
#time_stationary = np.linspace(0, res[0][0][0][bound], len(g2_stationary))
time_stationary = np.linspace(0, len(g2_stationary) * range_time / len(g2), len(g2_stationary))
plt.plot(time_stationary, g2_stationary)
plt.ylim(bottom=0)
plt.xlabel("tau (ms)")
plt.ylabel("|g2(tau)|")
plt.title("Average of |g2(t1,t2)| from t1 = 1.25 ms to t1 = 1.45 ms")
plt.show()