from matplotlib import pyplot as plt
import csv
import numpy as np

path_fst = "/home/nev/Documents/fastdiffamp/data-noise-set/first-try/"
path_sec = "/home/nev/Documents/fastdiffamp/data-noise-set/2nd-more-precision/"

def extract_data_from_csv_file(filename: str):
    file = open(filename)
    lines = file.readlines()
    if len(lines) < 2:
        raise Exception("Csv file error")
    elems = lines[1].split(",")
    start = float(elems[1])
    end = start + float(elems[2]) * (len(lines) - 3)
    y = []
    for i in range(2, len(lines)):
        y.append(float(lines[i].split(",")[0]))
    return [start, end, y]

def get_variance(data):
    n = len(data)
    mean = sum(data) / n
    variance = sum((x - mean) ** 2 for x in data) / (n - 1)
    return variance

data_list = []
for i in range(1, 11):
    data_list.append(extract_data_from_csv_file(path_sec + "NewFile" + str(i) + ".csv"))
time = np.linspace(data_list[0][0], data_list[0][1], len(data_list[0][2]))

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(18.5, 10.5)
var_list = []
for data in data_list:
    var_list.append(get_variance(data[2]))
    ax.plot(time, data[2])
ax.legend()
ax.set_xlabel('time')
ax.set_ylabel('correlations')

print(var_list)

plt.show()