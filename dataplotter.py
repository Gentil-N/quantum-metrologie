from matplotlib import pyplot as plt
import csv
import numpy as np

path_fst = "/home/nev/Documents/fastdiffamp/data-noise-set/first-try/"
path_sec = "/home/nev/Documents/fastdiffamp/data-noise-set/2nd-more-precision/"
path_thd = "/home/nev/Documents/fastdiffamp/data-noise-set/3th - 100khz low pow filter/"
path_fourth = "/home/nev/Documents/fastdiffamp/data-noise-set/4th - 3 files/"

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
for i in range(1, 28):
    data_list.append(extract_data_from_csv_file(path_fourth + "NewFile" + str(i) + ".csv"))
time = np.linspace(data_list[0][0], data_list[0][1], len(data_list[0][2]))

#fig, ax = plt.subplots(1, 1)
#fig.set_size_inches(18.5, 10.5)
#var_list = []
#for data in data_list:
#    var_list.append(get_variance(data[2]))
#    ax.plot(time, data[2])
#ax.legend()
#ax.set_xlabel('time')

var_list = []
for data in data_list:
    var_list.append(get_variance(data[2]))

power_data = [5.80, 5.16, 4.40, 3.85, 3.11, 2.19, 1.40, 0.89, 0.42]
pw_dt_cpy = power_data.copy()
varlist_fst = []
varlist_sec = []
varlist_thd = []
for i in range(0, len(var_list), 3):
    print(i + 2)
    varlist_fst.append(var_list[i])
    varlist_sec.append(var_list[i + 1])
    varlist_thd.append(var_list[i + 2])

pw_dt_cpy.extend(pw_dt_cpy)
pw_dt_cpy.extend(power_data.copy())
print(pw_dt_cpy)
var_tot = varlist_fst.copy()
var_tot.extend(varlist_sec)
var_tot.extend(varlist_thd)
z_fst = np.polyfit(pw_dt_cpy, var_tot, 1)
#z_sec = np.polyfit(power_data, varlist_sec, 1)
#z_thd = np.polyfit(power_data, varlist_thd, 1)

print(z_fst)

x = np.linspace(power_data[0], power_data[-1], 100)
y_fst = z_fst[1] + x * z_fst[0]
#y_sec = z_sec[1] + x * z_sec[0]
#y_thd = z_thd[1] + x * z_thd[0]

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(18.5, 10.5)
ax.plot(power_data, varlist_fst, ".", color="green")
ax.plot(x, y_fst, "--", color="green")
ax.plot(power_data, varlist_sec, ".", color="blue")
#ax.plot(x, y_sec, "--", color="blue")
ax.plot(power_data, varlist_thd, ".", color="red")
#ax.plot(x, y_thd, "--", color="red")
ax.legend()
ax.set_xlabel('time')

plt.show()