from matplotlib import pyplot as plt
import csv
import numpy as np

path_fst = "/home/nev/Documents/fastdiffamp/data-noise-set/first-try/"
path_sec = "/home/nev/Documents/fastdiffamp/data-noise-set/2nd-more-precision/"
path_thd = "/home/nev/Documents/fastdiffamp/data-noise-set/3th - 100khz low pow filter/"
path_fourth = "/home/nev/Documents/fastdiffamp/data-noise-set/4th - 3 files/"
path_fifth = "/home/nev/Documents/fastdiffamp/data-noise-set/5th/"
path_sixth = "/home/nev/Documents/fastdiffamp/data-noise-set/6th/"
path_seventh = "/home/nev/Documents/fastdiffamp/data-noise-set/7th/"
path_eighth = "/home/nev/Documents/fastdiffamp/data-noise-set/8th/"

path_weak_balanced = "/home/nev/Documents/fastdiffamp/data-noise-set/weak-balanced/"
path_strong_balanced = "/home/nev/Documents/fastdiffamp/data-noise-set/strong-balanced/"
path_strong_unb = "/home/nev/Documents/fastdiffamp/data-noise-set/strong-unb/"
path_str_wk_balanced = "/home/nev/Documents/fastdiffamp/data-noise-set/str-wk-balanced/"
path_str_wk_unb = "/home/nev/Documents/fastdiffamp/data-noise-set/str-wk-unb/"

path_unb_left = "/home/nev/Documents/fastdiffamp/data-noise-set/unb-left-9to925/"
path_unb_right = "/home/nev/Documents/fastdiffamp/data-noise-set/unb-right-9to915/"
path_bal = "/home/nev/Documents/fastdiffamp/data-noise-set/bal-9to919/"

path_bal_demod = "/home/nev/Documents/fastdiffamp/data-noise-set/bal-demod/"
path_unb_left_demod = "/home/nev/Documents/fastdiffamp/data-noise-set/unb-left-demod/"
path_unb_right_demod = "/home/nev/Documents/fastdiffamp/data-noise-set/unb-right-demod/"
path_unb_left_demod_comp = "/home/nev/Documents/fastdiffamp/data-noise-set/unb-left-demod-comp/"
path_unb_right_demod_comp = "/home/nev/Documents/fastdiffamp/data-noise-set/unb-right-demod-comp/"

def extract_data_from_csv_file(filename: str):
    file = open(filename)
    lines = file.readlines()
    if len(lines) < 2:
        raise Exception("Csv file error")
    elems = lines[1].split(",")
    start = float(elems[2])
    end = start + float(elems[3]) * (len(lines) - 3)
    y = []
    for i in range(2, len(lines)):
        y.append(float(lines[i].split(",")[1]))
    return [start, end, y]

def get_variance(data):
    n = len(data)
    mean = sum(data) / n
    variance = sum((x - mean) ** 2 for x in data) / (n - 1)
    return variance

def get_mean(data):
    n = len(data)
    mean = sum(data) / n
    return mean

def get_data_list(path_to_file, up_to):
    data_list = []
    for i in range(1, up_to):
        data_list.append(extract_data_from_csv_file(path_to_file + "NewFile" + str(i) + ".csv"))
    return data_list

#data_list5 = get_data_list(path_fifth, 12)
#data_list6 = get_data_list(path_sixth, 12)
#data_list7 = get_data_list(path_seventh, 12)
#data_list8 = get_data_list(path_eighth, 12)
#time = np.linspace(data_list[0][0], data_list[0][1], len(data_list[0][2]))
#data_weak_balanced = get_data_list(path_weak_balanced, 11)
#data_strong_balanced = get_data_list(path_strong_balanced, 11)
#data_strong_unb = get_data_list(path_strong_unb, 11)
#data_str_wk_balanced = get_data_list(path_str_wk_balanced, 11)
#data_str_wk_unb = get_data_list(path_str_wk_unb, 11)
#data_unb_left = get_data_list(path_unb_left, 12)
#data_unb_right = get_data_list(path_unb_right, 12)
#data_bal = get_data_list(path_bal, 12)
data_bal_demod = get_data_list(path_bal_demod, 19)
data_unb_left_demod = get_data_list(path_unb_left_demod_comp, 23)
data_unb_right_demod = get_data_list(path_unb_right_demod_comp, 24)

#for i in range(0, 11):
#    data_list.append(extract_data_from_csv_file(path_sixth + "NewFile" + str(i) + ".csv"))

#fig, ax = plt.subplots(1, 1)
#fig.set_size_inches(18.5, 10.5)
#var_list = []
#for data in data_list:
#    ax.plot(time, data[2])
#ax.legend()
#ax.set_xlabel('time')

def get_var_list(data_list):
    var_list = []
    for data in data_list:
        var_list.append(get_variance(data[2]))
    return var_list

def get_mean_list(data_list):
    mean_list = []
    for data in data_list:
        mean_list.append(get_mean(data[2]))
    return mean_list

#var_list5 = get_var_list(data_list5)
#var_list6 = get_var_list(data_list6)
#var_list7 = get_var_list(data_list7)
#var_list8 = get_var_list(data_list8)
#var_weak_balanced = get_var_list(data_weak_balanced)
#var_strong_balanced = get_var_list(data_strong_balanced)
#var_strong_unb = get_var_list(data_strong_unb)
#var_str_wk_balanced = get_var_list(data_str_wk_balanced)
#var_str_wk_unb = get_var_list(data_str_wk_unb)
#var_unb_left = get_var_list(data_unb_left)
#var_unb_right = get_var_list(data_unb_right)
#var_bal = get_var_list(data_bal)
var_bal_demod = get_var_list(data_bal_demod)
bal_mean = [-0.22, -0.20, -0.18, -0.17, -0.15, -0.3, -0.11, -0.09, -0.07, -0.07, -0.06, -0.06, -0.06, -0.06, -0.06, -0.06, -0.06, -0.06]
var_unb_left_demod = get_var_list(data_unb_left_demod)
unb_left_mean = [-7.81, -7.04, -6.23, -5.31, -4.51, -4.32, -4.12, -3.95, -3.76, -3.59, -2.67, -1.78, -0.91, -0.75, -0.58, -0.39, -0.31, -0.22, -0.14, -0.12, -0.10, -0.09]
var_unb_right_demod = get_var_list(data_unb_right_demod)
unb_right_mean = [7.61, 6.83, 6.01, 5.13, 4.69, 4.36, 4.01, 3.67, 3.28, 2.95, 2.59, 2.41, 1.65, 0.794, 0.625, 0.470, 0.279, 0.197, 0.109, 0.017, 0.003, -0.011, -0.020]

#mean_weak_balanced = get_mean_list(data_weak_balanced)
#mean_strong_balanced = get_mean_list(data_strong_balanced)
#mean_strong_unb = get_mean_list(data_strong_unb)
#mean_str_wk_balanced = get_mean_list(data_str_wk_balanced)
#mean_str_wk_unb = get_mean_list(data_str_wk_unb)

#power_data5 = [0.780, 0.750, 0.715, 0.590, 0.485, 0.367, 0.259, 0.153, 0.059, 0.026, 0.025] # for the 5th
#power_data6 = [0.0988, 0.0958, 0.089, 0.075, 0.0617, 0.046, 0.0328, 0.016, 0.0065, 0.0025, 0.0037] # for the 6th
#power_data56 = [0.780, 0.750, 0.715, 0.590, 0.485, 0.367, 0.259, 0.153, 0.059, 0.026, 0.025, 0.0988, 0.0958, 0.089, 0.075, 0.0617, 0.046, 0.0328, 0.016, 0.0065, 0.0025, 0.0037] # for the 5th and 6th
#power_data7 = [0.771, 0.765, 0.720, 0.606, 0.517, 0.354, 0.267, 0.125, 0.053, 0.007, 0.020] # for the 7th
#power_data8 = [0.780, 0.768, 0.725, 0.610, 0.512, 0.353, 0.210, 0.112, 0.035, 0.008, 0.019] # for the 8th
#power_weak_balanced = [0.016, 0.015, 0.0128, 0.0105, 0.0076, 0.0044, 0.0021, 0.0013, 0.00097, 0.00008]
# old str wk power_strong_balanced = [1.77, 1.66, 1.52, 1.32, 0.99, 0.69, 0.43, 0.22, 0.14, 0.12]
# old str wk power_strong_unb = [1.77, 1.69, 1.48, 1.31, 1.02, 0.70, 0.43, 0.26, 0.15, 0.12]
# old str wk power_str_wk_balanced = [0.664, 0.560, 0.465, 0.332, 0.200, 0.100, 0.060, 0.051]
# old str wk power_str_wk_unb = [0.692, 0.570, 0.442, 0.333, 0.195, 0.103, 0.064, 0.051]
#power_strong_balanced = [2.0, 1.9, 1.7, 1.45, 1.2, 0.8, 0.47, 0.288, 0.195, 0.130]
#power_strong_unb = [2.0, 1.9, 1.73, 1.45, 1.14, 0.86, 0.49, 0.295, 0.160, 0.131]
#power_strong_unb = [x / 2.0 for x in power_strong_unb]
#power_str_wk_balanced = [0.936, 0.905, 0.830, 0.700, 0.540, 0.370, 0.239, 0.126, 0.071, 0.058]
#power_str_wk_unb = [0.936, 0.936, 0.820, 0.670, 0.546, 0.380, 0.244, 0.135, 0.075, 0.061]
#power_unb_left = [0.009, 0.018, 0.037, 0.070, 0.111, 0.2, 0.330, 0.495, 0.600, 0.872, 0.925]
#power_unb_right = [0.009, 0.018, 0.042, 0.080, 0.113, 0.195, 0.323, 0.479, 0.648, 0.854, 0.915]
#power_bal = [0.009, 0.020, 0.044, 0.075, 0.110, 0.196, 0.311, 0.470, 0.625, 0.870, 0.919]
#power_str_wk_unb = [x / 2.0 for x in power_str_wk_unb]
power_bal = [0.890, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.08, 0.06, 0.04, 0.03, 0.02, 0.01, 0.008, 0.006, 0.0048]
power_unb_left = [0.890, 0.8, 0.7, 0.6, 0.5, 0.48, 0.46, 0.44, 0.42, 0.40, 0.3, 0.2, 0.1, 0.08, 0.06, 0.04, 0.03, 0.02, 0.01, 0.008, 0.006, 0.0048]
power_unb_right = [0.890, 0.8, 0.7, 0.6, 0.54, 0.5, 0.46, 0.42, 0.38, 0.34, 0.3, 0.28, 0.2, 0.1, 0.08, 0.06, 0.04, 0.03, 0.02, 0.01, 0.008, 0.006, 0.0048]

def get_curve_fit(x_points, y_points, order):
    z = np.polyfit(x_points, y_points, order)
    x = np.linspace(x_points[0], x_points[-1], 100)
    y_points = z[order].copy()
    for i in range(1, order + 1):
        y_points += x**i * z[order - i]
    return [x, y_points]

def get_exp_fit(x_points, y_points):
    z = np.polyfit(x_points, np.log(y_points), 2)
    x = np.linspace(x_points[0], x_points[-1], 100)
    return [x, np.exp(z[2]) * np.exp(z[0] * x**2)]

#fit5 = get_curve_fit(power_data5, var_list5)
#fit6 = get_curve_fit(power_data6, var_list6)
#fit7 = get_curve_fit(power_data7, var_list7)
#fit8 = get_curve_fit(power_data8, var_list8)
#fit_weak_balanced = get_curve_fit(power_weak_balanced, var_weak_balanced, 1)
#fit_strong_balanced = get_curve_fit(power_strong_balanced, var_strong_balanced, 1)
#fit_strong_unb_0 = get_curve_fit(power_strong_unb[0:5], var_strong_unb[0:5], 1)
#fit_strong_unb_1 = get_curve_fit(power_strong_unb[5:10], var_strong_unb[5:10], 1)
#fit_str_wk_balanced = get_curve_fit(power_str_wk_balanced, var_str_wk_balanced, 1)
#fit_str_wk_unb_0 = get_curve_fit(power_str_wk_unb[0:3], var_str_wk_unb[0:3], 1)
#fit_str_wk_unb_1 = get_curve_fit(power_str_wk_unb[2:10], var_str_wk_unb[2:10], 1)
fit_bal = get_curve_fit(power_bal, var_bal_demod, 1)
fit_unb_left = get_curve_fit(power_unb_left[10:22], var_unb_left_demod[10:22], 1)
fit_unb_right = get_curve_fit(power_unb_right[10:23], var_unb_right_demod[10:23], 1)

fig2, ax2 = plt.subplots(1, 1)
fig2.set_size_inches(18.5, 10.5)
#ax2.plot(power_data5, var_list5, ".", label="var 5")
#ax2.plot(fit5[0], fit5[1], label="fit 5")
#ax2.plot(power_data6, var_list6, ".", label="var 6")
#ax2.plot(fit6[0], fit6[1], label="fit 6")
#ax2.plot(power_data7, var_list7, ".", label="var 7")
#ax2.plot(fit7[0], fit7[1], label="fit 7")
#ax2.plot(power_data8, var_list8, ".", label="var 8")
#ax2.plot(fit8[0], fit8[1], label="fit 8")
#ax2.plot(power_weak_balanced, var_weak_balanced, "+", color="fuchsia", label="balanced detection (0.8 µW to 16.0 µW)")
#ax2.plot(fit_weak_balanced[0], fit_weak_balanced[1], color="pink", label="fit (0.8 µW to 16.0 µW)")
#ax2.plot(power_strong_balanced, var_strong_balanced, "+", color="blue", label="balanced detection (0.130 mW to 2.0 mW)")
#ax2.plot(fit_strong_balanced[0], fit_strong_balanced[1], color="cyan", label="fit (0.130 mW to 2.0 mW)")
#ax2.plot(power_str_wk_balanced, var_str_wk_balanced, "+", color="darkgreen", label="balanced detection (58.0 µW to 0.936 mW)")
#ax2.plot(fit_str_wk_balanced[0], fit_str_wk_balanced[1], color="forestgreen", label="fit (58 µw to 0.936 mW)")
#ax2.plot(power_strong_unb, var_strong_unb, "+", color="peru", label="unbalanced detection (0.131 mW to 2.0 mW)")
#ax2.plot(fit_strong_unb_0[0], fit_strong_unb_0[1], color="orange", label="fit (0.131 mW to 2.0 mW)")
#ax2.plot(fit_strong_unb_1[0], fit_strong_unb_1[1], color="orange")
#ax2.plot(power_str_wk_unb, var_str_wk_unb, "+", color="red", label="unbalanced detection (61 µW to 0.936 mW)")
#ax2.plot(fit_str_wk_unb_0[0], fit_str_wk_unb_0[1], color="tomato", label="fit (61 µW to 0.936 mW)")
#ax2.plot(fit_str_wk_unb_1[0], fit_str_wk_unb_1[1], color="tomato")
#ax2.plot(power_unb_left, var_unb_left, ".", color="black", label="unbalanced detection left (9 µW to 0.925 mW)")
#ax2.plot(power_unb_right, var_unb_right, ".", color="blue", label="unbalanced detection right (9 µW to 0.915 mW)")
#ax2.plot(power_bal, var_bal, ".", color="red", label="balanced detection (9 µW to 0.919 mW)")
ax2.plot(power_bal, var_bal_demod, "+", label="balanced detection")
ax2.plot(power_unb_left, var_unb_left_demod, "+", label="unbalanced detection (left diode ON)")
ax2.plot(power_unb_right, var_unb_right_demod, "+", label="unbalanced detection (right diode ON)")
ax2.plot(fit_bal[0], fit_bal[1], label="fit balanced detection")
ax2.plot(fit_unb_left[0], fit_unb_left[1], label="fit unbalanced detection (left diode ON)")
ax2.plot(fit_unb_right[0], fit_unb_right[1], label="fit unbalanced detection (right diode ON)")
ax2.legend()
ax2.set_xlabel('power (mW)')
ax2.set_ylabel('Variance of the output tension (V)')

#fit_mean_weak_balanced = get_curve_fit(power_weak_balanced, mean_weak_balanced, 1)
#fit_mean_strong_balanced = get_curve_fit(power_strong_balanced, mean_strong_balanced, 1)
#fit_mean_strong_unb = get_curve_fit(power_strong_unb, mean_strong_unb, 1)
#fit_mean_str_wk_balanced = get_curve_fit(power_str_wk_balanced, mean_str_wk_balanced, 1)
#fit_mean_str_wk_unb_0 = get_curve_fit(power_str_wk_unb[0:3], mean_str_wk_unb[0:3], 1)
#fit_mean_str_wk_unb_1 = get_curve_fit(power_str_wk_unb[2:10], mean_str_wk_unb[2:10], 1)
#
fig3, ax3 = plt.subplots(1, 1)
fig3.set_size_inches(18.5, 10.5)
#ax3.plot(power_weak_balanced, mean_weak_balanced, "+", color="fuchsia", label="balanced detection (0.8 µW to 16.0 µW)")
#ax3.plot(fit_weak_balanced[0], fit_mean_weak_balanced[1], color="pink", label="fit (0.8 µW to 16.0 µW)")
#ax3.plot(power_strong_balanced, mean_strong_balanced, "+", color="blue", label="balanced detection (0.130 mW to 2.0 mW)")
#ax3.plot(fit_strong_balanced[0], fit_mean_strong_balanced[1], color="cyan", label="fit (0.130 mW to 2.0 mW)")
#ax3.plot(power_str_wk_balanced, mean_str_wk_balanced, "+", color="darkgreen", label="balanced detection (58.0 µW to 0.936 mW)")
#ax3.plot(fit_mean_str_wk_balanced[0], fit_mean_str_wk_balanced[1], color="forestgreen", label="fit (58 µw to 0.936 mW)")
#ax3.plot(power_strong_unb, mean_strong_unb, "+", color="peru", label="unbalanced detection (0.131 mW to 2.0 mW)")
#ax3.plot(fit_mean_strong_unb[0], fit_mean_strong_unb[1], color="orange", label="fit (0.131 mW to 2.0 mW)")
#ax3.plot(power_str_wk_unb, mean_str_wk_unb, "+", color="red", label="unbalanced detection (61 µW to 0.936 mW)")
#ax3.plot(fit_mean_str_wk_unb_0[0], fit_mean_str_wk_unb_0[1], color="tomato", label="fit (61 µW to 0.936 mW)")
#ax3.plot(fit_mean_str_wk_unb_1[0], fit_mean_str_wk_unb_1[1], color="tomato")
ax3.plot(power_bal, bal_mean, "+", label="balanced")
ax3.plot(power_unb_left, unb_left_mean, "+", label="unbalanced left")
ax3.plot(power_unb_right, unb_right_mean, "+", label="unbalanced right")
ax3.legend()
ax3.set_xlabel('power (mW)')
ax3.set_ylabel('Mean output tension (V)')

plt.show()
exit()

pw_dt_cpy = power_data.copy()
varlist_fst = []
varlist_sec = []
varlist_thd = []
for i in range(0, len(var_list), 3):
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
