from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import j0
from scipy.special import j1

amplitude = np.linspace(25, 100, num=16)
ratio = [0.23, 0.41, 0.46, 1.05, 1.67, 2.67, 3.67, 5.06, 6.16, 8.31, 10, 12.37, 14.19, 17.16, 18.50, 21.16]

amplitude_ext = np.linspace(25, 105, num=17)
ratio_ext = [0.23, 0.41, 0.46, 1.05, 1.67, 2.67, 3.67, 5.06, 6.16, 8.31, 10, 12.37, 14.19, 17.16, 18.50, 21.16, 24]

# Exponential approximation
z = np.polyfit(amplitude_ext, np.log(ratio_ext), 3)
def fit_exp(x):
    return np.exp(z[3]) * np.exp(z[2] * x) * np.exp(z[1] * np.power(x, 2)) * np.exp(z[0] * np.power(x, 3))
xfit = np.linspace(0, 100, 1000)
yfit = fit_exp(xfit) #z[2] + z[1] * xfit + z[0] * xfit**2

def theo_ratio(x, a, b):
    return a * (j1(x / 100.0 * b)/j0(x / 100.0 * b))**2
popt, pcov = curve_fit(theo_ratio, amplitude[0:len(amplitude) - 9], ratio[0:len(ratio) - 9])
x = np.linspace(0, 60, num=1000)

def alpha(x):
    return x / (1.0 + 2.0 * x)

print(yfit)
alpha_list = alpha(yfit / 10.0) # one divides by 10 because yfit is in percent

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(18.5, 10.5)
ax.plot(amplitude, ratio, "+", label="experimental data")
ax.plot(xfit, yfit, label="exponential fit (3rd degree)")
print(popt[0], " ", popt[1])
ax.plot(x, theo_ratio(x, popt[0], popt[1]), label="bessel fit")
ax2 = ax.twinx()
ax2.plot(xfit, alpha_list, label="alpha")
ax2.set_ylabel('Alpha coefficient: Ps = α * Ptot')
ax.set_xlabel('Amplitude (%)')
ax.set_ylabel('Sidebands ratio (%)')
ax.legend()
ax2.legend()
ax.grid()

fig3, ax3 = plt.subplots(1, 1)
fig3.set_size_inches(18.5, 10.5)
amplitude_by_ten = np.array([10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100])
dBm_diff = [1, 4, 8, 14, 18, 28, 39, 41, 42, 45, 45, 45]
ax3.plot(amplitude_by_ten, dBm_diff, "+", color="blue")
ax3.set_ylabel("Power difference (43 MHz beat - noise) (dB)", color="blue")
ax4 = ax3.twinx()
sideband_intensity = 10.0 * alpha(fit_exp(amplitude_by_ten) / 10)
ax4.plot(amplitude_by_ten, sideband_intensity, "+", color="orange")
ax4.set_ylabel("Sideband intensity (nW)", color="orange")
ax3.set_xlabel('Amplitude (%)')
ax3.grid()
ax3.set_title("Amplitude variation at minimum total power (0.01µW from the cavity)")

fig5, ax5 = plt.subplots(1, 1)
fig5.set_size_inches(18.5, 10.5)
amplitude_by_ten = np.array([10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100])
sideband_power = 10.0 * alpha(fit_exp(amplitude_by_ten) / 10)
dBm_diff = [1, 4, 8, 14, 18, 28, 39, 41, 42, 45, 45, 45]
ax5.plot(sideband_power, dBm_diff, "+")
ax5.set_ylabel("Power difference (43 MHz beat - noise) (dB)", color="blue")
ax5.set_xlabel('Sideband Power (nW)')
ax5.set_xscale('log')
ax5.grid()
ax5.set_title("Amplitude variation at minimum total power (0.01µW from the cavity)")

plt.show()