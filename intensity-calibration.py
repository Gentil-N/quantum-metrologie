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
xfit = np.linspace(0, 100, 1000)
yfit = np.exp(z[3]) * np.exp(z[2] * xfit) * np.exp(z[1] * xfit**2) * np.exp(z[0] * xfit**3) #z[2] + z[1] * xfit + z[0] * xfit**2

def theo_ratio(x, a, b):
    return a * (j1(x / 100.0 * b)/j0(x / 100.0 * b))**2
popt, pcov = curve_fit(theo_ratio, amplitude, ratio)
x = np.linspace(0, 100, num=1000)

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(18.5, 10.5)
ax.plot(amplitude, ratio, "+", label="experimental data")
ax.plot(xfit, yfit, label="exponential fit (3rd degree)")
print(popt[0], " ", popt[1])
ax.plot(x, theo_ratio(x, popt[0], popt[1]), label="bessel fit")
ax.set_xlabel('Maximum  generator operation (%)')
ax.set_ylabel('Sidebands ratio (%)')
ax.legend()
ax.grid()
plt.show()