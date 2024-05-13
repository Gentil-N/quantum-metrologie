from matplotlib import pyplot as plt
import numpy as np

def airy(phi, r1, r2):
    return (1-r1**2)*(1-r2**2)*(1.0 / ((1 - r1*r2)**2 + 4*r1*r2*np.sin(phi)**2))

def airy_dist(phi, r1, r2):
    return  airy(phi, r1, r2) / airy(0, r1, r2)

def ax_plot(x, r1, r2, label=""):
    ax1.plot(x, airy_dist(np.pi*(x+1), r1, r2), label=label)

def linewidth(r1, r2):
    return 2.0 / np.pi * np.arcsin((1-r1*r2)/(2*np.sqrt(r1*r2)))

def linewidthsec(r12):
    return 2.0 / np.pi * np.arcsin((1-r12)/(2*np.sqrt(r12)))

x = np.linspace(-5, 5, 1000)

fig1 = plt.figure(num=1)
ax1 = fig1.subplots(nrows=1, ncols=1)
ax_plot(x, 0.9, 0.9, label="r1=r2=0.9")
ax_plot(x, 0.8, 0.8, label="r1=r2=0.8")
ax_plot(x, 0.7, 0.7, label="r1=r2=0.7")
ax1.set(xlabel="Ф [rads]", ylabel="A'(Ф) normalized [no units]")
ax1.legend()

fig2 = plt.figure(num=2)
ax2 = fig2.subplots(nrows=1, ncols=1)
x = np.linspace(0, 1, 1000)
ax2.plot(x, linewidthsec(x))
ax2.set(xlabel="r1r2", ylabel="FMWH")
ax2.legend()

print(linewidth(0.9, 0.9))
print(linewidth(0.7, 0.7))
print(linewidth(0.74, 0.74))
print(linewidth(0.4, 0.4))

plt.show()
