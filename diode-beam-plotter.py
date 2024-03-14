from matplotlib import pyplot as plt
import numpy as np


position_list = np.linspace(88, 110, 23)
print(position_list)
position_list /= 10.0
position_list -= 10.0
print(position_list)

left_diode = [-0.726, -0.732, -0.750, -0.836, -1.29, -4.08, -8.85, -10.9, -10.6, -10.7, -10.7, -10.6, -11.1, -10.6, -10.9, -10.7, -10.7, -9.80, -4.79, -1.32, -0.741, -0.735, -0.717]
right_diode = [-0.533, -0.407, 0.982, 0.736, 1.20, 9.02, 10.7, 10.8, 11.0, 11.1, 11.1, 11.1, 10.7, 10.6, 10.5, 10.2, 7.04, 1.56, -0.100, -0.492, -0.580, -0.620, -0.640]

left_diode_sec = [-0.282, -0.266, -0.605, -1.55, -3.14, -4.37, -5.09, -5.26, -5.27, -5.44, -5.27, -5.37, -5.38, -5.31, -4.99, -4.12, -2.74, -1.32, -0.515, -0.265, -0.173, -0.145, -0.186]
right_diode_sec = [-0.106, -0.331, 0, 0.524, 1.57, 3.16, 4.32, 4.75, 4.78, 4.84, 4.81, 4.81, 4.77, 4.74, 4.47, 3.53, 2.11, 0.924, 0.392, 0.100, 0, -0.106, -0.097]

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(18.5, 10.5)
ax.plot(position_list, left_diode, color = 'red', label = "left diode (2.3 mW input power on main beam)")
ax.plot(position_list, right_diode, linestyle='dashed', color = 'red', label = "right diode (2.3 mW input power on main beam)")
ax.plot(position_list, left_diode_sec, color = 'orange', label = "left diode (1.25 mW input power on main beam)")
ax.plot(position_list, right_diode_sec, linestyle='dashed', color = 'orange', label = "right diode (1.25 mW input power on main beam)")
ax.legend()
ax.grid()
ax.set_xlabel('relative position diode/beam (mm)')
ax.set_ylabel('average volt output signal (V)')

plt.show()
