import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from occultation import *


adr = 'C:/Moi/_py/Astronomy/Solar System/kernels/'
files = ['naif0012.tls', 'pck00010.tpc', 'earth_latest_high_prec.bpc', 'de440_2030.bsp']
kernels = [adr+i for i in files]
bodies = [301, 4]


obs_loc = (7.744083817548831, 48.58313582900411, 140)
t_ini = datetime(2022, 12, 7)
t0 = find_exact_t0(t_ini, bodies, obs_loc, kernels)

time_window, ind_im, ind_em, pos_moon, pos_mars, size_moon, size_mars = \
             get_occultation(bodies, obs_loc, t0, kernels)

t_im = time_window[ind_im]
t_em = time_window[ind_em]

az_mars = [i[1] for i in pos_mars]
alt_mars = [i[2] for i in pos_mars]
az_moon = [i[1] for i in pos_moon]
alt_moon = [i[2] for i in pos_moon]


x = np.array(az_mars) - np.array(az_moon)
y = np.array(alt_mars) - np.array(alt_moon)

time_str = [str(i)[11:16] for i in time_window]

fig, ax = plt.subplots()
ax.scatter([0], [0], s=20, c='b')
moon_circle = plt.Circle((0, 0), size_moon.mean()/2, color='blue')
for i in range(len(x)):
    mars_circle = plt.Circle((x[i], y[i]), size_mars[i]/2, color='red')
    ax.add_patch(mars_circle)
ax.add_patch(moon_circle)
ax.axvline(x=x[ind_im], c='green', ls='-.')
ax.axvline(x=x[ind_em], c='brown', ls='-.')
ax.set_aspect('equal')
ax.set_xticks([x[0], x[ind_im], x[ind_em], x[-1]])
ax.set_xticklabels([time_str[0], time_str[ind_im], time_str[ind_em], time_str[-1]])
ax.set_xlabel('UTC Time (HH:MM)')
ax.set_ylabel('Degree')
plt.show()

# p = np.array(pos_mars)
# az, alt = p[:,1], p[:,2]
