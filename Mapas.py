import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse
# This program calculates the tangential and sagittal curvature maps. 

plt.rc('text', usetex=True)
diop = 333  # curvature to diopter  (n-1)*C
N = 60  # pupil
Div = 150
dr = N / Div
N2 = 2 * np.pi
dt = N2 / 150
r, tt = np.mgrid[slice(-N, N + dr, dr),
                 slice(dt, N2 + dt, dt)]
x = r * np.cos(tt)
y = r * np.sin(tt)
c = 1/N     #np.divide(1, N, where=r != 0)  # 1/r
s2 = x**2 + y**2
ro2 = x**2 + y**2
# Slopes for curvature calculation
factor = np.sqrt(1 - (c ** 2) * s2)
factor[np.isnan(factor)] = 0
print(factor)
fx = np.divide(c * x, factor, where=factor != 0)
fy = np.divide(c * y, factor, where=factor != 0)
fxx = np.divide(c * (1 - (c ** 2) * (y ** 2)), factor ** (3), where=factor ** 3 != 0)
fyy = np.divide(c * (1 - (c ** 2) * (x ** 2)), factor ** (3), where=factor ** 3 != 0)
fxy = np.divide( (c ** 2) * (x * y), factor ** (3), where=factor ** 3 != 0)
f_x = fx
f_xx = fxx
f_y = fy
f_yy = fyy
f_xy = fxy
dig = 4
# Tangential
ctagNum = 0.5 * (f_xx + f_yy) + 0.5 * (f_xx - f_yy) * ((x ** 2 - y ** 2) * (1/ro2)) + 2 * f_xy * (x * y * (1/ro2))
ctagNum[np.isnan(ctagNum)] = 0
ctagDen = ((1 + (1/ro2)*(f_x * x + f_y * y) ** 2) ** (3 / 2)) * (
            (1 + (1/ro2)*(f_x * y - f_y * x) ** 2) ** (1 / 2))
ctagDen[np.isnan(ctagDen)] = 0
ctag = ctagNum / ctagDen
ctag[np.isnan(ctag)] = 0
# Sagittal
csagNum = 0.5 * (f_xx + f_yy) - 0.5 * (f_xx - f_yy) * ((x ** 2 - y ** 2) * (1/ro2) ) - 2 * f_xy * (x * y * (1/ro2) )
csagDen = ((1 + (1/ro2)*(f_x * y  - f_y * x ) ** 2) ** (3 / 2)) * ((1 + (1/ro2)*(f_x * x  + f_y * y ) ** 2) ** (1 / 2))

csag = csagNum / csagDen

z = ctag

z_min, z_max = np.min(z), np.max(z)

fig, ax = plt.subplots(1, 1)

c = ax.pcolormesh(x, y, z, cmap='jet', vmin=z_min, vmax=z_max, shading='nearest')

fig.colorbar(c, ax=ax)
#fig.suptitle(r'$Curvature \ tangential$', fontsize=20)
#ax.tick_params(labelsize=16, width=2)
#ax.grid(True, linestyle='-.')
#ax.set_ylabel(r'$Pupil \ (mm)$', fontsize = 16)

plt.show()
