import numpy as np
import matplotlib.pyplot as plt
from da1 import da1
from da2_snake import da2_snake

##
# read multiple data files
lx = 1.0
ly = 1.0

data_path = r'D:\numerical_simulation\share\mcl_complex\tt\\'

m = 128
n = 128

# Use da2_snake to process and plot the data sequence
x_lst, y_lst, mass_lst, Area_lst, t_lst, u_c_lst = da2_snake(
    0, 1, 3, m, n, data_path, '', lx, ly, False
)


##
# read single data files
lx = 1.0
ly = 1.0

data_path = r'D:\numerical_simulation\share\mcl_complex\tt\\'

m = 128
n = 128
r, z, u, v, p, f, f2, t = da1(data_path + 'data003', n + 1, m + 1) # , m + 1, n + 1)

# Plot
plt.figure()
plt.title(f't = {t}')
plt.contour(r, z, f, levels=[0.5], colors='m')
plt.contour(r, z, f2, levels=[0.5], colors='m')
plt.axis('equal')
plt.title(f't = {t}')

ip = 5
scale_factor = 0.1  # Adjust arrow size

# Downsample for quiver plot
r_q = r[::ip, ::ip]
z_q = z[::ip, ::ip]
u_q = u[::ip, ::ip] * scale_factor
v_q = v[::ip, ::ip] * scale_factor

plt.quiver(r_q, z_q, u_q, v_q, color='r', angles='xy', scale_units='xy', scale=1)
#plt.quiver(r, z, u, v)
plt.show()
