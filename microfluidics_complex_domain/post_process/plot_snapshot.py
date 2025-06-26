import numpy as np
import matplotlib.pyplot as plt
from da1 import da1
from da2_snake import da2_snake


# --- Part 2: read single data file ---
lx = 3.5
ly = 1.5
data_path = r'./'
m = 420
n = 180

r, z, u, v, p, f, f2, t = da1(data_path + 'data050', n + 1, m + 1) #m + 1, n + 1)


# Plot
plt.figure()
plt.title(f't = {t}')
plt.contour(r, z, f, levels=[0.5], colors='m')
plt.contour(r, z, f2, levels=[0.5], colors='m')
plt.title(f't = {t}')
plt.axis('equal')
plt.axis([0, lx, 0, ly])

ip = 5
scale_factor = 0.1  # Adjust arrow size

# Downsample for quiver plot
r_q = r[::ip, ::ip]
z_q = z[::ip, ::ip]
u_q = u[::ip, ::ip] * scale_factor
v_q = v[::ip, ::ip] * scale_factor

plt.quiver(r_q, z_q, u_q, v_q, color='r', angles='xy', scale_units='xy', scale=1)
plt.savefig('snake_complex.png', dpi=300, bbox_inches='tight')
#plt.quiver(r, z, u, v)
plt.show()
