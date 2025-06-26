import numpy as np
import matplotlib.pyplot as plt
import os
from da1 import da1

def da2_snake(N_start, dN, N_end, M, N, data_path, indicator, lx, ly, is_close=True):
    is_plot = True

    if is_close:
        plt.close('all')

    x_lst = []
    y_lst = []
    Area_lst = []
    mass_lst = []
    t_lst = []
    u_c_lst = []

    # Read boundary phase
    r, z, u, v, p, f, f2, t = da1(data_path + 'databdr', N + 1, M + 1) #, M + 1, N + 1)
    bdr_ph = f

    # Prepare video writer if needed (optional, not implemented here)
    # import cv2 or imageio for video writing if required

    ch = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    for i in range(N_start, N_end + 1, dN):
        k3 = i // 100
        k2 = (i - k3 * 100) // 10
        k1 = i % 10

        fname = f"data{ch[k3]}{ch[k2]}{ch[k1]}"
        if indicator:
            fname = f"{fname}{indicator}"
        print(os.path.join(data_path, fname))

        r, z, u, v, p, f, f2, t = da1(os.path.join(data_path, fname), N + 1, M + 1) #, M + 1, N + 1)
        ff = f - bdr_ph

        dx = r[0, 1] - r[0, 0]
        dy = z[1, 0] - z[0, 0]
        mass = np.sum(f2) * dx * dy
        mass_lst.append(mass)
        t_lst.append(t)

        if is_plot:
            plt.figure(1)
            plt.clf()
            plt.title(f't = {t}')
            plt.contour(r, z, f, levels=[0.5], colors='m')
            plt.contour(r, z, f2, levels=[0.5], colors='m')
            plt.axis('equal')

            ip = 5
            scale_factor = 0.1
            plt.quiver(
                r[::ip, ::ip], z[::ip, ::ip],
                u[::ip, ::ip] * scale_factor, v[::ip, ::ip] * scale_factor,
                color='r', angles='xy', scale_units='xy', scale=1
            )
            plt.pause(0.5)

    # Video writing not implemented; use matplotlib.animation or imageio if needed

    return np.array(x_lst), np.array(y_lst), np.array(mass_lst), np.array(Area_lst),np.array(t_lst),np.array(u_c_lst)