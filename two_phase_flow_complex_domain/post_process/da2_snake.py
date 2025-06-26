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

        r, z, u, v, p, f, f2, t = da1(os.path.join(data_path, fname), M + 1, N + 1)
        ff = f - bdr_ph

        # Placeholder for contour-to-area/centroid calculation
        # You need to implement Contour2Area or use skimage.measure
        try:
            from skimage import measure

            # Find contours at level 0.5
            C = measure.find_contours(ff, 0.5)
            if C:
                # Take the largest contour
                contour = max(C, key=lambda x: x.shape[0])
                # Convert pixel indices to coordinates
                contour_x = np.interp(contour[:, 1], np.arange(ff.shape[1]), r[0, :])
                contour_y = np.interp(contour[:, 0], np.arange(ff.shape[0]), z[:, 0])
                # Area using shoelace formula
                area = 0.5 * np.abs(np.dot(contour_x, np.roll(contour_y, 1)) - np.dot(contour_y, np.roll(contour_x, 1)))
                centroid_x = np.mean(contour_x)
                centroid_y = np.mean(contour_y)
                x_lst.append(centroid_x)
                y_lst.append(centroid_y)
                Area_lst.append(area)
            else:
                area = np.nan
                centroid_x = np.nan
                centroid_y = np.nan
        except Exception as e:
            print(f"Error in contour/area calculation: {e}")

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
