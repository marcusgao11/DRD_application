import numpy as np
import os
import datetime
import struct


def da1(file, m, n):
    # 显示文件最后修改时间（类似 MATLAB 中的 ff_info.date）
    mod_time = os.path.getmtime(file)
    print("Last modified:", datetime.datetime.fromtimestamp(mod_time))

    with open(file, 'rb') as fid:
        fid.seek(4, 1)  # Skip 4 bytes
        r = np.fromfile(fid, dtype=np.float64, count=m * n).reshape((n, m))

        fid.seek(8, 1)
        z = np.fromfile(fid, dtype=np.float64, count=m * n).reshape((n, m))

        fid.seek(8, 1)
        u = np.fromfile(fid, dtype=np.float64, count=m * n).reshape((n, m))

        fid.seek(8, 1)
        v = np.fromfile(fid, dtype=np.float64, count=m * n).reshape((n, m))

        fid.seek(8, 1)
        p = np.fromfile(fid, dtype=np.float64, count=m * n).reshape((n, m))

        fid.seek(8, 1)
        f = np.fromfile(fid, dtype=np.float64, count=m * n).reshape((n, m))

        fid.seek(8, 1)
        f2 = np.fromfile(fid, dtype=np.float64, count=m * n).reshape((n, m))

        fid.seek(8, 1)
        t = struct.unpack('d', fid.read(8))[0]  # Read a single double

    return r, z, u, v, p, f, f2, t