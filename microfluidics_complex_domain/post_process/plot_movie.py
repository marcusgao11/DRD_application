import numpy as np
import matplotlib.pyplot as plt
from da1 import da1
from da2_snake import da2_snake

# --- Part 1: snake case ---
lx = 3.5
ly = 1.5
data_path = r'./'
m = 420
n = 180

# Placeholder for da2_snake function
x_lst, y_lst, mass_lst, Area_lst, t_lst, u_c_lst = da2_snake(0, 10, 200, m, n, data_path, '', lx, ly, False)