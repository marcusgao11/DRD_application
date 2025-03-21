
% case 2: snake
clc; clear; close all;

lx = 2.0; ly = 1.0;

data_path = 'E:\numerical_simulation\share\mcl_complex_fiber_modify_u_1.0_largeB_degenerateM\data\';

m=420;n=180;
[x_lst, y_lst, mass_lst, Area_lst, t_lst, u_c_lst] = da2_snake(0,1,12, m,n, data_path,'', lx,ly,false);



    