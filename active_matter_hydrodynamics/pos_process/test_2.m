



clc; clear; close all;
m=200;n=200;

data_path = 'E:\numerical_simulation\share\swimmer\';

quiver_params.scale  = 0.02; quiver_params.on_off = true;

quiver_params.is_plot = true;
quiver_params.is_write_vedio = true;
quiver_params.movie_nm = 'moive_0322';

quiver_params.is_pause = false;

[x_lst, y_lst, Area_lst, t_lst, u_c_lst, r,z,u,v] = da2_new(0,1,20, m,n, data_path,'', quiver_params);

aa = [t_lst, x_lst, y_lst]; 
save swm_fluid_particle.txt -ascii aa

disp('finished');

