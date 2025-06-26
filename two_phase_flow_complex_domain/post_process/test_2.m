

%%
% read multiple data files
clc; clear; close all;

lx = 1.0; ly = 1.0;

data_path = 'D:\numerical_simulation\share\mcl_complex\tt\';

m=128;n=128;
[x_lst, y_lst, mass_lst, Area_lst, t_lst, u_c_lst] = da2_snake(0,1,12, m,n, data_path,'', lx,ly,false);



%%
% read single data file
clc; clear; close all;

lx = 1.0; ly = 1.0;
data_path = 'D:\numerical_simulation\share\mcl_complex\tt\';

m=128;n=128;

[r,z,u,v,p,f,f2,t]=da1([data_path 'data003'],m+1,n+1);
% plot
figure; hold on;
title(['t=' num2str(t)]);
contour(r,z,f,[0.5 0.5], 'm'); hold on;
contour(r,z,f2,[0.5 0.5], 'm');
axis equal;
title(['t = ' num2str(t)])

ip = 5;
scale_factor = 0.1; % 缩放系数，调整箭头的大小
quiver(r(1:ip:end,1:ip:end),z(1:ip:end,1:ip:end),...
    u(1:ip:end,1:ip:end)*scale_factor,v(1:ip:end,1:ip:end)*scale_factor,'r','AutoScale', 'off');
%quiver(r,z,u,v);