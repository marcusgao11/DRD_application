

function [x_lst, y_lst, mass_lst, Area_lst, t_lst, u_c_lst] = da2_snake(N_start,dN,N_end,M,N,data_path, indicator, lx,ly, is_close)


is_plot = true;

if nargin < 8 || isempty(is_close)
    is_close = true; % 设置默认值
end
    
if is_close
    close all;
end

x_lst = []; y_lst = []; Area_lst = []; mass_lst = [];
t_lst=[];
u_c_lst = [];

[r,z,u,v,p,f,t]=da1([data_path 'databdr'],M+1,N+1);
bdr_ph = f;


outputVideo = VideoWriter('SNAKE_2_20240622.avi');
outputVideo.FrameRate = 0.8;  % 设置帧率
% 打开文件进行写入
open(outputVideo);

        
ch = ['0' '1' '2' '3' '4' '5' '6' '7' '8' '9'];
for i=N_start:dN:N_end
       %k1=floor(i/10);
       %k2=i-k1*10;
    
    k3 = floor(i/100);
    k2 = floor((i-k3*100)/10);
    k1 = mod(i,10);
    
    fname = ['data' ch(k3+1) ch(k2+1) ch(k1+1)];
    if ~isempty(indicator)
        fname = [fname indicator];
    end
    disp([data_path fname]);
    
    [r,z,u,v,p,f,f2,t]=da1([data_path fname],M+1,N+1);
   
    ff = f - bdr_ph;    
    
    try
        [Area,Centroid,IN]=Contour2Area(C);
        Area
        Centroid  
        x_lst = [x_lst; Centroid(1)]; Area_lst = [Area_lst; Area];
        y_lst = [y_lst; Centroid(2)]; 

    cathch ME
        disp(['发生了一个错误: ', ME.message]);
    end
    
    dx = r(1,2)-r(1,1); dy=z(2,1)-z(1,1);
    mass_lst = [mass_lst, sum(sum(f2))*dx*dy];

    if is_plot
        figure(1); hold off;  title(['t=' num2str(t)]);
        contour(r,z,f,[0.5 0.5], 'm'); hold on;
        contour(r,z,f2,[0.5 0.5], 'm');
        axis equal;
        
        hold on; %,contour(x,y,r');
        %axis([0 Lx 0 Ly]); axis equal;
        ip = 5;

        scale_factor = 0.1; % 缩放系数，调整箭头的大小
        quiver(r(1:ip:end,1:ip:end),z(1:ip:end,1:ip:end),...
            u(1:ip:end,1:ip:end)*scale_factor,v(1:ip:end,1:ip:end)*scale_factor,'r','AutoScale', 'off');
        %quiver(r,z,u,v);

        hold off;
        
        % 捕捉当前帧并写入视频
        frame = getframe(gcf);
        writeVideo(outputVideo, frame);
        
    end
   
    
    i,t
    pause(0.5);
    
end

close(outputVideo);

end

