%%plot the result as time evolves

%m=256; n =128;
%da2(1,1,100,m,n, 'E:\numerical_simulation\tmp\', '')

% da2(1,1,100,m,n, 'E:\numerical simulation\sav_q\','');

% da2(1,1,100,m,n, 'E:\numerical simulation\share\','');

% m=256; n =512;
% da2(1,1,100,m,n, 'E:\numerical simulation\share\','');

function [x_lst, y_lst, mass_lst, Area_lst, t_lst, u_c_lst] = da2_snake(N_start,dN,N_end,M,N,data_path, indicator, lx,ly, is_close)


is_plot = true;

if nargin < 8 || isempty(is_close)
    is_close = true; % 设置默认值
end
    
if is_close
    close all;
end
% M = 257
% N = 513

x_lst = []; y_lst = []; Area_lst = []; mass_lst = [];
t_lst=[];
u_c_lst = [];

[r,z,u,v,p,f,t]=da1([data_path 'databdr'],M+1,N+1);
bdr_ph = f;

figure(3); hold on;
contour(r,z,bdr_ph,[0.1 0.5 0.9], 'k'); 
        

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
    
%     pcolor(r,z,f);
%     axis equal
%     shading interp;
    
        %contour(r,z,f,[1.0E-10 1.0E-10], 'm');
    %
    %
    

    ff = f - bdr_ph;    
    %C = contour(r,z,f2,[0.5 0.5], 'm');    
    
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
%     dx = geom.dx; dy = geom.dy;
     mass_lst = [mass_lst, sum(sum(f2))*dx*dy];

    
    if is_plot
        figure(1); hold off;  title(['t=' num2str(t)]);
        contour(r,z,f,[0.5 0.5], 'm'); hold on;
        contour(r,z,f2,[0.5 0.5], 'm');
        axis equal;
        title(['t = ' num2str(t)])

%         % analytical results
%         rad = 0.25*0.95; %0.25; 
%         theta = 45*pi/180; xc = lx/2.0; yc = ly/4.0;
%         Ana_radius = rad * sqrt( pi/2.0 / (theta - sin(theta)*cos(theta)) );
%         Ana_H = Ana_radius*(1 - cos(theta));
%         Ana_radius*sin(theta)
% 
%         Ana_xc = xc;
%         Ana_yc = yc - (Ana_radius - Ana_H);
% 
%         rectangle('Position', [Ana_xc - Ana_radius, Ana_yc - Ana_radius, 2*Ana_radius, 2*Ana_radius], ...
%               'Curvature', [1, 1], ...
%               'EdgeColor', 'b', 'LineWidth', 2);
%           axis equal;
          

%           A0 = pi*rad^2/2
%           ep = 0.005;
%           del_A = 5*ep*rad*2
%           del_A/A0
%           sqrt( (A0*0.9842)*2/pi )
          
        hold on; %,contour(x,y,r');
        %axis([0 Lx 0 Ly]); axis equal;
        ip = 5;

        scale_factor = 0.1; % 缩放系数，调整箭头的大小
        quiver(r(1:ip:end,1:ip:end),z(1:ip:end,1:ip:end),...
            u(1:ip:end,1:ip:end)*scale_factor,v(1:ip:end,1:ip:end)*scale_factor,'r','AutoScale', 'off');
        %quiver(r,z,u,v);

        hold off;
        
        figure(2); hold off;  title(['t=' num2str(t)]);
        contour(r,z,f2); hold on;
        contour(r,z,f,[0.5 0.5], 'm');
        %contour(r,z,bdr_ph,[0.1 0.5 0.9], 'm'); 
        
        axis equal;
        
        figure(3); hold on;  title(['t=' num2str(t)]);
        contour(r,z,f2,[0.5 0.5], 'm');
        %contour(r,z,bdr_ph,[0.1 0.5 0.9], 'k'); 
     
        figure(4); hold off;
        mesh(r,z,u); title(['t=' num2str(t)]);
        
        
%         figure(5); hold off;  title(['t=' num2str(t)]);
%         %f_tol = f2*4 + bdr_ph - (f - bdr_ph);
%         f_tol = bdr_ph + (1-bdr_ph)*2.0 + f2*2.0; hold on;
%         pcolor(r,z,f_tol); shading interp; %colormap(othercolor('parula'));  %colormap([jet(128); 1 0 0]);
%         colorbar;
%         
%         contour(r,z, f - bdr_ph, [0.5 0.5], 'r');
%         contour(r,z,bdr_ph, [0.5 0.5], 'm', 'LineWidth', 2);
        
        
        figure(5); hold off;
        
        f_tol = bdr_ph + (1-bdr_ph)*2.0 + f2*2.0;
        
        set(gcf,  'OuterPosition',[602.333333333333 266 1515.66666666667 892]);
                
        % 创建 contour
        contourf(r,z, f_tol,'LineWidth',2.5,'EdgeColor',[0 0 0],...
        'FaceColor',[0.301960784313725 0.745098039215686 0.933333333333333],...
        'LevelList',3);
   
        % 创建 ylabel
        ylabel('$z/H$','Interpreter','latex');

        % 创建 xlabel
        xlabel('$x/H$','Interpreter','latex');
        
        set(gca, 'FontSize',45);
        
        
        hold on;
        % 创建 contour
        contourf(r,z,f-bdr_ph, 'LineWidth',1,'EdgeColor',[0 0 0],...
            'FaceColor',[1 0.411764705882353 0.16078431372549],...
                    'LevelList',0.5);

        ip = 10; scale_factor = 0.1; % 缩放系数，调整箭头的大小
        quiver(r(1:ip:end,1:ip:end),z(1:ip:end,1:ip:end),...
            u(1:ip:end,1:ip:end)*scale_factor,v(1:ip:end,1:ip:end)*scale_factor,'k','AutoScale', 'off',...
            'LineWidth',2);
 
        % 创建 contour
        contourf(r,z, bdr_ph, 'LineWidth',3,'EdgeColor',[0 0 0],...
        'FaceColor',[0.8 0.8 0.8],...
        'LevelList',0.5);

        axis([0 3.5 0 1.5]);
        title(['t=' sprintf('%.2f', t)]);
        
             % 捕捉当前帧并写入视频
            frame = getframe(gcf);
            writeVideo(outputVideo, frame);
        
    end
    
    low_bdr = ly/4; 
    %[theta_fit] = fit_circle(r,z,f2, low_bdr, low_bdr, true, 10, 'b')
    

    u_c = sum( sum( ff .* u )) /  sum( sum( ff ));
    
    u_c_lst = [u_c_lst; u_c];
 
    max(max(abs(u)))
    t_lst = [t_lst; t];
    
    i,t
    pause(0.5);

    %clear r z u v p f t;
    
end

close(outputVideo);

end

function [] = no_use()
    figure
    plot(t_lst, y_lst);
    hold on;
    
%     % 1% position inside
%     percent = 0.01;
%     i_index = find(y_lst > y_lst(end)*(1.0+percent), 1, 'last') + 1;
%     
%     plot([0 80], [y_lst(end)*(1.0+percent) y_lst(end)*(1.0+percent)]);
%     
%     num_of_round = 0;
%     for i = 2:i_index
%         if x_lst(i) < x_lst(i-1)
%             num_of_round = num_of_round +1;
%         end
%     end
%     
%     % every roudtrip is 2.0 and last point added x_lst(i_index) - 1.0, the
%     % start point is 0.5
%     l_e = num_of_round*1.0 + x_lst(i_index) - 1.0 + 0.5;
    
    
    figure
    plot(t_lst, x_lst);
    
    figure;
    mesh(r,z,u);
    title([ 't=' num2str(t)])
    
    figure; hold on;
    contour(r,z,u, [0.1:0.1:0.7 0.7:0.02:1.0]);
    C = contour(r,z,f,[0.5 0.5], 'm');
    axis equal;

    figure; hold on;
    contour(r,z,u, 0.:0.02:1.0);
    C = contour(r,z,f,[0.5 0.5], 'm');
    axis equal;
    
    
    
%     h_e = [0.429, 0.3998, 0.393, 0.386] - 0.35/2.0;
%     l_e = [6.4339, 2.7086/4.5058, 2.8369,  1.3472] 

    % time stop
    i_index = find(t_lst < 3, 1, 'last') + 1;
       
    num_of_round = 0;
    for i = 2:i_index
        if x_lst(i) < x_lst(i-1)
            num_of_round = num_of_round +1;
        end
    end
    
    % every roudtrip is 2.0 and last point added x_lst(i_index) - 1.0, the
    % start point is 0.5
    l_e = num_of_round*3.5/2.0 + x_lst(i_index) - (3.5/2.0-0.5) - (3.5/2.0-0.5)

    
    dr = r(1,2)-r(1,1); dz = z(2,1)-z(1,1);

pp = (p(1:end-2,2:end-1) + p(3:end,2:end-1) - 2.0*p(2:end-1,2:end-1))/dr^2.0 ...
+ (p(2:end-1, 1:end-2) + p(2:end-1,3:end) - 2.0*p(2:end-1,2:end-1))/dz^2.0;


div = (u(2:end,2:end) - u(1:end-1,2:end))/dr ...
    + (v(2:end,2:end)- v(2:end,1:end-1))/dz;

end

