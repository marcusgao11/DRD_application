

%%plot the result as time evolves

%m=256; n =128;
%da2(1,1,100,m,n, 'E:\numerical_simulation\tmp\', '')

% da2(1,1,100,m,n, 'E:\numerical simulation\sav_q\','');

% da2(1,1,100,m,n, 'E:\numerical simulation\share\','');

% m=256; n =512;
% da2(1,1,100,m,n, 'E:\numerical simulation\share\','');

function [t_lst, h_lst]=da2_2(N_start,dN,N_end,M,N,data_path, indicator, is_plot)

    if nargin < 8 || isempty(is_plot)
        is_plot = false; % 将defaultValue替换为你想要的默认值
    end
    
% M = 257
% N = 513

h_lst = []; t_lst = [];

outputVideo = VideoWriter('MICP_20240622.avi');
outputVideo.FrameRate = 0.8;  % 设置帧率
% 打开文件进行写入
open(outputVideo);

i_count = 0;   colors = {'b', 'g', 'r', 'c', 'm', 'k'};
ch = ['0' '1' '2' '3' '4' '5' '6' '7' '8' '9'];

i = 0;
k3 = floor(i/100);
k2 = floor((i-k3*100)/10);
k1 = mod(i,10);
fname = ['data' ch(k3+1) ch(k2+1) ch(k1+1)];
[r0,z0,u0,v0,p0,f0,c0,t0] = da1([data_path fname],M+1,N+1);

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
    
    try
        [r,z,u,v,p,f,c,t] = da1([data_path fname],M+1,N+1);
    catch ME
        fprintf('捕获到错误: %s\n', ME.message);
        break;
    end
         
%     pcolor(r,z,f);
%     axis equal
%     shading interp;
    
    if is_plot
        i_count = i_count + 1;
        
        figure(1);
        mesh(r,z,c);
        title(['c, t=' num2str(t)]);

        figure(3)
        contour(r,z,f);

        figure(2); hold on;
        %contour(r,z,f,[1.0E-10 1.0E-10], 'm');
        if i == 0
            [C h] = contour(r,z,f,[1.0E-10 1.0E-10], 'k');
        else
            [C h] = contour(r,z,f,[1.0E-10 1.0E-10], colors{mod(i_count, length(colors))+1});
        end

        axis equal;
        
        figure(4);  
        set(gcf,  'OuterPosition',[602.333333333333 266 1515.66666666667 892]);

        hold off;
        
        %mesh(r,z,u);
        contour(r,z,f,'LineWidth',3,'LineColor',[0 0 0],...
            'LevelList',1e-06);
        hold on;
               
        contour(r,z,f0,'LineWidth',3,'LineStyle',':',...
            'LineColor',[0 0 0],...
            'LevelList',1e-06);
        
        % 创建 ylabel
        ylabel('$z/H$','Interpreter','latex');

        % 创建 xlabel
        xlabel('$x/H$','Interpreter','latex');
        set(gca, 'FontSize',45);

        cc = c;
        cc(f > 0.0) = NaN;
        pcolor(r,z,cc); axis equal;
        shading interp; colorbar; hold on;
        caxis([0.8 0.9]);
        axis equal; colormap([jet(128); 1 0 0]);
        %axis off
 
        inter = 10;
        scale_factor = 0.02;
        quiver(r(1:inter:end,1:inter:end),z(1:inter:end,1:inter:end), ...
          u(1:inter:end,1:inter:end)*scale_factor,v(1:inter:end,1:inter:end)*scale_factor, ...
          'Color', 'k', 'AutoScale', 'off', 'LineWidth',2);
        %'AutoScale', 'on', 'AutoScaleFactor', 1.5, 'LineWidth',2);
        %title(['t=' num2str(t)]);
        title(['t=' sprintf('%.2f', t)]);
      
        axis([0 2 0 1]);
  
          % 捕捉当前帧并写入视频
            frame = getframe(gcf);
            writeVideo(outputVideo, frame);

    %      cc = c(1,:); c0 = 0.8; dz = z(1,2) - z(1,1);
    %      dc_diff = ((cc / c0).^2 - 1) * (1 - c);
    %      dc = zeros(size(cc));
    %      dc(2:end-1) = - (cc(3:end)-cc(1:end-2))/2.0/dz;
    %      figure(5); hold on;
    %      plot(z(1,:), dc_diff);
    %      plot(z(1,:), dc);     

        %,contour(x,y,r');
        figure(5); hold off;
        contour(r,z,f,[1.0E-10 1.0E-10], 'm'); hold on;
        %axis([0 Lx 0 Ly]); axis equal;
        %quiver(r,z,u,v,'r'); 
        inter = 5;
        quiver(r(1:inter:end,1:inter:end),z(1:inter:end,1:inter:end),u(1:inter:end,1:inter:end),v(1:inter:end,1:inter:end));
        axis equal;
        title(['t=' num2str(t)]);
        
        
        figure(6);
        cc = c;
        cc(f > 0.0) = NaN;
        pcolor(r,z,cc); axis equal;
        shading interp; colorbar;
        title('c');

        %hold off;
    else
        fig = figure('visible', 'off');
        [C h] = contour(r,z,f,[1.0E-10 1.0E-10], 'm');
    end
    
    try
        h_lst = [h_lst; C(2,2)]; 
    catch ME
        fprintf('捕获到错误: %s\n', ME.message);
        break;
    end
      
    if is_plot
        pause(0.5);
    end
    
    t_lst = [t_lst; t];
    
    i,t
    clear r z u v p f t;
end

close(outputVideo);

end
