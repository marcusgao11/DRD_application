

function [x_lst, y_lst, Area_lst, t_lst, u_c_lst, r,z,u,v]=da2_new(N_start,dN,N_end,M,N,data_path, indicator, quiver_params)
%close all;
% M = 257
% N = 513

is_text = true;
is_plot = quiver_params.is_plot;
is_pause = quiver_params.is_pause;

x_lst = []; y_lst = []; Area_lst = []; t_lst=[];
u_c_lst = [];


is_write_vedio = quiver_params.is_write_vedio;
if (is_write_vedio)
    outputVideo = VideoWriter(quiver_params.movie_nm);
    outputVideo.FrameRate = 50;  % 设置帧率
    % 打开文件进行写入
    open(outputVideo);
end

screenSize = get(0, 'ScreenSize');
% screenSize的格式为[left, bottom, width, height]

ch = ['0' '1' '2' '3' '4' '5' '6' '7' '8' '9'];
pic_num  = 1;

for i=N_start:dN:N_end
       %k1=floor(i/10);
       %k2=i-k1*10;
    
    k4 = floor(i/1000);
    k3 = floor( (i-k4*1000) /100);
    k2 = floor( (i- k4*1000 - k3*100) /10 );
    k1 = mod(i,10);
    
    fname = ['data' ch(k4+1) ch(k3+1) ch(k2+1) ch(k1+1)];
    if ~isempty(indicator)
        fname = [fname indicator];
    end
    disp([data_path fname]);
    
    [r,z,u,v,p,f, f_2, pf1, pf2,t]=da1([data_path fname],M+1,N+1);
    
%     pcolor(r,z,f);
%     axis equal
%     shading interp;
    
if (is_plot) 
    figure(1);
    %figHandle = findobj('Type', 'figure', 'Number', 1);
    set(gcf,  'OuterPosition', [927.4 42.6 1100 1100]);
    %set(figHandle, 'Position', [screenSize(3)/4*2.5, 100, 800, 600]);

%     %contour(r,z,f,[1.0E-10 1.0E-10], 'm');
%     C = contour(r,z,f,[0.5 0.5], 'm');
%     axis equal;
    
    %contour(r,z,f,[1.0E-10 1.0E-10], 'm');
    hold off;

    %%C = contour(r,z,f_2,[0.5 0.5], 'k'); hold on;
    

    %contourf(r,z, f_2*0.8/0.5, [0.5 0.5]*0.8/0.5);  hold on;  % swimmers
    contourf(r,z, f_2*0.8/0.5, 'LineWidth',2,'FaceColor',[0.8 0.8 0.8],'LevelList',0.8); hold on;
    
    contourf(r,z,(f-f_2),'LineWidth',2,'EdgeColor',[0 0 0],'FaceColor',[1 0 0],'LevelList',0.5);
    
    %C = contourf(r,z,(f-f_2),[0.5 0.5], 'm'); %colormap([0 0 1; 0 0 1;]);  hold on; % fibers
    colormap jet; %caxis([0 1]); 
    %colorbar;      
        
    %hold on; %,contour(x,y,r');
    %axis([0 Lx 0 Ly]); axis equal;
    ip = 4;
    
    scale_factor = quiver_params.scale/max(max(sqrt(u.^2+v.^2)));     % 缩放系数，调整箭头的大小    
    on_off = quiver_params.on_off;
    
    % 创建 ylabel
    ylabel('$z/H$','Interpreter','latex');

    % 创建 xlabel
    xlabel('$x/H$','Interpreter','latex');
    
    if (on_off) 
        % 创建 quiver
        quiver(r(1:ip:end,1:ip:end),z(1:ip:end,1:ip:end), u(1:ip:end,1:ip:end),v(1:ip:end,1:ip:end),...
            'MaxHeadSize',1,'LineWidth',0.8,'Color',[0 0 1], 'AutoScaleFactor',2,'AlignVertexCenters','on');
        %quiver(r(1:ip:end,1:ip:end),z(1:ip:end,1:ip:end),...
        %    u(1:ip:end,1:ip:end)*scale_factor,v(1:ip:end,1:ip:end)*scale_factor,'b','AutoScale', 'off');
        hold on; %axis([0 1 0 1]);    
        %contourf(r,z,f,[0.5 0.5], 'm'); %colormap([0 0 1; 0 0 1;]);  hold on; cirles
    end
    
    %axis([2,4,2,4]);
end

    if (is_text)
        % plot arrow
        swm = [data_path fname '_swm.txt'];
        swm = load(swm);

        if (~isempty(swm))
            ii_swm = 1:size(swm,1);
            x = swm(ii_swm,2); y = swm(ii_swm,3);
            x_u = cos(swm(ii_swm,4)); y_v = sin(swm(ii_swm,4));
            if (is_plot)
                %quiver(x,y,x_u,y_v, 0.1, 'k');
                %plot(x,y,'ro');
                quiver(x,y,x_u,y_v,...
                    'MaxHeadSize',1,'LineWidth',2, 'Color',[0 0 0], 'AutoScaleFactor',0.2);
            end
        end

        fib = [data_path fname '_fib.txt'];
        fib = load(fib);
        if (~isempty(fib))
            x = fib(:,2); y = fib(:,3);
            
            if (is_plot)
                plot(x,y,'*');
            end
            
            x_lst = [x_lst;  x(floor(length(x)/2)+1)];
            y_lst = [y_lst;  y(floor(length(x)/2)+1)];
        end
        
    end
    
%     %axis equal; 
%     axis([0 1 0 1]);    
%     xticks([0 0.2 0.4 0.6 0.8 1]); xticklabels({'0','0.2','0.4','0.6','0.8','1'})
%     set(gca, 'FontSize',45);

    % 取消以下行的注释以保留坐标区的 X 范围
    xlim(gca,[0 1]);
    % 取消以下行的注释以保留坐标区的 Y 范围
    ylim(gca,[0 1]);
    box(gca,'on');
    hold(gca,'off');
    % 设置其余坐标区属性
    set(gca,'BoxStyle','full','DataAspectRatio',[1 1 1],'FontName','Serif',...
        'FontSize',45,'Layer','top','LineWidth',3,'PlotBoxAspectRatio',[1 1 2],...
        'TickLength',[0.005 0.025],'XTick',[0 0.2 0.4 0.6 0.8 1],'XTickLabel',...
        {'0','0.2','0.4','0.6','0.8','1'},'YTick',[0 0.2 0.4 0.6 0.8 1],'YTickLabel',{'0','0.2','0.4','0.6','0.8','1'});



if (is_plot)
    title(['t/' 't' '\fontsize{25}{0}' '\fontsize{45} =' sprintf('%.2f', t)]);
    %title(['t=' num2str(t) ' max(abs(u)), max(abs(v))='  num2str(max(max((u)))) ', ' num2str(max(max((v)))) ])
    %quiver(r,z,u,v);
    
    hold off;
end
    
    if (is_write_vedio)
        F=getframe(gcf);
        I=frame2im(F);
        [I,map]=rgb2ind(I,8);
        if(pic_num == 1)
            imwrite(I,map,'test.gif','gif','Loopcount',inf,'DelayTime',0.02);
        else
            imwrite(I,map,'test.gif','gif','WriteMode','append','DelayTime',0.02);
        end
        pic_num = pic_num + 1;
        
%             % 捕捉当前帧并写入视频
%             frame = getframe(gcf);
%             writeVideo(outputVideo, frame);
    end
            
% if (is_plot)
%     
%     figure(2);
%     mesh(r,z,u); title('u');
% 
%     figure(3);
%     mesh(r,z,v); title('v');
%     
%     figure(4);
%     mesh(r,z,pf1);
%     
%     figure(5);
%     mesh(r,z,pf2);
%     
% end
%     
    %u_c = sum( sum( f .* u )) /  sum( sum( f ));
    
    % 对于swimmer找到前面的active part， 先找到中间位置，再分割
    split_pos = sum( sum( f .* r )) /  sum( sum( f ));
    
    ff = f;
    ff(r<split_pos) = NaN;
    u_c = nansum( nansum( ff .* u )) /  nansum( nansum( ff ));
    
    %x_lst = [x_lst; nansum( nansum( ff .* r )) /  nansum( nansum( ff ))];
    
    u_c_lst = [u_c_lst; u_c]; t_lst = [t_lst; t];
    
    %max(max(abs(u))) + max(max(abs(v)))
    
    sum(sum(u))
      
    if (is_pause)
       pause(0.5);
    end
    
    i,t
    %clear r z u v p f t;
    
end

if (is_write_vedio)
    close(outputVideo);
end

end


