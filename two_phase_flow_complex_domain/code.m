% 11-4-2010 

close all;

M=129;
N=129;   

M=257;
N=257; 
% % 
% M=385;
% N=385; 
% 
M=513;
N=513;
% % % % 
% M=1025;
% N=1025;
%   
% M=513
% N=257

fname='data010';


cd ~
cd CouettleFlow/implicitscheme/pressure_cor_2nd/periodic_bc/

cd .. 
cd ..

%cd 2ndorderGNBC/periodic_bc%/small_dt_data/
cd pressure_cor_2nd/periodic_bc/
%cd pressure_cor_2nd/periodic_bc/


[r,z,u,v,p,f,t]=da1(fname,M,N);

t


cd .. 
cd ..

%cd ..


%fname='data010';

cd 2ndorderGNBC/periodic_bc%/small_dt_data/
[r,z,u2,v2,p2,f2,t2]=da1(fname,M,N);

't-t2'
t-t2
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)*3/5 scrsz(3)/2 scrsz(4)*2/5])

mesh(u2-u)
% mesh(u-u2)
title('u2-u');

figure('Position',[scrsz(3)/2 scrsz(4)*3/5 scrsz(3)/2 scrsz(4)*2/5])

mesh(v2-v)
% mesh(v-v2)
title('v2-v');


figure('Position',[1 1 scrsz(3)/2 scrsz(4)*2/5])

mesh(p2-p)
% mesh(p-p2)
title('p2-p');


figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)*2/5])
mesh(f2-f)
% mesh(f-f2)
title('\phi2-\phi');
pause


% dr=r(1,2)-r(1,1)
% dz=z(2,1)-z(1,1)
% 
% 
% M = M-1;
% N = N-1;
% 
% Dxp=zeros(N+1,M+1);
% Dyp=zeros(N+1,M+1);
% Dxp2=zeros(N+1,M+1);
% Dyp2=zeros(N+1,M+1);
% 
% for j=2:N
%     for i=2:M
%         Dxp(j,i) = ( p(j,i+1)-p(j,i-1) )/2/dr;
%         Dyp(j,i) = ( p(j+1,i)-p(j-1,i) )/2/dz;
%         Dxp2(j,i) = ( p2(j,i+1)-p2(j,i-1) )/2/dr;
%         Dyp2(j,i) = ( p2(j+1,i)-p2(j-1,i) )/2/dz;
%     end
% end
% 
% mesh(r,z,Dxp-Dxp2);
% pause 
% mesh(r,z,Dyp-Dyp2);





% slen=0.4;
% shig=0.4; 
% 
% 
% ue= slen*sin(t)*(sin(pi*r/slen)).^2.*sin(2*pi*z/shig);
%  
% ve= -shig*sin(t)*sin(2*pi*r/slen).*(sin(pi*z/shig)).^2;
% 
% pe=  (slen*shig*sin(t)*cos(2*pi*z/shig).*sin(2*pi*r/slen))/2/pi;
%  
% fe = cos((pi*z)/shig)*sin(t).*sin((2*pi*r)/slen);
% 
% 
% 
% scrsz = get(0,'ScreenSize');
% figure('Position',[1 scrsz(4)*3/5 scrsz(3)/2 scrsz(4)*2/5])
% 
% mesh(u-ue)
% % mesh(u-u2)
% title('u-ue');
% 
% 
% 
% figure('Position',[scrsz(3)/2 scrsz(4)*3/5 scrsz(3)/2 scrsz(4)*2/5])
% 
% mesh(v-ve)
% % mesh(v-v2)
% title('v-ve');
% 
% 
% 
% figure('Position',[1 1 scrsz(3)/2 scrsz(4)*2/5])
% 
% mesh(p-pe)
% % mesh(p-p2)
% title('p-pe');
% 
% 
% 
% figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)*2/5])
% mesh(f-fe)
% % mesh(f-f2)
% title('\phi-\phi_e');
%pause

close all;
