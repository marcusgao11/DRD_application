% m=256; n =128;
% [r,z,u,v,p,f,t]=da1('E:\numerical_simulation\tmp\data001',m+1,n+1);
% [r,z,u,v,p,f,t]=da1(['E:\numerical simulation\share\data070'],m+1,n+1);

function [r,z,u,v,p,f,c,t]=da1(file,m,n)
%close all

ff_info = dir(file);
ff_info.date

fid=fopen(file,'r');
fseek(fid,4,'cof');
r=fread(fid,[n,m],'double');
fseek(fid,8,'cof');
z=fread(fid,[n,m],'double');
fseek(fid,8,'cof');
u=fread(fid,[n,m],'double');
fseek(fid,8,'cof');
v=fread(fid,[n,m],'double');
fseek(fid,8,'cof');
p=fread(fid,[n,m],'double');
fseek(fid,8,'cof');
f=fread(fid,[n,m],'double');
fseek(fid,8,'cof');

c=fread(fid,[n,m],'double');
fseek(fid,8,'cof');

t=fread(fid,1,'double');

r = r'; z = z'; u = u';
v = v'; p = p'; f = f'; c = c';

fclose(fid);



% h=z(2,1)-z(1,1);
% x=r(2,:);
% p1=p(2,:);
% u1=u(2,:);
% f1=f(2,:);
% for i=2:m-2
% px(i)=(p1(i+1)-p1(i-1))/(x(i+1)-x(i-1));
% fx(i)=(f1(i+1)-f1(i-1))/(x(i+1)-x(i-1));
% ux(i)=(u1(i+1)-u1(i-1))/(x(i+1)-x(i-1));
%         h1=x(i+1)-x(i);
%         h2=x(i)-x(i-1);
%         h3=x(i+2)-x(i);
%         a= 2.*(h2-h3)/(h1*(h1+h2)*(h1-h3));
%         b= 2.*(h1+h3)/(h2*(h1+h2)*(h2+h3));
%         c= 2.*(h1-h2)/(h3*(h3+h2)*(h1-h3));
%       uux(i)=u1(i).*u1(i);
%       uxx(i)=a*u1(i+1)+b*u1(i-1)+c*u1(i+2)-(a+b+c)*u1(i);
%       uyy(i)=(u(3,i)-2*u(2,i)+u(1,i))/(h*h);
%       fy(i)=(f(3,i)-f(1,i))/(2.*h);
%       uy(i)=(u(3,i)-u(1,i))/(2.*h);
% end
%subplot(1,2,1)
%plot(x(2:m-2),abs(px(2:m-2)))
%hold
%plot(x(2:m-2),abs(uxx(2:m-2)),'r')
%plot(x(2:m-2),abs(uyy(2:m-2)),'g')
%plot(x(2:m-2),abs(0.03*uux(2:m-2)),'k')
%figure
%subplot(1,2,2)
%plot(x(2:m-2),abs(uy(2:m-2)),'g')
%hold
%plot(x(2:m-2),abs(12.*fx(2:m-2).*fy(2:m-2)),'k')


