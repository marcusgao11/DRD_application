function [e1,e2,e3,e4,e11,e22,e33,e44]...
    =accuracy(M,N,r1,z1,u1,v1,p1,f1,t1,r,z,u,v,p,f,t)
%[r,z,u1,v1,p1,f1,t1,u,v,p,f,t]=

%M=512;
%N=M/2;
%M=1025;
%N=513;
kk = 2;
close all;


%ch = ['1' '2' '4' '8' '16' '32'];

%    fname='data001';

    dr = 0.5/(M);
    dz = 1.0/(N);

%    format short e;

%    cd at16/

%    [r1,z1,u1,v1,p1,f1,t1]=da1(fname2,2*M+1,2*N+1);
    
    f1=f1(1:kk:2*N+1,1:kk:2*M+1);
    u1=u1(1:kk:2*N+1,1:kk:2*M+1);
    v1=v1(1:kk:2*N+1,1:kk:2*M+1);
    p1=p1(1:kk:2*N+1,1:kk:2*M+1);

%    cd ..

%    cd at8/

%    [r,z,u,v,p,f,t]=da1(fname1,M+1,N+1);

%    cd ..

     t
     t1



     ff = (f-f1); %/max(max(abs(f)));
     disp('\phi max');
     e1 = max(max(abs(ff)))

     
     disp('L2')
     e11 = L2_error(M,N,dr,dz,f-f1,f)

     ff = ff(1:kk:N,1:kk:M);
    % mesh(ff)
    % title('\phi')




     uu=(u-u1); %/max(max(abs(u)));
     disp('u max');
     e2 = max(max(abs(uu)))

     disp('L2')
     e22 = L2_error(M,N,dr,dz,u-u1,u)

     uu = uu(1:kk:N,1:kk:M);
     %mesh(uu)
     %title('u')
     %pause



     vv=(v-v1); %/max(max(abs(v)));
     disp('v max');
     e3 = max(max(abs(vv)))

     disp('L2')
     e33 = L2_error(M,N,dr,dz,v-v1,v)

     vv = vv(1:kk:N,1:kk:M);
     %(vv)
     %title('v')
%     pause

     Cx = ones(M,1);
     Cz = ones(N,1);

     Cx(1) = 0.5;
     Cx(M) = 0.5;
     Cz(1) = 0.5;
     Cz(N) = 0.5;

          
% figure(1)
% mesh(p)
% figure(2)
% mesh(p1)
% pause

     pp = p-p1;
     sum1 = 0;
     sum2 = 0;
     sum3 = 0;

     for j=1:N
         for i=1:M
             sum1 = sum1 + Cx(i)*Cz(j);
             sum2 = sum2 + Cx(i)*Cz(j)*pp(j,i);
             sum3 = sum3 + Cx(i)*Cz(j)*p(j,i);
         end
     end

     sum2 = sum2/sum1;
     sum3 = sum3/sum1;


     pp = (p-p1-sum2); %/max(max(abs(p-sum3)));

     disp('p max');
     e4 = max(max(abs(pp)))
     
     

     disp('L2')
     e44 = L2_error(M,N,dr,dz,(p-p1-sum2),p-sum3)

     %pp = pp(1:kk:N,1:kk:M);
     %mesh(pp);
     %title('p')
     %pause
     %close
    figure(1)
    %mesh(p);
    plot(p(1,:));
    figure(2)
    %mesh(p1);
    plot(p1(1,:));
    pause     
    mesh(pp);
%      contour(f,1.0E-8);
%      hold on
%      contour(f1,1.0E-8);
t
t1
     pause




    function L2=L2_error(M,N,dr,dz,u1,u2)

     Cx = ones(M,1);
     Cz = ones(N,1);

     Cx(1) = 0.5;
     Cx(M) = 0.5;
     Cz(1) = 0.5;
     Cz(N) = 0.5;

     sum2 = 0;
     sum3 = 0;

     M
     N
     size(u2)
     size(u2)
     for j=1:N
         for i=1:M
             sum2 = sum2 + Cx(i)*Cz(j)*u1(j,i)^2;
             sum3 = sum3 + Cx(i)*Cz(j)*u2(j,i)^2;
         end
     end

     sum2 = sqrt(sum2*dr*dz);
     sum3 = sqrt(sum3*dr*dz);

     L2 = sum2;%/sum3;
     

% close all;
% clear all;
% M=64;
% N=32;
% 
% cd at2/
% [r,z,u,v,p,f,t]=da1('data001',2*M+1,2*N+1);
% 
% cd ..
% cd at1/ 
% [r,z,u1,v1,p1,f1,t1]=da1('data001',M+1,N+1);
% 
% cd ..
% 
% ff=f(1:2:2*N+1,1:2:2*M+1);
% uu=u(1:2:2*N+1,1:2:2*M+1);
% vv=v(1:2:2*N+1,1:2:2*M+1);
% pp=p(1:2:2*N+1,1:2:2*M+1);
