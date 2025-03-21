%%plot the result as time evolves
function []=da2(M,N,dN)
close all;
% M = 257
% N = 513

ch = ['0' '1' '2' '3' '4' '5' '6' '7' '8' '9'];
for i=1:dN:999
       %k1=floor(i/10);
       %k2=i-k1*10;
    
    k3 = floor(i/100);
    k2 = floor((i-k3*100)/10);
    k1 = mod(i,10);
    [r,z,u,v,p,f,t]=da1(['data' ch(k3+1) ch(k2+1) ch(k1+1)],M,N);
    pcolor(r,z,f);
    axis equal
    shading interp;
    pause(0.5);
    i,t
    clear r z u v p f t;
end

end
