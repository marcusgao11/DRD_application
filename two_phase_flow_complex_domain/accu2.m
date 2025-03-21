fid1=fopen('test1.txt','w');
fid2=fopen('test2.txt','w');

format short e;

close all;

M = 128;
N = 64;

fname = 'smal001';

fprintf(fid1,'%s\n','(nr,nz) & \phi & rate& u & rate& v & rate& p & rate\\');
fprintf(fid2,'%s\n','(nr,nz) & \phi & rate& u & rate& v & rate& p & rate\\');


fprintf(fid1,'%i%s',M,' & '); 
fprintf(fid2,'%i%s',M,' & ');


%[e1,e2,e3,e4,e11,e22,e33,e44]=accuracy(128,'data002','data004');
%[e1,e2,e3,e4,e11,e22,e33,e44]=accuracy(128,'d_np002','d_np004');

cd ac1/
[r,z,u,v,p,f,t]=da1(fname,M+1,N+1);

cd ..
cd ac2/

[r1,z1,u1,v1,p1,f1,t1]=da1(fname,2*M+1,2*N+1);
cd ..


[e1,e2,e3,e4,e11,e22,e33,e44]...
    =accuracy(M,N,r1,z1,u1,v1,p1,f1,t1,r,z,u,v,p,f,t)

fprintf(fid1,'%7.2e',e1);
fprintf(fid1,'%s \t\t %s',' & ',' & ');

fprintf(fid1,'%7.2e',e2);
fprintf(fid1,'%s \t\t %s',' & ',' & ');

fprintf(fid1,'%7.2e',e3);
fprintf(fid1,'%s \t\t %s',' & ',' & ');

fprintf(fid1,'%7.2e',e4);
fprintf(fid1,'%s \t\t \n',' &\\ ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
fprintf(fid2,'%7.2e',e11);
fprintf(fid2,'%s \t\t %s',' & ',' & ');
        
fprintf(fid2,'%7.2e',e22);
fprintf(fid2,'%s \t\t %s',' & ',' & ');

fprintf(fid2,'%7.2e',e33);
fprintf(fid2,'%s \t\t %s',' & ',' & ');

fprintf(fid2,'%7.2e',e44);
fprintf(fid2,'%s \t\t \n',' &\\');



 fprintf(fid1,'%s','(256,128)-(512,256)   & ');
 fprintf(fid2,'%s','(256,128)-(512,256)   & ');
 
 
cd ac4/
[r,z,u,v,p,f,t]=da1(fname,M*4+1,N*4+1);

cd ..
[ee1,ee2,ee3,ee4,ee11,ee22,ee33,ee44]...
    =accuracy(2*M,2*N,r,z,u,v,p,f,t,r1,z1,u1,v1,p1,f1,t1)

%[ee1,ee2,ee3,ee4,ee11,ee22,ee33,ee44]=accuracy(256,'data004','data008');
%[ee1,ee2,ee3,ee4,ee11,ee22,ee33,ee44]=accuracy(256,'d_np004','d_np008');

fprintf(fid1,'%7.2e',ee1);
fprintf(fid1,'%s %4.2f %s',' & ',log2(e1/ee1),' & ');

fprintf(fid1,'%7.2e',ee2);
fprintf(fid1,'%s %4.2f %s',' & ',log2(e2/ee2),' & ');

fprintf(fid1,'%7.2e',ee3);
fprintf(fid1,'%s %4.2f %s',' & ',log2(e3/ee3),' & ');

fprintf(fid1,'%7.2e',ee4);
fprintf(fid1,'%s %4.2f %s \n',' & ',log2(e4/ee4),'\\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid2,'%7.2e',ee11);
fprintf(fid2,'%s %4.2f %s',' & ',log2(e11/ee11),' & ');

fprintf(fid2,'%7.2e',ee22);
fprintf(fid2,'%s %4.2f %s',' & ',log2(e22/ee22),' & ');

fprintf(fid2,'%7.2e',ee33);
fprintf(fid2,'%s %4.2f %s',' & ',log2(e33/ee33),' & ');

fprintf(fid2,'%7.2e',ee44);
fprintf(fid2,'%s %4.2f %s \n',' & ',log2(e44/ee44),'\\');



fprintf(fid1,'%s','(512,256)-(1024,512)  & ');
fprintf(fid2,'%s','(512,256)-(1024,512)  & ');
% 
% [e1,e2,e3,e4,e11,e22,e33,e44]=accuracy(512,'data008','data016');
% %[e1,e2,e3,e4,e11,e22,e33,e44]=accuracy(512,'d_np008','d_np016');

cd ac8/
[r1,z1,u1,v1,p1,f1,t1]=da1(fname,M*8+1+2,N*8+1);

cd ..
[e1,e2,e3,e4,e11,e22,e33,e44]...
    =accuracy(4*M+1,4*N,r1,z1,u1,v1,p1,f1,t1,r,z,u,v,p,f,t)

fprintf(fid1,'%7.2e',e1);
fprintf(fid1,'%s %4.2f %s',' & ',log2(ee1/e1),' & ');

fprintf(fid1,'%7.2e',e2);
fprintf(fid1,'%s %4.2f %s',' & ',log2(ee2/e2),' & ');

fprintf(fid1,'%7.2e',e3);
fprintf(fid1,'%s %4.2f %s',' & ',log2(ee3/e3),' & ');

fprintf(fid1,'%7.2e',e4);
fprintf(fid1,'%s %4.2f \n',' & ',log2(ee4/e4),'\\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid2,'%7.2e',e11);
fprintf(fid2,'%s %4.2f %s',' & ',log2(ee11/e11),' & ');

fprintf(fid2,'%7.2e',e22);
fprintf(fid2,'%s %4.2f %s',' & ',log2(ee22/e22),' & ');

fprintf(fid2,'%7.2e',e33);
fprintf(fid2,'%s %4.2f %s',' & ',log2(ee33/e33),' & ');

fprintf(fid2,'%7.2e',e44);
fprintf(fid2,'%s %4.2f \n',' &',log2(ee44/e44),'\\');
% 
% 
% 
% fprintf(fid1,'%s','(1024,512)-(2048,1024)& ');
% fprintf(fid2,'%s','(1024,512)-(2048,1024)& ');
% 
% [ee1,ee2,ee3,ee4,ee11,ee22,ee33,ee44]=accuracy(1024,'data016','data032');
% %[ee1,ee2,ee3,ee4,ee11,ee22,ee33,ee44]=accuracy(1024,'d_np016','d_np032');
% 
% fprintf(fid1,'%7.2e',ee1);
% fprintf(fid1,'%s %4.2f %s',' & ',log2(e1/ee1),' & ');
% 
% fprintf(fid1,'%7.2e',ee2);
% fprintf(fid1,'%s %4.2f %s',' & ',log2(e2/ee2),' & ');
% 
% fprintf(fid1,'%7.2e',ee3);
% fprintf(fid1,'%s %4.2f %s',' & ',log2(e3/ee3),' & ');
% 
% fprintf(fid1,'%7.2e',ee4);
% fprintf(fid1,'%s %4.2f %s \n',' & ',log2(e4/ee4),'\\');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf(fid2,'%7.2e',ee11);
% fprintf(fid2,'%s %4.2f %s',' & ',log2(e11/ee11),' & ');
% 
% fprintf(fid2,'%7.2e',ee22);
% fprintf(fid2,'%s %4.2f %s',' & ',log2(e22/ee22),' & ');
% 
% fprintf(fid2,'%7.2e',ee33);
% fprintf(fid2,'%s %4.2f %s',' & ',log2(e33/ee33),' & ');
% 
% fprintf(fid2,'%7.2e',ee44);
% fprintf(fid2,'%s %4.2f %s \n',' & ',log2(e44/ee44),'\\');




fclose(fid1);
fclose(fid2);

open('test1.txt');
open('test2.txt');