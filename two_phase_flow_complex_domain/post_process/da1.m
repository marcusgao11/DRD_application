function [r,z,u,v,p,f,f2, t]=da1(file,m,n)
    
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
    f2=fread(fid,[n,m],'double');

    fseek(fid,8,'cof');
    t=fread(fid,1,'double');
    fclose(fid);

end