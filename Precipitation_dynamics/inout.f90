    subroutine output(M,N,x,y,u,v,p,f,c,t,fname)
    implicit none
    integer M,N,i,k
    real x(0:M,0:N),y(0:M,0:N),u(-1:M+1,-1:N),v(-1:M,0:N)
    ! careful: here p is not the same format as u,v,p
    real p(-1:M,-1:N),f(-1:M,-1:N),c(-1:M,-1:N)
    real t
    character fname*7
    open(11,file=fname, status='unknown',form='unformatted')
    write(11)((x(i,k),k=0,N),i=0,M)
    write(11)((y(i,k),k=0,N),i=0,M)
    write(11)((u(i,k),k=0,N),i=0,M)
    write(11)((v(i,k),k=0,N),i=0,M)
    write(11)((p(i,k),k=0,N),i=0,M)
    write(11)((f(i,k),k=0,N),i=0,M)
    write(11)((c(i,k),k=0,N),i=0,M)
    write(11)t
    close(11)
    return
    end

    subroutine output_regular(M,N,x,y,u,v,p,f,t,fname)
    implicit none
    integer M,N,i,k
    real x(0:M,0:N),y(0:M,0:N),u(-1:M+1,-1:N+1),v(-1:M+1,-1:N+1),f(-1:M+1,-1:N+1),p(-1:M+1,-1:N+1)
    ! careful: here p is not the same format as u,v,p
    real t
    character fname*7
    open(11,file=fname, status='unknown',form='unformatted')
    write(11)((x(i,k),k=0,N),i=0,M)
    write(11)((y(i,k),k=0,N),i=0,M)
    write(11)((u(i,k),k=0,N),i=0,M)
    write(11)((v(i,k),k=0,N),i=0,M)
    write(11)((p(i,k),k=0,N),i=0,M)
    write(11)((f(i,k),k=0,N),i=0,M)
    write(11)t
    close(11)
    return
    end

    subroutine input(M,N,x,y,u,v,p,f,c,t,fname)
    implicit none
    integer M,N,i,k
    real x(0:M,0:N),y(0:M,0:N),u(-1:M+1,-1:N),v(-1:M,0:N)
    ! careful: here p is not the same format as u,v,p
    real p(-1:M,-1:N),f(-1:M,-1:N),c(-1:M,-1:N)
    real t
    character fname*7

    open(31,file=fname, status='old',form='unformatted')
    read(31)((x(i,k),k=0,N),i=0,M)
    read(31)((y(i,k),k=0,N),i=0,M)
    read(31)((u(i,k),k=0,N),i=0,M)
    read(31)((v(i,k),k=0,N),i=0,M)
    read(31)((p(i,k),k=0,N),i=0,M)
    read(31)((f(i,k),k=0,N),i=0,M)
    read(31)((c(i,k),k=0,N),i=0,M)
    read(31)t
    close(31)
    return
    end


    subroutine output_var(M,N,var,fname)
    implicit none
    integer M,N,i,k
    real var(1:M,1:N)
    character fname*7
    open(11,file=fname, status='unknown',form='unformatted')
    write(11)((var(i,k),k=1,N),i=1,M)
    close(11)
    return
    end

