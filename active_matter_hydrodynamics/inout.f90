    subroutine output(M,N,x,y,u,v,p,f,t,fname)
    implicit none
    integer M,N,i,k
    real x(0:M,0:N),y(0:M,0:N)
    real u(-1:M,-1:N),v(-1:M,-1:N)
    ! careful: here p is not the same format as u,v,p
    real p(-1:M,-1:N),f(-1:M,-1:N)
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


    subroutine output2(M,N,x,y,u,v,p,f,f_2, pf1, pf2, t,fname)
    implicit none
    integer M,N,i,k
    real x(0:M,0:N),y(0:M,0:N)
    real u(-1:M,-1:N),v(-1:M,-1:N)
    ! careful: here p is not the same format as u,v,p
    real p(-1:M,-1:N),f(-1:M,-1:N),f_2(-1:M,-1:N)
    real pf1(-1:M,-1:N), pf2(-1:M,-1:N)
    real t
    character fname*7
    open(11,file=fname, status='unknown',form='unformatted')
    write(11)((x(i,k),k=0,N),i=0,M)
    write(11)((y(i,k),k=0,N),i=0,M)
    write(11)((u(i,k),k=0,N),i=0,M)
    write(11)((v(i,k),k=0,N),i=0,M)
    write(11)((p(i,k),k=0,N),i=0,M)
    write(11)((f(i,k),k=0,N),i=0,M)
    write(11)((f_2(i,k),k=0,N),i=0,M)
    write(11)((pf1(i,k),k=0,N),i=0,M)
    write(11)((pf2(i,k),k=0,N),i=0,M)
    write(11)t
    close(11)
    return
    end



    subroutine output_fib_swm(np,pos, np_swm, pos_swm, swm_ang, t, fname)
        implicit none
        integer np, np_swm, i,k
        real pos(0:np-1,0:1), pos_swm(0:np_swm-1, 0:1), swm_ang(0:np_swm-1)
        real t
        character fname*7

        ! 打开文件
        open(unit=10, file=fname // '_fib.txt', status='replace', action='write')
        ! 写入数组
        do i = 0, np-1
            write(10, '(I10, 3F10.5)') i, pos(i,0), pos(i,1), t
        end do
        ! 关闭文件
        close(10)

        ! 打开文件
        open(unit=10, file=fname // '_swm.txt', status='replace', action='write')
        ! 写入数组
        do i = 0, np_swm-1
            write(10, '(I10, 3F10.5)') i, pos_swm(i,0), pos_swm(i,1), swm_ang(i)
        end do
        ! 关闭文件
        close(10)

        return
        end




    subroutine output_regular(M,N,x,y,u,v,p,f,t,fname)
    implicit none
    integer M,N,i,k
    real x(0:M,0:N),y(0:M,0:N)
    real u(-1:M,-1:N),v(-1:M,-1:N)
    real f(-1:M,-1:N),p(-1:M,-1:N)
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

    subroutine input(M,N,x,y,u,v,p,f,t,fname)
    implicit none
    integer M,N,i,k
    real x(0:M,0:N),y(0:M,0:N),u(-1:M+1,-1:N),v(-1:M,0:N)
    ! careful: here p is not the same format as u,v,p
    real p(-1:M,-1:N),f(-1:M,-1:N)
    real t
    character fname*7

    open(31,file=fname, status='old',form='unformatted')
    read(31)((x(i,k),k=0,N),i=0,M)
    read(31)((y(i,k),k=0,N),i=0,M)
    read(31)((u(i,k),k=0,N),i=0,M)
    read(31)((v(i,k),k=0,N),i=0,M)
    read(31)((p(i,k),k=0,N),i=0,M)
    read(31)((f(i,k),k=0,N),i=0,M)
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


    subroutine read_var(M,N,var,fname)
    implicit none
    integer M,N,i,k
    real var(1:M,1:N)
    character fname*40
    open(31,file=fname, status='old',form='unformatted')
    read(31)((var(i,k),k=1,N),i=1,M)
    close(31)
    return
    end

