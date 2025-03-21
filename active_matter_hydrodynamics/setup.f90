    subroutine setup()
    use global
    implicit none
    integer i,j,k
    real temp(0:nr,0:nz),t1
    real xlsB_ave, xlsT_ave

!!!!!coefficients for differentiation
    do i=0,nr
       xa1(i) = -1.0/2.0/dr
       xb1(i) = 0.0
       xc1(i) = -xa1(i)
    enddo

    do j=0,nz
       ya1(j) = -1.0/2.0/dz
       yb1(j) = 0.0
       yc1(j) = -ya1(j)
    enddo

    do i=0,nr
       xa2(i) = 1.0/dr2
       xb2(i) = -2.0*xa2(i)
       xc2(i) = xa2(i)
    enddo

    do j=0,nz
       ya2(j) = 1.0/dz2
       yb2(j) = -2.0*ya2(j)
       yc2(j) = ya2(j)
    enddo

!!!!!!!!!!!!!!!!!!!!!!   setup FFT staff for periodic BCs

!!!!!!!!!!!! in x direction

  i=0
     t1 = pi*i/nr
     fmx(i) = 4.0*sin(t1)**2/dr2
  do i=1,nr/2
     t1 = pi*i/nr
     fmx(i*2-1) = 4.0*sin(t1)**2/dr2
     fmx(i*2)=fmx(i*2-1)
  enddo

  call VRFFTI(nr,wsavexf)

!!!!!!!!!!! in y direction
   j=0
   t1 = pi*j/nz
   fmy(j) = 4.0*sin(t1)**2/dz2

   do j=1,nz/2
       t1 = pi*j/nz
       fmy(j*2-1) = 4.0*sin(t1)**2/dz2
       fmy(j*2) = fmy(j*2-1)
   enddo

  call VRFFTI(nz,wsaveyf)

  do j=0,nz
     do i=0,nr
        temp(i,j)=fmx(i)+fmy(j)
      !   Ph_Lpp(i,j)= 1. + xld*dt*s_implicit*temp(i,j) + xld*I_thickness**2*dt*temp(i,j)**2
      !   Ph_Lpp_2(i,j)= 1.5 + xld*dt*s_implicit*temp(i,j) + xld*I_thickness**2*dt*temp(i,j)**2
      !   Ve_Lpp(i,j)= re + dt*temp(i,j)
      !   Ve_Lpp_2(i,j)= 1.5*re + dt*temp(i,j)
        Pre_Lpp(i,j) = -temp(i,j)
     enddo
  enddo
  Pre_Lpp(0,0) = -1.0


!!!!!!!!!!!!!!!!!!!!!!!   setup CFT staff for periodic in x direction
!!! \partial_n \phi=0 in z direction
!
!   call VRFFTI(nr,wsavexf)
!
!   call VCOSTI(nz+1,wsaveyc)
!
!!!!!!!!!!!!! in x direction
!   i=0
!      t1 = pi*i/nr
!      fmx(i) = 4.0*sin(t1)**2/dr2
!   do i=1,nr/2
!      t1 = pi*i/nr
!      fmx(i*2-1) = 4.0*sin(t1)**2/dr2
!      fmx(i*2)=fmx(i*2-1)
!   enddo
!
!!!!!!!!!!!! in y direction
!   do j=0,nz
!      t1 = pi*j/nz/2.
!      fmy(j) = 4.0*sin(t1)**2/dz2
!   enddo
!
!   do j=0,nz
!      do i=0,nr
!         temp(i,j)=fmx(i)+fmy(j)
!         Ph_Lcf(i,j)= 1. + xld*dt*s_implicit*temp(i,j)/I_thickness + xld*I_thickness*dt*temp(i,j)**2
!         Ph_Lcf_2(i,j)= 1.5 + xld*dt*s_implicit*temp(i,j)/I_thickness + xld*I_thickness*dt*temp(i,j)**2
!      enddo
!   enddo



!!!!!!!! setup for costine transform in two directions
!
   do i=0,nr
      t1 = pi*i/2.0/nr
      fmx(i) = 4.0*sin(t1)**2/dr2
   enddo

   call VCOSTI(nr+1,wsavexc)

   do j=0,nz
      t1 = pi*j/2.0/nz
      fmy(j) = 4.0*sin(t1)**2/dz2
   enddo

   call VCOSTI(nz+1,wsaveyc)

   do j=0,nz
      do i=0,nr
         temp(i,j)=fmx(i)+fmy(j)
         Ve_Lcc(i,j)= re + dt*temp(i,j)
         Ve_Lcc_2(i,j)= re*1.5 + dt*temp(i,j)
         Pre_Lcc(i,j) = -temp(i,j)
      enddo
   enddo

	Pre_Lcc(0,0) = -1.0


  ! setup for phi
  do i=0,nr-1
      t1 = pi*i/2.0/(nr-1)
      fmx_ph(i) = 4.0*sin(t1)**2/dr2
   enddo

   call VCOSTI(nr,wsavexc_ph)

   do j=0,nz-1
      t1 = pi*j/2.0/(nz-1)
      fmy_ph(j) = 4.0*sin(t1)**2/dz2
   enddo

   call VCOSTI(nz,wsaveyc_ph)

   do j=0,nz-1
      do i=0,nr-1
         temp(i,j)=fmx_ph(i)+fmy_ph(j)
         Ph_Lcc(i,j)= 1. + xld*dt*s_implicit*temp(i,j)/I_thickness + xld*I_thickness*dt*temp(i,j)**2
         Ph_Lcc_2(i,j)= 1.5 + xld*dt*s_implicit*temp(i,j)/I_thickness + xld*I_thickness*dt*temp(i,j)**2
      enddo
   enddo

!!!!!!!!!!!! discrete cosine transform only in x direction, and Robin BC in the other direction

   t1 = re*(1.0 + lambda_rho) / (1.0 + lambda_eta)
   xlsB_ave = xlsB*(1.0+slipratio)/2.0
   xlsT_ave = xlsT*(1.0+slipratio)/2.0

   xlsB_ave = (2.0*xlsB_ave - dz) / (2.0*xlsB_ave + dz)
   xlsT_ave = (2.0*xlsT_ave - dz) / (2.0*xlsT_ave + dz)

   do i=0,nr
      t1 = pi*i/nr/2.
      fmx(i) = 4.0*dt*sin(t1)**2/dr2
   enddo

   do k=0,nr
      i=1
      ua(i,k)=0.0
      ub(i,k)= t1 + ( (2.0-xlsB_ave) *dt/dz2 + fmx(k)) ! +2.*dt/dz/xlsB_ave)
      uc(i,k)=(-1.0*dt/dz2)

      ua_2(i,k)=0.0
      ub_2(i,k)= t1*1.5 + ( (2.0-xlsB_ave) *dt/dz2 + fmx(k)) ! +2.*dt/dz/xlsB_ave)
      uc_2(i,k)=(-1.0*dt/dz2)

      do i=2,nz-1
          ua(i,k)=(-dt/dz2)
          ub(i,k)= t1 + (2.0*dt/dz2 + fmx(k))
          uc(i,k)=(-dt/dz2)

          ua_2(i,k)=(-dt/dz2)
          ub_2(i,k)= t1*1.5 + (2.0*dt/dz2 + fmx(k))
          uc_2(i,k)=(-dt/dz2)
      end do

      i=nz
      ua(i,k)=(-1.0*dt/dz2)
      ub(i,k)= t1 + ( (2.0-xlsT_ave) *dt/dz2 + fmx(k)) ! +2.*dt/dz/xlsT_ave)
      uc(i,k)=0.0

      ua_2(i,k)=(-1.0*dt/dz2)
      ub_2(i,k)= t1*1.5 + ( (2.0-xlsT_ave) *dt/dz2 + fmx(k)) ! +2.*dt/dz/xlsT_ave)
      uc_2(i,k)=0.0
   end do

!!!!!!!!!!!! discrete cosine transform only in x direction, and Dirichlet BC in the other direction
   do i=0,nr-1
      t1 = pi*i/(nr-1)/2.
      fmx(i) = 4.0*dt*sin(t1)**2/dr2
   enddo
   call VCOSTI(nr,wsavexc_v)

   do k=0,nr-1
      i=1
      va(i,k)=0.0
      vb(i,k)= t1 + 2.0*dt/dz2 + fmx(k)
      vc(i,k)=-1.0*dt/dz2

      va_2(i,k)=0.0
      vb_2(i,k)= t1*1.5 + 2.0*dt/dz2 + fmx(k)
      vc_2(i,k)=-1.0*dt/dz2
      do i=2,nz-2
          va(i,k)=-dt/dz2
          vb(i,k)= t1 + 2.0*dt/dz2 + fmx(k)
          vc(i,k)=-dt/dz2

          va_2(i,k)=-dt/dz2
          vb_2(i,k)= t1*1.5 + 2.0*dt/dz2 + fmx(k)
          vc_2(i,k)=-dt/dz2
      end do

      i=nz-1
      va(i,k)=-1.0*dt/dz2
      vb(i,k)= t1 + 2.0*dt/dz2 + fmx(k)
      vc(i,k)=0.0

      va_2(i,k)=-1.0*dt/dz2
      vb_2(i,k)= t1*1.5 + 2.0*dt/dz2 + fmx(k)
      vc_2(i,k)=0.0
  end do


!!!!!!!!!!!!!!!!!!!!!!!   setup for Splittine_PhSolve
!
!t1 = 2*I_thickness/sqrt(xld*dt)
!
!s_implicit = max(t1,s_implicit)
!
!alpha_s = ( -s_implicit + sqrt(s_implicit**2-t1**2) )/2
!
!
!!!!!!!! setup for costine transform in two directions
!
!   do i=0,nr
!      t1 = pi*i/2.0/nr
!      fmx(i) = 4.0*sin(t1)**2/dr2
!   enddo
!
!   do j=0,nz
!      t1 = pi*j/2.0/nz
!      fmy(j) = 4.0*sin(t1)**2/dz2
!   enddo
!
!   do j=0,nz
!      do i=0,nr
!         temp(i,j)=fmx(i)+fmy(j)
!         rmu_s_Lc(i,j) = temp(i,j)*xld - alpha_s*xld/(I_thickness**2)
!         Ph_s_Lc(i,j) = temp(i,j)*(I_thickness**2) + alpha_s + s_implicit
!      enddo
!   enddo

    end subroutine

    subroutine dif1x(M,N,ua,ub,uc,u,du)   !!can be used for dif1x dif2x
    implicit none
    real ua(0:M),ub(0:M),uc(0:M)
    real u(-1:M+1,-1:N+1),du(-1:M+1,-1:N+1)
    integer i,j,M,N

    do j=0,N
       do i=0,M
          du(i,j) = ua(i)*u(i-1,j) + ub(i)*u(i,j) + uc(i)*u(i+1,j)
       enddo
    enddo

    end subroutine

    subroutine dif2x(M,N,ua,ub,uc,u,du)   !!can be used for dif1x dif2x
    implicit none
    real ua(0:M),ub(0:M),uc(0:M)
    real u(-1:M+1,-1:N+1),du(-1:M+1,-1:N+1)
    integer i,j,M,N

    do j=0,N
       do i=0,M
          du(i,j) = ua(i)*u(i-1,j) + ub(i)*u(i,j) + uc(i)*u(i+1,j)
       enddo
    enddo

    end subroutine

    subroutine dif1y(M,N,ua,ub,uc,u,du)  !!can be used for dif1y dif2y
    implicit none
    real ua(0:N),ub(0:N),uc(0:N)
    real u(-1:M+1,-1:N+1),du(-1:M+1,-1:N+1)
    integer i,j,M,N

    do j=0,N
       do i=0,M
          du(i,j) = ua(j)*u(i,j-1) + ub(j)*u(i,j) + uc(j)*u(i,j+1)
       enddo
    enddo

    end subroutine

    subroutine dif2y(M,N,ua,ub,uc,u,du)  !!can be used for dif1y dif2y
    implicit none
    real ua(0:N),ub(0:N),uc(0:N)
    real u(-1:M+1,-1:N+1),du(-1:M+1,-1:N+1)
    integer i,j,M,N

    do j=0,N
       do i=0,M
          du(i,j) = ua(j)*u(i,j-1) + ub(j)*u(i,j) + uc(j)*u(i,j+1)
       enddo
    enddo

    end subroutine

!!!! Regular tridiagonal solver

       subroutine rfords(n,ff,aa,bb,cc)
       implicit none
       integer n,j
       real ff(1:n),aa(1:n),bb(1:n),cc(1:n)
       real f(1:n),a(1:n),b(1:n),c(1:n)
       real r
       a=aa
       b=bb
       c=cc
       f=ff

        do j=2,n
           r = - a(j)/b(j-1)
           b(j) = b(j) + r*c(j-1)
           f(j) = f(j) + r*f(j-1)
       enddo

       if(abs(b(n)) <= 1.0E-12) then
          ff(n)=0.0
       else
          ff(n)=f(n)/b(n)
       endif

       do j=n-1,1,-1
          ff(j) = ( f(j) - c(j)*ff(j+1) )/b(j)
       enddo
       return
       end subroutine
