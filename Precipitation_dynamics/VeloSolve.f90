  subroutine VeloSolve(psi,rmu_after,ph_after,t)
  use global
  implicit none
  integer i,k
  real psi(-1:nr,-1:nz)
  real ph_after(-1:nr+1,-1:nz+1),rmu_after(-1:nr+1,-1:nz+1),t

  real,external::force1,force2

  real w1(-1:nr+1,-1:nz+1),w2(-1:nr+1,-1:nz+1),w3(-1:nr+1,-1:nz+1),w4(-1:nr+1,-1:nz+1)
  real fu(-1:nr+1,-1:nz+1),fv(-1:nr+1,-1:nz+1)

  real temp(-1:nr+1,-1:nz+1),tempp(-1:nr+1,-1:nz+1)
  real temp_p(0:nr+1,0:nz+1)

!  real rhs((nr-1)*(nz-1)),ComputedSolution((nr-1)*(nz-1))

  real gT(0:nr),gB(0:nr),cf

  real rhs_s((nr+1)*(nz+1)),ComputedSolution_s((nr+1)*(nz+1)) !s for slip case
  real rhs((nr+1)*(nz-1)),ComputedSolution((nr+1)*(nz-1))

  external MatrixVector_u_s,PSolve_u_s
  external MatrixVector_u,PSolve_u

  external MatrixVector_u_s_CG,PSolve_u_s_CG
  external MatrixVector_u_CG,PSolve_u_CG
  integer ITER,INFO
  real work((nr+1)*(nz+1),7),RESID
!        call  dif1x(nr,nz,xa1,xb1,xc1,ph,w1)
!        call  dif1y(nr,nz,ya1,yb1,yc1,ph,w2)
!
!      do k=0,nz
!         do i=0,nr
!            fu(i,k)=ca*dt*rmu(i,k)*w1(i,k)
!            fv(i,k)=ca*dt*rmu(i,k)*w2(i,k)
!         enddo
!      enddo


        call  dif1x(nr,nz,xa1,xb1,xc1,ph_after,w1)
        call  dif1y(nr,nz,ya1,yb1,yc1,ph_after,w2)

      do k=0,nz
         do i=0,nr
            fu(i,k)=ca*rmu_after(i,k)*w1(i,k)
            fv(i,k)=ca*rmu_after(i,k)*w2(i,k)
         enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CALL  dif1x(nr,nz,xa1,xb1,xc1,u,w1)
    call  dif1y(nr,nz,ya1,yb1,yc1,u,w2)
    call  dif1x(nr,nz,xa1,xb1,xc1,v,w3)
    call  dif1y(nr,nz,ya1,yb1,yc1,v,w4)

    do k=0,nz
          do i=0,nr
            fu(i,k)=fu(i,k)-re*rho(i,k)*u(i,k)*w1(i,k)-re*rho(i,k)*v(i,k)*w2(i,k)
            fv(i,k)=fv(i,k)-re*rho(i,k)*u(i,k)*w3(i,k)-re*rho(i,k)*v(i,k)*w4(i,k)
          end do
    end do

    CALL  dif1x(nr,nz,xa1,xb1,xc1,eta,temp)
    call  dif1y(nr,nz,ya1,yb1,yc1,eta,tempp)

    fu = fu + (2.*temp*w1 + tempp*(w2+w3))
    fv = fv + (temp*(w2+w3) + 2.*tempp*w4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pressure related terms
    !call dif1x(nr,nz,xa1,xb1,xc1,p+psi,w3)
    !call dif1y(nr,nz,ya1,yb1,yc1,p+psi,w4)

    temp_p = p + psi
    do k=0,nz
      do i=0,nr
        w3(i,k) = (temp_p(i+1, k) - temp_p(i, k) + temp_p(i+1, k+1) - temp_p(i, k+1))/2.0/dr
        w4(i,k) = (temp_p(i, k+1) - temp_p(i, k) + temp_p(i+1, k+1) - temp_p(i+1, k))/2.0/dz
      enddo
    enddo

!!$OMP PARALLEL DO PRIVATE(k,i) DEFAULT(SHARED)
      do k=0,nz
         do i=0,nr
            fu(i,k)=fu(i,k) - w3(i,k) + u(i,k)*re*rho(i,k)/dt + rho(i,k)*ggx
            fv(i,k)=fv(i,k) - w4(i,k) + v(i,k)*re*rho(i,k)/dt + rho(i,k)*ggz
         end do
      enddo
!!$OMP END PARALLEL DO

! external force

!   do k=0,nz
!      do i=0,nr
!         fu(i,k) = fu(i,k) + force1( r(i,k),z(i,k),t )
!         fv(i,k) = fv(i,k) + force2( r(i,k),z(i,k),t )
!      enddo
!   enddo


fu = fu*dt/eta
fv = fv*dt/eta


if(xlsB >= 1.0E-8) then
!!!! the value of u at boundary Robin BCs

      do i=0,nr
         cf=sqrt(2.0)*cos(angleB(i)*pi/180.0)/3.0*pi/2.
         alphaB(i) = 0.5*(1.0+ph_after(i,0))/xlsB + 0.5*(1.0-ph_after(i,0))/(xlsB*slipratio)
         gB(i) =  uwB*alphaB(i) + ca*( I_thickness*(ph_after(i,-1)-ph_after(i,1))/(2.*dz)&
& - cf*cos(pi*ph_after(i,0)/2.) + alpha_s*(ph_after(i,0)-ph(i,0)) )*(ph_after(i+1,0)-ph_after(i-1,0))/(2.*dr)/eta(i,0)
      end do

      do i=0,nr
         cf=sqrt(2.0)*cos(angleT(i)*pi/180.0)/3.0*pi/2.
         alphaT(i) = 0.5*(1.0+ph_after(i,nz))/xlsT+0.5*(1.0-ph_after(i,nz))/(xlsT*slipratio)
         gT(i) = uwT*alphaT(i) + ca*( I_thickness*(ph_after(i,nz+1)-ph_after(i,nz-1))/(2.*dz)&
& - cf*cos(pi*ph_after(i,nz)/2.) + alpha_s*(ph_after(i,nz)-ph(i,nz)) )*(ph_after(i+1,nz)-ph_after(i-1,nz))/(2.*dr)/eta(i,nz)
      end do

     do i=0,nr
        fu(i,0) = fu(i,0) + 2.*dt*gB(i)/dz
        fu(i,nz) = fu(i,nz) + 2.*dt*gT(i)/dz
     enddo

!     call ESOLP_x(nr,nz,ua,ub,uc,fu,u,wsavexf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!    transform the RHS to a single vector by calling AtoB
!!!         initial guess for the solution using the value at last time step

    call Atob(nr,nz,fu,rhs_s)

    call Atob(nr,nz,u,ComputedSolution_s)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  call FFT based preconditioned GMRES method to solve linear system of PH
!!  here 150 is the restart number of GMRES Iteration

       ! call output(nr,nz,r,z,u,v,fu,fv,0.0,'testin5')
       ! print *,'output testin5'
!	pause
!	call fgmres((nr+1)*(nz+1),150,ComputedSolution_s,rhs_s,tol,MatrixVector_u_s,PSolve_u_s)

ITER = 500
RESID = tol
    call BICGSTAB((nr+1)*(nz+1),rhs_s,ComputedSolution_s,work,(nr+1)*(nz+1),ITER,RESID,MatrixVector_u_s_CG,PSolve_u_s_CG,INFO)
	write(21,*) 'u 1 BICGSTAB',ITER
	if(INFO .ne. 0) then
		print *,'error'
		stop
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! transform the solution back to matrix form

    call btoA(nr,nz,u,ComputedSolution_s)

!!!! ghost points for u
   do k=0,nz
       u(-1,k) = u(1,k)
       u(nr+1,k) = u(nr-1,k)
   enddo
     do i=0,nr
        u(i,-1) = u(i,1) + 2.*dz*(gB(i)-u(i,0)*alphaB(i))
        u(i,nz+1) = u(i,nz-1) + 2.*dz*(gT(i)-u(i,nz)*alphaT(i))
     enddo

else

!! ~~~~~~~~~~~~~~~~~~~~~~~~~ need to be changed
! need to change MatrixVector_u
!!!!  no-slip along T.B. boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!    transform the RHS to a single vector by calling AtoB
!!!         initial guess for the solution using the value at last time step

    call Atob_v(nr,nz,fu,rhs)

    call Atob_v(nr,nz,u,ComputedSolution)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  call FFT based preconditioned GMRES method to solve linear system of PH
!!  here 150 is the restart number of GMRES Iteration
    call fgmres_u((nr+1)*(nz-1),150,ComputedSolution,rhs,tol)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! transform the solution back to matrix form

    call btoA_v(nr,nz,u,ComputedSolution)

  do i=0,nr
     u(i,0) = uwB
     u(i,nz) = uwT
  enddo
  do k=0,nz
       u(-1,k) = u(nr-1,k)
       u(nr+1,k) = u(1,k)
  enddo

endif

!!!!!  discrete transform for velocity v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!    transform the RHS to a single vector by calling AtoB
!!!         initial guess for the solution using the value at last time step

    call Atob_v(nr,nz,fv,rhs)

    call Atob_v(nr,nz,v,ComputedSolution)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  call FFT based preconditioned GMRES method to solve linear system of PH
!!  here 150 is the restart number of GMRES Iteration
!    call fgmres((nr+1)*(nz-1),150,ComputedSolution,rhs,tol,MatrixVector_u,PSolve_u)

ITER = 500
RESID = tol
    call BICGSTAB((nr+1)*(nz-1),rhs,ComputedSolution,work,(nr+1)*(nz+1),ITER,RESID,MatrixVector_u_CG,PSolve_u_CG,INFO)
	write(21,*) 'v 1 BICGSTAB',ITER
	if(INFO .ne. 0) then
		print *,'error'
		stop
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! transform the solution back to matrix form

    call btoA_v(nr,nz,v,ComputedSolution)

	do i=0,nr
		v(i,-1)=v(i,1)
		v(i,nz+1)=v(i,nz-1)
	enddo

!!!!!!!!!!!!!!!!!!!!! no flux of incompressiblity condition   \partiakl_n (\grad \bf u) = 0 at BD
!u(-1,nz+1)=u(nr-1,nz+1)
!u(-1,-1)=u(nr-1,-1)
!u(nr+1,nz+1)=u(1,nz+1)
!u(nr+1,-1)=u(1,-1)
!
!if(xlsB >= 1.0E-8) then
!  do i=0,nr
!     k=0
!     v(i,k-1) = 2*v(i,k) - v(i,k+1) - dz/2*( (u(i+1,k+1)-u(i-1,k+1))/dr/2 - (u(i+1,k-1)-u(i-1,k-1))/2/dr )
!     k=nz
!     v(i,nz+1) = 2*v(i,k) - v(i,nz-1) - dz/2*( (u(i+1,k+1)-u(i-1,k+1))/dr/2 - (u(i+1,k-1)-u(i-1,k-1))/2/dr )
!  enddo
!else   !!!!!!!!!!!!added
!  do i=0,nr
!     k=0
!     v(i,k-1) = 2*v(i,k) - v(i,k+1) - dz*( 4*(u(i+1,k+1)-u(i-1,k+1)) - (u(i+1,k+2)-u(i-1,k+2)) - 3*(u(i+1,k)-u(i-1,k)) )/2/dr/2
!     k=nz
!     v(i,k+1) = 2*v(i,k) - v(i,k-1) - dz*( (u(i+1,k-2)-u(i-1,k-2)) - 4*(u(i+1,k-1)-u(i-1,k-1)) + 3*(u(i+1,k)-u(i-1,k))  )/2/dr/2
!  enddo
!endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!         for velocity v
!do i=0,nr
!   gB(i) = (u(i+1,0)-u(i-1,0))/2./dr
!   fv(i,0) = fv(i,0) + 2.*dt*gB(i)/dz
!
!   gT(i) = -(u(i+1,nz)-u(i-1,nz))/2./dr
!   fv(i,nz) = fv(i,nz) + 2.*dt*gT(i)/dz
!enddo

!     call DCT(nr,nz,fv,v,Ve_Lc,wsavexc,wsaveyc)
!    !! Top BD&Bottom BD
!
!    do i=0,nr
!       v(i,-1)=v(i,1) + 2.*dz*gB(i)
!       v(i,nz+1)=v(i,nz-1) + 2.*dz*gT(i)
!    enddo

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Matrix operator of linear system for CH equation
  subroutine  MatrixVector_u(xx,yy)
  use global
  implicit none
  real xx((nr+1)*(nz-1))
  real yy((nr+1)*(nz-1))
  real f1(-1:nr+1,-1:nz+1),f2(-1:nr+1,-1:nz+1)
  integer i,j


  call btoA_v(nr,nz,f1,xx)

  do j=0,nz
     f1(-1,j) = f1(1,j)
     f1(nr+1,j) = f1(nr-1,j)
  enddo

  do j=1,nz-1
     do i=0,nr
        f2(i,j) = f1(i,j)*re*rho(i,j)/eta(i,j) &
& - dt*( (f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)-2*f1(i,j))/dz2 )
      enddo
  enddo

    call Atob_v(nr,nz,f2,yy)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        FFT based preconditioner for GMRES

  subroutine PSolve_u(xx,yy)
  use global
  implicit none
  real xx((nr+1)*(nz-1)),yy((nr+1)*(nz-1))
  real f1(-1:nr+1,-1:nz+1),fvv(-1:nr+1,-1:nz-1),vv(-1:nr+1,-1:nz-1)
  integer i,k

!    call btoA_v(nr,nz,f1,xx)
!
!     do k=0,nz-2
!        do i=0,nr
!           fvv(i,k) = f1(i,k+1)
!        enddo
!     enddo
!
!     call ESOLP_x(nr,nz-2,va,vb,vc,fvv,vv,wsavexc)
!
!     do k=1,nz-1
!        do i=0,nr
!           f1(i,k) = vv(i,k-1)
!        enddo
!     enddo
!
!!!!!!!!!!!!!!!!!!!!!!  bcd for v
!    call Atob_v(nr,nz,f1,yy)

	yy=xx
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Matrix operator of linear system for CH equation
  subroutine  MatrixVector_u_s(xx,yy)
  use global
  implicit none
  real xx((nr+1)*(nz+1)),yy((nr+1)*(nz+1))
  real f1(-1:nr+1,-1:nz+1),f2(-1:nr+1,-1:nz+1)
  integer i,j


  call btoA(nr,nz,f1,xx)

  do i=0,nr
     f1(i,-1) = f1(i,1)
     f1(i,nz+1) = f1(i,nz-1)
  enddo
  do j=0,nz
     f1(-1,j) = f1(1,j)
     f1(nr+1,j) = f1(nr-1,j)
  enddo

  do j=0,nz
     do i=0,nr
        f2(i,j) = f1(i,j)*re*rho(i,j)/eta(i,j) &
& - dt*( (f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)-2*f1(i,j))/dz2 )
      enddo
  enddo

  do i=0,nr
  	j=0
	f2(i,j) = f2(i,j) + f1(i,j)*2*dt*alphaB(i)/dz
  	j=nz
	f2(i,j) = f2(i,j) + f1(i,j)*2*dt*alphaT(i)/dz
  enddo

    call Atob(nr,nz,f2,yy)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        FFT based preconditioner for GMRES

  subroutine PSolve_u_s(xx,yy)
  use global
  implicit none
  real xx((nr+1)*(nz+1)),yy((nr+1)*(nz+1))
  real f1(-1:nr+1,-1:nz+1)

    call btoA(nr,nz,f1,xx)

    call ESOLP_x(nr,nz,ua,ub,uc,f1,f1,wsavexc)

    call Atob(nr,nz,f1,yy)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Matrix operator of linear system for CH equation
  subroutine  MatrixVector_u_CG(alpha,xx,beta,yy)
  use global
  implicit none
  real xx((nr+1)*(nz-1)),yy((nr+1)*(nz-1)),alpha,beta
  real f1(-1:nr+1,-1:nz+1),f2(-1:nr+1,-1:nz+1)
  integer i,j


  call btoA_v(nr,nz,f1,xx)

  do j=0,nz
     f1(-1,j) = f1(1,j)
     f1(nr+1,j) = f1(nr-1,j)
  enddo

  do j=1,nz-1
     do i=0,nr
        f2(i,j) = f1(i,j)*re*rho(i,j)/eta(i,j) &
& - dt*( (f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)-2*f1(i,j))/dz2 )
      enddo
  enddo

  call btoA_v(nr,nz,f1,yy)

    call Atob_v(nr,nz,alpha*f2+beta*f1,yy)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        FFT based preconditioner for GMRES

  subroutine PSolve_u_CG(yy,xx)
  use global
  implicit none
  real xx((nr+1)*(nz-1)),yy((nr+1)*(nz-1))
  real f1(-1:nr+1,-1:nz+1),fvv(-1:nr+1,-1:nz-1),vv(-1:nr+1,-1:nz-1)
  integer i,k

    call btoA_v(nr,nz,f1,xx)

     do k=0,nz-2
        do i=0,nr
           fvv(i,k) = f1(i,k+1)
        enddo
     enddo

     call ESOLP_x(nr,nz-2,va,vb,vc,fvv,vv,wsavexc)

     do k=1,nz-1
        do i=0,nr
           f1(i,k) = vv(i,k-1)
        enddo
     enddo

!!!!!!!!!!!!!!!!!!!!!  bcd for v
    call Atob_v(nr,nz,f1,yy)

  end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Matrix operator of linear system for CH equation
  subroutine  MatrixVector_u_s_CG(alpha,xx,beta,yy)
  use global
  implicit none
  real xx((nr+1)*(nz+1)),yy((nr+1)*(nz+1)),alpha,beta
  real f1(-1:nr+1,-1:nz+1),f2(-1:nr+1,-1:nz+1)
  integer i,j


  call btoA(nr,nz,f1,xx)

  do j=0,nz
     f1(-1,j) = f1(1,j)
     f1(nr+1,j) = f1(nr-1,j)
  enddo
  do i=0,nr
     f1(i,-1) = f1(i,1)
     f1(i,nz+1) = f1(i,nz-1)
  enddo

  do j=0,nz
     do i=0,nr
        f2(i,j) = f1(i,j)*re*rho(i,j)/eta(i,j) &
& - dt*( (f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)-2*f1(i,j))/dz2 )
      enddo
  enddo

  do i=0,nr
  	j=0
	f2(i,j) = f2(i,j) + f1(i,j)*2*dt*alphaB(i)/dz
  	j=nz
	f2(i,j) = f2(i,j) + f1(i,j)*2*dt*alphaT(i)/dz
  enddo

  call btoA(nr,nz,f1,yy)

    call Atob(nr,nz,alpha*f2+beta*f1,yy)

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        FFT based preconditioner for GMRES

  subroutine PSolve_u_s_CG(yy,xx)
  use global
  implicit none
  real xx((nr+1)*(nz+1)),yy((nr+1)*(nz+1))
  real f1(-1:nr+1,-1:nz+1)

    call btoA(nr,nz,f1,xx)

    call ESOLP_x(nr,nz,ua,ub,uc,f1,f1,wsavexc)

    call Atob(nr,nz,f1,yy)

  end subroutine
