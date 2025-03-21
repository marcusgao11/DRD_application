  subroutine PreSolve_2(plan,plan_back,psi, coef)
  use global
  implicit none
  integer*8 plan,plan_back
  real psi(-1:nr,-1:nz)
  real coef

  !real uin(-1:nr+1,-1:nz+1)
  !real w1(-1:nr+1,-1:nz+1),w2(-1:nr+1,-1:nz+1)
  !real pnb(0:nr),pnt(0:nr),pnl(0:nz),pnr(0:nz)

  integer i,k
  real uin(0:nr-1,0:nz-1)


  uin = 0.0
  do k=0,nz-1
     do i=0,nr-1
      uin(i,k) = ( u(i+1,k)-u(i,k) )/dr + ( v(i,k+1)-v(i,k) )/dz
     enddo
  enddo

  print *,'output SUM SUM', SUM(SUM(uin,1))/nr/nz


!TODO: check which is faster
!TODO: solve the linear system
   !call DCT(nr,nz,uin,psi,Pre_Lcc,wsavexc,wsaveyc)
   !psi(1:nr, 1:nz) = uin


!  ! test dct solver
!  do k=-1,nz
!     do i=-1,nr
!      p(i,k) = cos(20* pi*(i+0.5)*dr/ slen )*cos(10* pi*(k+0.5)*dz/shig)! + 1.0E-3
!     enddo
!  enddo
!
!  uin = 0.0
!  do k=0,nz-1
!     do i=0,nr-1
!      uin(i,k) = ( p(i+1,k) + p(i-1,k) - 2.0*p(i,k) )/dr2 +&
!       & ( p(i,k-1) + p(i,k+1) - 2.0*p(i,k))/dz2
!     enddo
!  enddo

  call dct_fftw(nr,nz, plan,plan_back, uin, psi, Pre_Lcc, 4.0*nr*nz)


!  print *, MAXVAL(MAXVAL(ABS(p(0:nr-1,0:nz-1) - psi(0:nr-1,0:nz-1) ),1))
!  print *, MINVAL(MINVAL(ABS(p(0:nr-1,0:nz-1) - psi(0:nr-1,0:nz-1) ),1))


!  ! for test
!  ph = 0.0
!  do k=0,nz-1
!     do i=0,nr-1
!      ph(i,k) = ( psi(i+1,k) + psi(i-1,k) - 2.0*psi(i,k) )/dr2 +&
!       & ( psi(i,k-1) + psi(i,k+1) - 2.0*psi(i,k))/dz2
!     enddo
!  enddo
!
!  print *, MAXVAL(MAXVAL(ABS(uin(0:nr-1,0:nz-1) - ph(0:nr-1,0:nz-1) ),1))
!  print *, MINVAL(MINVAL(ABS(uin(0:nr-1,0:nz-1) - ph(0:nr-1,0:nz-1) ),1))
!
!  call output(nr,nz,r,z,u,v,psi,ph,0,'testin2')
!  print *,'output testin2'



!   ! update velocity if needed
!   ! TODO: update velocity
!   do k=0,nz-1 ! Correct the u-velocity
!      do i=0,nr
!       u(i,k) = u(i,k) - (psi(i,k)-psi(i-1,k))/dr;
!      enddo
!   enddo

!   do k=1,nz-1 ! Correct the v-velocity
!     do i=0,nr-1
!       v(i,k) = v(i,k) - (psi(i,k)-psi(i,k-1))/dz;
!     enddo
!   enddo

!   ! apply boundary conditions
!   do k=0,nz-1
!      u(-1,k) = 2.0*uL(k) - u(1,k)
!   enddo

!   do k=0,nz-1
!      u(nr+1,k) = 2.0*uR(k) - u(nr-1,k)
!   enddo

!   do i=-1,nr+1
!      u(i,-1) = uwB*2.0 - u(i,0)
!      u(i,nz) = uwT*2.0 - u(i,nz-1)
!   enddo

!   do i=-1,nr
!      v(i,0) = 0.0
!      v(i,nz) = 0.0
!   enddo

!   do k=1,nz-1
!      v(-1,k) = v(0,k)
!      v(nr,k) = v(nr-1,k)
!   enddo

!  ! for test
!  ph = 0.0
!  do k=0,nz-1
!   do i=0,nr-1
!      ph(i,k) = ( u(i+1,k)-u(i,k) )/dr  + ( v(i,k+1)-v(i,k) )/dz
!   enddo
!  enddo

!  print *, "div", MAXVAL(MAXVAL(ABS(ph(0:nr-1,0:nz-1)),1))
!  print *,"div", MINVAL(MINVAL(ABS(ph(0:nr-1,0:nz-1)),1))
!
!  call output(nr,nz,r,z,u,v,psi,ph,0,'testin3')
!  print *,'output testin3'
!	pause

  ! re-scale of pressure
  psi = psi*coef

  ! update velocity if needed

  ! verify
  end subroutine


