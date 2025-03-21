  subroutine PhSolve(rmu_after,ph_after,t)
  use global
  implicit none
  real ph_after(-1:nr,-1:nz),rmu_after(-1:nr,-1:nz),t

  integer i,k

  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
  real ff(-1:nr,-1:nz)

  real pnt(0:nr-1),pnb(0:nr-1),pnl(0:nz-1),pnr(0:nz-1)
  real cf, tmp

  real rhs(1:nr*nz),Computed_Solution(1:nr*nz)
  real,external::ch_bulk_energy

!  external MATRIX_VECTOR,PSolve

  external MATRIX_VECTOR_CG,PSolve_CG
  integer ITER,INFO
  real work(nr*nz,7),RESID

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  calcuating the RHS of CH equation

  f2 = ph

  !TODO: should be u\grad\phi
  do k=0,nz-1
    do i=0,nr-1
       f2(i,k) = f2(i,k) - dt*(u(i+1,k)+u(i,k))/2.0 * (ph(i+1,k)-ph(i-1,k))/2.0/dr&
         & - dt*(v(i,k+1)+v(i,k))/2.0 * (ph(i,k+1)-ph(i,k-1))/2.0/dz
    enddo
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do k=0,nz-1
     do i=0,nr-1
        f1(i,k) = s_implicit*ph(i,k)/I_thickness - ch_bulk_energy(ph(i,k), bdr_ph(i,k))/I_thickness &
          - ph(i,k)*(ph(i,k)-1.0)*bdr_ph_grad(i,k)*cos(angle*pi/180)/sqrt(2.0)
     enddo
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! calculate right hand side as -RHS1 + [(-\epsilon^2)\grad^2 + sI]RHS2
  call op_double_gradient_Neum(nr, nz, dr2, dz2, f2, ph_M_inner, ff)
  f1 = -f1 + s_implicit*f2/I_thickness - I_thickness*ff

!   do j=0,nz-1
!      do i=0,nr-1
!          f1(i,j) = -f1(i,j) + s_implicit*f2(i,j)/I_thickness - I_thickness*( (f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 &
! &+ (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
!      enddo
!   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  transform the RHS to a single vector by calling AtoB
!!! initial guess for the solution using the value at last time step

    call Atob(nr-1,nz-1,f1,rhs)
    call Atob(nr-1,nz-1,rmu,Computed_Solution)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  call FFT based preconditioned GMRES method to solve linear system of PH
!!  here 150 is the restart number of GMRES Iteration

!    call fgmres((nr+1)*(nz+1),150,Computed_Solution,rhs,tol,MATRIX_VECTOR,PSolve)
    
ITER = 500
RESID = tol
    call BICGSTAB(nr*nz,rhs,Computed_Solution,work,nr*nz,ITER,RESID,MATRIX_VECTOR_CG,PSolve_CG,INFO)
	write(*,*) '\phi 1 BICGSTAB',ITER
	if(INFO .ne. 0) then
		print *,'error',INFO
		!stop
	endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! transform the solution back to matrix form

    call btoA(nr-1,nz-1,rmu_after,Computed_Solution)

!!!!!!!!!!!!!!!!!!!! \partial_n \mu = 0
  do i=0,nr-1
     rmu_after(i,-1)=rmu_after(i,0)
     rmu_after(i,nz)=rmu_after(i,nz-1)
  enddo

  do k=-1,nz
      rmu_after(-1,k) = rmu_after(0,k)
      rmu_after(nr,k) = rmu_after(nr-1,k)
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \ph
  call op_double_gradient_Neum(nr, nz, dr2, dz2, rmu_after, ph_M_0, ff)
  ph_after = f2 + dt*xld*ff

!   do k=0,nz-1
!      do i=0,nr-1
!          ph_after(i,k) = f2(i,k) + dt*xld*( (rmu_after(i+1,k)+rmu_after(i-1,k)-2*rmu_after(i,k))/dr2 &
! &+ (rmu_after(i,k+1)+rmu_after(i,k-1)-2*rmu_after(i,k))/dz2 )
!      enddo
!   enddo


k=0
do i=0,nr-1
  ph_after(i,k-1) = ph_after(i,k)
enddo

k=nz-1
do i=0,nr-1
    ph_after(i,k+1) = ph_after(i,k)
enddo

do k=-1,nz
    ph_after(-1,k) = ph_after(0,k)
    ph_after(nr,k) = ph_after(nr-1,k)
enddo

  end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MATRIX_VECTOR ans PSolve are used by fgmres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Matrix operator of linear system for CH equation
!  subroutine  MATRIX_VECTOR(xx,yy)
!  use global
!  implicit none
!  real xx((nr+1)*(nz+1))
!  real yy((nr+1)*(nz+1))
!  real f1(-1:nr+1,-1:nz+1),f2(-1:nr+1,-1:nz+1)
!  real ff1(-1:nr+1,-1:nz+1),ff2(-1:nr+1,-1:nz+1)
!  integer i,j
!
!  call btoA(nr,nz,f1,xx)
!
!  do j=0,nz
!     f1(-1,j)=f1(1,j)
!     f1(nr+1,j)=f1(nr-1,j)
!  enddo
!
!  do i=0,nr
!     f1(i,-1)=f1(i,1)
!     f1(i,nz+1)=f1(i,nz-1)
!  enddo
!
!  do j=0,nz
!     do i=0,nr
!        ff1(i,j) = -dt*xld*( (f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)-2*f1(i,j))/dz2)
!      enddo
!  enddo
!
!  do j=0,nz
!     ff1(-1,j)=ff1(1,j)
!     ff1(nr+1,j)=ff1(nr-1,j)
!  enddo
!
!  do i=0,nr
!     ff1(i,-1)=ff1(i,1)
!     ff1(i,nz+1)=ff1(i,nz-1)
!  enddo
!
!
!  do j=0,nz
!     do i=0,nr
!        ff2(i,j) = s_implicit*ff1(i,j)/I_thickness - I_thickness*( (ff1(i+1,j)+ff1(i-1,j)-2*ff1(i,j))/dr2 &
!&+ (ff1(i,j+1)+ff1(i,j-1)-2*ff1(i,j))/dz2 )
!     enddo
!  enddo
!
!  do i=0,nr
!     j=0
!     ff2(i,j) = ff2(i,j) + 2*alpha_s*ff1(i,j)/dz
!     j=nz
!     ff2(i,j) = ff2(i,j) + 2*alpha_s*ff1(i,j)/dz
!  enddo
!
!
!f2 = f1
!
!  do i=0,nr
!     j=0
!     f2(i,j)= f2(i,j) - 2*xld/vs/dz*( (f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)&
!&-2*f1(i,j))/dz2 )
!     j=nz
!     f2(i,j)= f2(i,j) - 2*xld/vs/dz*( (f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)&
!&-2*f1(i,j))/dz2 )
!  enddo
!
!f1 = f2 + ff2
!
!    call Atob_2(nr,nz,f1,yy)
!
!
!! diag = (-I_thickness**2*(-2./dr2 - 2./dz2) + s_implicit)
!
!!do i=1,(nr+1)*(nz+1)
!!   yy(i) = yy(i)/diag
!!enddo
!
!! diag =(-dt*xld*(-2./dr2 - 2./dz2))
!
!!do i=(nr+1)*(nz+1)+1,NN
!!   yy(i) = yy(i)/diag
!!enddo
!
!  end subroutine
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        FFT based preconditioner for GMRES
!
!  subroutine PSolve(xx,yy)
!  use global
!  implicit none
!  real xx((nr+1)*(nz+1)),yy((nr+1)*(nz+1))
!  real f1(-1:nr+1,-1:nz+1),f2(-1:nr+1,-1:nz+1)
!
!
!  call btoA(nr,nz,f1,xx)
!
!  call DCT(nr,nz,f1,f2,Ph_Lcc,wsavexc,wsaveyc)
!
!  call Atob(nr,nz,f2,yy)
!
!  end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Matrix operator of linear system for CH equation
  subroutine  MATRIX_VECTOR_CG(alpha,xx,beta,yy)
  use global
  implicit none
  real xx(nr*nz),yy(nr*nz)
  real alpha,beta
  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
  real ff1(-1:nr,-1:nz),ff2(-1:nr,-1:nz)
  integer i,j

  call btoA(nr-1,nz-1,f1,xx)

  call op_double_gradient_Neum(nr, nz, dr2, dz2, f1, ph_M_0, ff1)
  ff1 = -dt*xld*ff1

  ! do j=0,nz-1
  !    do i=0,nr-1
  !       ff1(i,j) = -dt*xld*( (f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)-2*f1(i,j))/dz2)
  !     enddo
  ! enddo


  call op_double_gradient_Neum(nr, nz, dr2, dz2, ff1, ph_M_inner, ff2)
  ff2 = s_implicit*ff1/I_thickness - I_thickness*ff2

!   ! -\epsilon\laplace + s/\epsilon
!   do j=0,nz-1
!      do i=0,nr-1
!         ff2(i,j) = s_implicit*ff1(i,j)/I_thickness - I_thickness*( (ff1(i+1,j)+ff1(i-1,j)-2*ff1(i,j))/dr2 &
! &+ (ff1(i,j+1)+ff1(i,j-1)-2*ff1(i,j))/dz2 )
!      enddo
!   enddo


  f1 = f1 + ff2

  call btoA(nr-1,nz-1,f2,yy)

  call Atob(nr-1,nz-1,alpha*f1+beta*f2,yy)


! diag = (-I_thickness**2*(-2./dr2 - 2./dz2) + s_implicit)

!do i=1,(nr+1)*(nz+1)
!   yy(i) = yy(i)/diag
!enddo

! diag =(-dt*xld*(-2./dr2 - 2./dz2))

!do i=(nr+1)*(nz+1)+1,NN
!   yy(i) = yy(i)/diag
!enddo

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        FFT based preconditioner for GMRES
  subroutine PSolve_CG(yy,xx)
  use global
  implicit none
  real xx(nr*nz),yy(nr*nz)
  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)

  ! yy = xx
  ! return
  
  call btoA(nr-1,nz-1,f1,xx)

  !TODO: here just simply modify nr, may not correct
  ! instead of find half away pre-conditioner, direct use regular grids
  call DCT(nr-1,nz-1,f1,f2,Ph_Lcc,wsavexc_ph,wsaveyc_ph)

  call Atob(nr-1,nz-1,f2,yy)

  end subroutine

