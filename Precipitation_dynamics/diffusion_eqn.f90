! diffusion_eqn for concentration
! first order implicit scheme

  subroutine diffusion_eqn(c_after,t)
  use global
  implicit none
  real c_after(-1:nr,-1:nz),t

  integer i,k

  real f2(-1:nr,-1:nz)

  real cf, tmp, moL, moR, coef

  real rhs(1:nr*nz),Computed_Solution(1:nr*nz)

!  external MATRIX_VECTOR,PSolve

  external MATRIX_VECTOR_CG_C,PSolve_CG_C
  integer ITER,INFO
  real work(nr*nz,7),RESID


shift_delta = 1.0E-3

! initialize
! TODO: move this
coef = 0.00 !01
ph_trunc = ph
  do k=-1,nz
      do i=-1,nr
         if(ph(i,k) <= -1.0) then
            !ph_trunc(i,k) = -1.0
            D_inner(i,k) = 1.0
         elseif(ph(i,k) >= 1.0) then
            !ph_trunc(i,k) = 1.0
            D_inner(i,k) = coef
         else
          D_inner(i,k) = coef + (1.0 - coef) * (1.0 - ph(i,k))/2.0
         endif
      enddo
  enddo


!D_inner = 1.0 + (0.01 - 1.0) * (1.0 - ph)/2.0
!D_inner = 1.0 - (1.0 - ph)/2.0
!D_inner = D_inner + shift_delta


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  calcuating the RHS of CH equation

  f2 = conc

  !TODO: should be u\grad\phi
  do k=0,nz-1
    do i=0,nr-1
       f2(i,k) = f2(i,k) - dt*(u(i+1,k)+u(i,k))/2.0 * (conc(i+1,k)-conc(i-1,k))/2.0/dr&
         & - dt*(v(i,k+1)+v(i,k))/2.0 * (conc(i,k+1)-conc(i,k-1))/2.0/dz
    enddo
  enddo

 ! ph_grad for further use
 ! boundary condition
  do k=0,nz-1
    do i=0,nr-1

!      if (conc(i,k) >= C_0) then
!          tmp = conc(i,k)**2/C_0**2 - 1.0
!      else
!          tmp = 0.0
!      endif

       tmp = conc(i,k)**2/C_0**2 - 1.0

       f2(i,k) = f2(i,k) - dt*tmp*(1.0 - conc(i,k)) * ph_grad(i,k) !/ ( ( 1.0 + ph_trunc(i,k))/2.0 + shift_delta )

       !f2(i,k) = f2(i,k) - dt*1.0 * ph_grad(i,k) / ( ( 1.0 + ph_trunc(i,k))/2.0 + 1.0E-6 )
    enddo
  enddo


!  ! Dirichlet on top and bottom
!  do i=0,nr-1
!     k = 0
!
!     moL = (D_inner(i,k-1) + D_inner(i,k))/2.0
!
!     f2(i,k) = f2(i,k) + dt*2.0*D_0*moL*c_bot/dz2  !/ ( ( 1.0 + ph_trunc(i,k))/2.0 + shift_delta )
!
!     k = nz - 1
!
!     moR = (D_inner(i,k) + D_inner(i,k+1))/2.0
!     f2(i,k) = f2(i,k) + dt*2.0*D_0*moR*c_top/dz2  !/ ( ( 1.0 + ph_trunc(i,k))/2.0 + shift_delta )
!  enddo


  ! Dirichlet on left and homogenernous at right
  do k=0,nz-1
     i = 0

     moL = (D_inner(i-1,k) + D_inner(i,k))/2.0

     f2(i,k) = f2(i,k) + dt*2.0*D_0*moL*c_left/dr2  !/ ( ( 1.0 + ph_trunc(i,k))/2.0 + shift_delta )

     !i = nr - 1
     !moR = (D_inner(i+1,k) + D_inner(i,k))/2.0
     !f2(i,k) = f2(i,k)  ! + dt*2.0*D_0*moR*c_top/dz2  !/ ( ( 1.0 + ph_trunc(i,k))/2.0 + shift_delta )
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  transform the RHS to a single vector by calling AtoB
!!! initial guess for the solution using the value at last time step

    call Atob(nr-1,nz-1,f2,rhs)
    call Atob(nr-1,nz-1,conc,Computed_Solution)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  call FFT based preconditioned GMRES method to solve linear system of PH
!!  here 150 is the restart number of GMRES Iteration

!    call fgmres((nr+1)*(nz+1),150,Computed_Solution,rhs,tol,MATRIX_VECTOR,PSolve)

ITER = 500
RESID = tol
    call BICGSTAB(nr*nz,rhs,Computed_Solution,work,nr*nz,ITER,RESID,MATRIX_VECTOR_CG_C,PSolve_CG_C,INFO)
	write(*,*) '\c 1 BICGSTAB',ITER
	if(INFO .ne. 0) then
		print *,'error',INFO
		!stop
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! transform the solution back to matrix form
    call btoA(nr-1,nz-1,c_after,Computed_Solution)

! boundary condition
!  do i=0,nr-1
!     c_after(i,-1) = 2.0*c_bot - c_after(i,0)
!     c_after(i,nz) = 2.0*c_top - c_after(i, nz-1)
!  enddo

! periodic at t and b
  do i=0,nr-1
     c_after(i,-1) = c_after(i, nz-1)
     c_after(i,nz) = c_after(i, 0)
  enddo


  do k=-1,nz
     c_after(-1,k) = 2.0*c_left - c_after(0,k)
     c_after(nr,k) = c_after(nr-1,k)
  enddo


  end subroutine



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Matrix operator of linear system for CH equation
  subroutine  MATRIX_VECTOR_CG_C(alpha,xx,beta,yy)
  use global
  implicit none
  real xx(nr*nz),yy(nr*nz)
  real alpha,beta
  real ff1(-1:nr,-1:nz),ff2(-1:nr,-1:nz)


    call btoA(nr-1,nz-1,ff1,xx)

    !call laplace_nd(ff1,ff2,nr,nz,dr2,dz2, c_bot, c_top)
    !call laplace_nd(ff1,ff2,nr,nz,dr2,dz2, 0.0, 0.0)


    !call laplace_nd_var_coef(ff1,ff2,nr,nz,dr2,dz2, D_inner)
    call laplace_dp_var_coef(ff1,ff2,nr,nz,dr2,dz2, D_inner)

    ff2 = ff1 - dt*D_0*ff2  !/ ( ( 1.0 + ph_trunc)/2.0 + shift_delta )


    call btoA(nr-1,nz-1,ff1,yy) !yy not changed

    call Atob(nr-1,nz-1, alpha*ff2 + beta*ff1,yy)


  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        FFT based preconditioner for GMRES
  subroutine PSolve_CG_C(yy,xx)
  use global
  implicit none
  real xx(nr*nz),yy(nr*nz)
  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)

  yy = xx


!  call btoA(nr-1,nz-1,f1,xx)
!
!  !TODO: here just simply modify nr, may not correct
!  ! instead of find half away pre-conditioner, direct use regular grids
!  call DCT(nr-1,nz-1,f1,f2,Ph_Lcc,wsavexc_ph,wsaveyc_ph)
!
!  call Atob(nr-1,nz-1,f2,yy)

  end subroutine

