  subroutine PhSolve(rmu_after,ph_after,t)
  use global
  implicit none
  real ph_after(-1:nr,-1:nz),rmu_after(-1:nr,-1:nz),t

  integer i,k

  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
  real ff(-1:nr,-1:nz)

  real pnt(0:nr-1),pnb(0:nr-1),pnl(0:nz-1),pnr(0:nz-1)
  real cf, tmp, tmp1

  real rhs(1:nr*nz),Computed_Solution(1:nr*nz)

!  external MATRIX_VECTOR,PSolve

  external MATRIX_VECTOR_CG,PSolve_CG
  integer ITER,INFO
  real work(nr*nz,7),RESID

!  ph_after = ph
!  return

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   !
!!!!!  calcuating the RHS of CH equation

  f2 = ph

!  !TODO: should be u\grad\phi
!  do k=0,nz-1
!    do i=0,nr-1
!       f2(i,k) = f2(i,k) - dt*(u(i+1,k)+u(i,k))/2.0 * (ph(i+1,k)-ph(i-1,k))/2.0/dr&
!         & - dt*(v(i,k+1)+v(i,k))/2.0 * (ph(i,k+1)-ph(i,k-1))/2.0/dz
!    enddo
!  enddo



  do k=0,nz-1
    do i=0,nr-1
      if (conc(i,k) >= C_0) then
          tmp = conc(i,k)**2/C_0**2 - 1.0
      else
          tmp = 0.0
      endif

      !tmp = 0.1

      f2(i,k) = f2(i,k) + dt*tmp*ph_grad(i,k)*2.0  ! in order to make ph local * (1.0-ph(i,k)**2.0)
    enddo
  enddo



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do k=0,nz-1
     do i=0,nr-1
        f1(i,k) = ((1.+s_implicit)*ph(i,k) - ph(i,k)**3)/I_thickness
     enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! boundary TODO:
! first Neuman

  do i=0,nr-1
     k = 0
     cf = sqrt(2.0)*cos(angleB(i)*pi/180.0)/3.0*pi/2.
     tmp = (ph(i,k-1) + ph(i,k))/2.0
     pnb(i) = cf*cos(tmp*pi/2) + alpha_s*tmp
  enddo


if (bot_type_ph== 'd') then
  do i=0,nr-1
     k = 0

     tmp1 = (ph(i,k) + ph(i,k-1))/2.0

     cf = tmp1/dt - (u(i+1,k)+u(i,k) + u(i+1,k-1)+u(i,k-1))/4.0&
       *((ph(i+1,k-1) + ph(i+1,k))/2.0 - (ph(i-1,k-1) + ph(i-1,k))/2.0 )/2.0/dr

     pnb(i) = pnb(i) + cf/vs

     tmp1 = I_thickness/dz + 0.5* ( 1.0/dt/vs + alpha_s )
     coef_bot = (I_thickness/dz - 0.5* ( 1.0/dt/vs + alpha_s ) ) / tmp1

     pnb(i) = pnb(i)/tmp1

     f1(i,k) = f1(i,k) + pnb(i)*I_thickness/dz2
  enddo
else if (bot_type_ph== 'p') then
  continue

else
  do i=0,nr-1
     k = 0

     tmp1 = I_thickness/dz + 0.5* alpha_s
     coef_bot = (I_thickness/dz - 0.5* alpha_s) / tmp1
     pnb(i) = pnb(i)/tmp1

     f1(i,k) = f1(i,k) + pnb(i)*I_thickness/dz2
  enddo
endif


  do i=0,nr-1
     k = nz-1
     cf = sqrt(2.0)*cos(angleT(i)*pi/180.0)/3.0*pi/2.
     tmp = (ph(i,k) + ph(i,k+1))/2.0
     pnt(i) = cf*cos(tmp*pi/2) + alpha_s*tmp
  enddo


if (top_type_ph  == 'd') then
  do i=0,nr-1
     k = nz-1
     tmp1 = (ph(i,k) + ph(i,k+1))/2.0
     cf = tmp1/dt - (u(i+1,k)+u(i,k) + u(i+1,k+1)+u(i,k+1))/4.0&
       *((ph(i+1,k+1) + ph(i+1,k))/2.0 - (ph(i-1,k+1) + ph(i-1,k))/2.0 )/2.0/dr

     pnt(i) = pnt(i) + cf/vs

     tmp1 = I_thickness/dz + 0.5* ( 1.0/dt/vs + alpha_s )
     coef_top = (I_thickness/dz - 0.5* ( 1.0/dt/vs + alpha_s ) ) / tmp1

     pnt(i) = pnt(i)/tmp1

     f1(i,k) = f1(i,k) + pnt(i)*I_thickness/dz2
  enddo
else if (top_type_ph== 'p') then
  continue

else
  do i=0,nr-1
     k = nz-1

     tmp1 = I_thickness/dz + 0.5* alpha_s
     coef_top = (I_thickness/dz - 0.5* alpha_s) / tmp1

     pnt(i) = pnt(i)/tmp1

     f1(i,k) = f1(i,k) + pnt(i)*I_thickness/dz2
  enddo
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! -------------------------
if (top_type_ph== 'p') then
    call laplace_np(f1,ff,nr,nz,dr2,dz2)

else
    if (lr_bd_type == 'p') then
    call laplace_pn(f1,ff,nr,nz,dr2,dz2)
  else
    call laplace_nn(f1,ff,nr,nz,dr2,dz2)
  endif
endif

  f2 = f2 - dt*xld*ff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!  transform the RHS to a single vector by calling AtoB
!!! initial guess for the solution using the value at last time step

    call Atob(nr-1,nz-1,f2,rhs)
    call Atob(nr-1,nz-1,ph,Computed_Solution)


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

    call btoA(nr-1,nz-1,ph_after,Computed_Solution)

!!!!!!!!!!!!!!!!!!!! \partial_n \mu = 0
if (top_type_ph== 'p') then
    do i=0,nr-1
     ph_after(i,-1) = ph_after(i,nz-1)
     ph_after(i,nz) = ph_after(i,0)
    enddo
else
    do i=0,nr-1
     ph_after(i,-1) = ph_after(i,0) * coef_bot + pnb(i)
     ph_after(i,nz) = ph_after(i,nz-1) * coef_top + pnt(i)
    enddo
endif


  if (lr_bd_type == 'p') then
    do k=-1,nz
      ph_after(-1,k) = ph_after(nr-1,k)
      ph_after(nr,k) = ph_after(0,k)
    enddo
  else
    do k=-1,nz
       ph_after(-1,k) = ph_after(0,k)
       ph_after(nr,k) = ph_after(nr-1,k)
    enddo
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \ph
! chemical potential, here can be changed numerical approximation
  do k=0,nz-1
     do i=0,nr-1
        rmu_after(i,k) = -I_thickness*( (ph_after(i+1,k)+ph_after(i-1,k)-2*ph_after(i,k))/dr2 &
          + (ph_after(i,k+1)+ph_after(i,k-1)-2*ph_after(i,k))/dz2 )&
         & - ph_after(i,k)/I_thickness + ph_after(i,k)**3/I_thickness
     enddo
  enddo

if (top_type_ph== 'p') then
    do i=0,nr-1
     rmu_after(i,-1) = rmu_after(i,nz-1)
     rmu_after(i,nz) = rmu_after(i,0)
    enddo
else
  do i=0,nr-1
     rmu_after(i,-1) = rmu_after(i,0)
     rmu_after(i,nz) = rmu_after(i,nz-1)
  enddo
endif


  if (lr_bd_type == 'p') then
    do k=-1,nz
       rmu_after(-1,k) = rmu_after(nr-1,k)
       rmu_after(nr,k) = rmu_after(0,k)
    enddo
  else
    do k=-1,nz
       rmu_after(-1,k) = rmu_after(0,k)
       rmu_after(nr,k) = rmu_after(nr-1,k)
    enddo
  endif

!!!!!!!!!!!!! for chemical potential
!    call dif2x(nr,nz,xa2,xb2,xc2,ph_after,w1)
!    call dif2y(nr,nz,ya2,yb2,yc2,ph_after,w2)
!
!    do j=0,nz
!       do i=0,nr
!          rmu_after(i,j) = -I_thickness*w1(i,j)-I_thickness*w2(i,j) + (-ph_after(i,j)+ph_after(i,j)**3)/I_thickness
!       enddo
!    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! BCDs for \partial_n \mu =0
!
!   do i=0,nr
!       rmu_after(i,-1) = rmu_after(i,1)
!       rmu_after(i,nz+1) = rmu_after(i,nz-1)
!   enddo
!
!   do j=0,nz
!      rmu_after(-1,j) = rmu_after(1,j)
!      rmu_after(nr+1,j) = rmu_after(nr-1,j)
!   enddo


!        call output(nr,nz,r,z,f1,f2,ph_after,rmu_after,0.0,'testin5')
!        print *,'output testin5'

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
  real f1(-1:nr,-1:nz)
  real ff1(-1:nr,-1:nz),ff2(-1:nr,-1:nz)
  integer i,j

  call btoA(nr-1,nz-1,f1,xx)

  ! -\dt \ld\laplace\
if (top_type_ph== 'p') then
    call laplace_np(f1,ff1,nr,nz,dr2,dz2)

else
  if (lr_bd_type == 'p') then
    call laplace_pn(f1,ff1,nr,nz,dr2,dz2)
  else
    call laplace_nn(f1,ff1,nr,nz,dr2,dz2)
  endif

endif



  ff1 = -I_thickness*ff1

  ! modify operators to satisfy non-homo boundary conditions
  if (top_type_ph== 'p') then
    continue
  else
    do i=0,nr-1
       j=0
       ff1(i,j) = ff1(i,j) + I_thickness/dz2*( f1(i,j-1) - coef_bot*f1(i,j) )

       j=nz-1
       ff1(i,j) = ff1(i,j) + I_thickness/dz2*( f1(i,j+1) - coef_top*f1(i,j) )
    enddo
  endif

  ! s/\epsilon
  ff1 = ff1 + f1*s_implicit/I_thickness

if (top_type_ph== 'p') then
    call laplace_np(ff1,ff2,nr,nz,dr2,dz2)

else
  ! -\dt\xld\laplace
  if (lr_bd_type == 'p') then
    call laplace_pn(ff1,ff2,nr,nz,dr2,dz2)
  else
    call laplace_nn(ff1,ff2,nr,nz,dr2,dz2)
  endif
endif


  ff2 = -dt*xld*ff2

  ! \beta
  ff2 = f1 + ff2


  call btoA(nr-1,nz-1,ff1,yy) !yy not changed

  call Atob(nr-1,nz-1, alpha*ff2 + beta*ff1,yy)


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


  call btoA(nr-1,nz-1,f1,xx)

  !TODO: here just simply modify nr, may not correct
  ! instead of find half away pre-conditioner, direct use regular grids
  call DCT(nr-1,nz-1,f1,f2,Ph_Lcc,wsavexc_ph,wsaveyc_ph)

  call Atob(nr-1,nz-1,f2,yy)

  end subroutine

