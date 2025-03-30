  !  solve phi and mu together on staggered grids

  subroutine ph_mu_solve(ph_after, rmu_after, t)
  use global
  implicit none
  integer i,k
  real ph_after(-1:nr,-1:nz),rmu_after(-1:nr,-1:nz),t

  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)

  !real pnt(0:nr),pnb(0:nr) !--------------------------------

  real rhs_ph_mu(2*nr*nz),Computed_Solution_ph_mu(2*nr*nz)
  real work_ph_mu(2*nr*nz,7)
  real,external::ch_bulk_energy

  real tmp
  real pnt(0:nr-1),pnb(0:nr-1)
  real,external::ch_bdr_condition


  external MATRIX_VECTOR_PH_MU_CG,PSolve_PH_MU_CG
  integer ITER,INFO
  real RESID

  external MATRIX_VECTOR_PH_MU_fgmres,PSolve_PH_MU_fgmres

!  ! TODO:
!  real gT(0:nr),gB(0:nr),cf

  if (phi_mode == 1) then
    ph_M_0 = (1.0 - bdr_ph) * ph**2.0 * (ph - 1.0)**2.0 ! 1.0 - bdr_ph
  elseif (phi_mode == 2) then 
    ph_M_0 = (1.0 - bdr_ph) * (ph**2.0 - 1.0)**2.0 ! 1.0 - bdr_ph
  endif


  !===============================================================================
  ! right hand side for phi and rmu
  f2 = ph * (bdr_ph_used + bdr_ph_used_adjust)

  !TODO: should be u\grad\phi
  do k=0,nz-1
    do i=0,nr-1
      f2(i,k) = f2(i,k) - dt*(u(i+1,k)+u(i,k))/2.0 * (ph(i+1,k)*bdr_ph_used(i+1,k) - ph(i-1,k)*bdr_ph_used(i-1,k) )/2.0/dr&
      & - dt*(v(i,k+1)+v(i,k))/2.0 * (ph(i,k+1)*bdr_ph_used(i,k+1) - ph(i,k-1)*bdr_ph_used(i,k-1) )/2.0/dz

    !   f2(i,k) = f2(i,k) - dt*(u(i+1,k)+u(i,k))/2.0 * (ph(i+1,k)-ph(i-1,k))/2.0/dr&
    !     & - dt*(v(i,k+1)+v(i,k))/2.0 * (ph(i,k+1)-ph(i,k-1))/2.0/dz
    enddo
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (phi_mode == 1) then 
    do k=0,nz-1
      do i=0,nr-1
          ! f1(i,k) = s_implicit*ph(i,k)/I_thickness* bdr_ph_used(i,k) - ch_bulk_energy(ph(i,k), bdr_ph(i,k))/I_thickness &
          ! - ph(i,k)*(ph(i,k)-1.0)*bdr_ph_grad(i,k)*cos(angle*pi/180)/sqrt(2.0)

          f1(i,k) = s_implicit*ph(i,k)/I_thickness - ch_bulk_energy(ph(i,k), bdr_ph(i,k), phi_mode)/I_thickness &
          - ph(i,k)*(ph(i,k)-1.0)*bdr_ph_grad(i,k)*cos(angle*pi/180)/sqrt(2.0)
      enddo
    enddo
  elseif (phi_mode == 2) then

    do k=0,nz-1
      do i=0,nr-1
          ! f1(i,k) = s_implicit*ph(i,k)/I_thickness* bdr_ph_used(i,k) - ch_bulk_energy(ph(i,k), bdr_ph(i,k))/I_thickness &
          ! - ph(i,k)*(ph(i,k)-1.0)*bdr_ph_grad(i,k)*cos(angle*pi/180)/sqrt(2.0)

          f1(i,k) = s_implicit*ph(i,k)/I_thickness - ch_bulk_energy(ph(i,k), bdr_ph(i,k), phi_mode)/I_thickness &
          - (ph(i,k)**2.0-1.0)*bdr_ph_grad(i,k)*cos(angle*pi/180)/sqrt(2.0)
      enddo
    enddo

  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do i=0,nr
!
!   cf = sqrt(2.0)*cos(angleB(i)*pi/180.0)/3.0*pi/2.
!
!   pnb(i) = cf*cos(ph(i,0)*pi/2) + alpha_s*ph(i,0)
!
!
!   f1(i,0) = f1(i,0) + pnb(i)*2/dz
!
!
!   cf = sqrt(2.0)*cos(angleT(i)*pi/180.0)/3.0*pi/2.
!
!   pnt(i) = cf*cos(ph(i,nz)*pi/2) + alpha_s*ph(i,nz)
!
!   f1(i,nz) = f1(i,nz) + pnt(i)*2/dz
!enddo


!        call output(nr,nz,r,z,u,v,f1,f2,0.0,'testin5')
!        print *,'output test 5'
!
!        call output(nr,nz,r,z,u,v,ph,rmu,0.0,'testin6')
!        print *,'output test 5'
!
!        pause


  ! bottom bdr
  do i=0,nr-1
    k = 0
    tmp = (ph(i,k-1) + ph(i,k))/2.0

   !  cf = sqrt(2.0)*cos(angleB(i)*pi/180.0)/3.0*pi/2.
   !  pnb(i) = cf*cos(tmp*pi/2) + alpha_s*tmp

    pnb(i) = ch_bdr_condition(angleB(i), tmp, pi, phi_mode) !+ alpha_s*tmp
  enddo

  do i=0,nr-1
    k = 0

    tmp = I_thickness/dz
    coef_bot = (I_thickness/dz) / tmp
    pnb(i) = pnb(i)/tmp

    f1(i,k) = f1(i,k) + pnb(i)*I_thickness/dz2
  enddo

  ! top bdr
  do i=0,nr-1
    k = nz-1
    tmp = (ph(i,k) + ph(i,k+1))/2.0

   !  cf = sqrt(2.0)*cos(angleT(i)*pi/180.0)/3.0*pi/2.
   !  pnt(i) = cf*cos(tmp*pi/2) + alpha_s*tmp

    pnt(i) = ch_bdr_condition(angleT(i), tmp, pi, phi_mode) !+ alpha_s*tmp
  enddo

  do i=0,nr-1
    k = nz-1

    tmp = I_thickness/dz !+ 0.5* alpha_s
    coef_top = (I_thickness/dz ) / tmp

    pnt(i) = pnt(i)/tmp

    f1(i,k) = f1(i,k) + pnt(i)*I_thickness/dz2
  enddo  


    call Atob_ph_mu(nr,nz,f1,f2,rhs_ph_mu)
    call Atob_ph_mu(nr,nz,ph,rmu,Computed_Solution_ph_mu)


ITER = 5000
RESID = tol
    call BICGSTAB(2*nr*nz,rhs_ph_mu,Computed_Solution_ph_mu,work_ph_mu,2*nr*nz,&
      &ITER,RESID,MATRIX_VECTOR_PH_MU_CG,PSolve_PH_MU_CG,INFO)
	write(*,*) 'ph_mu_solve: \phi 1 BICGSTAB',ITER, RESID, tol
	if(INFO .ne. 0) then
		print *,'error', INFO, tol
		!stop
	endif


!RESID = tol
!call fgmres_general(2*(nr+1)*(nz+1),150,Computed_Solution_ph_mu,rhs_ph_mu,&
! &RESID,MATRIX_VECTOR_PH_MU_fgmres,PSolve_PH_MU_fgmres)
!write(*,*) 'fgmres'


    call btoA_ph_mu(nr,nz,ph_after,rmu_after,Computed_Solution_ph_mu)
    do i=0,nr-1
      ph_after(i,-1) = ph_after(i,0) * coef_bot + pnb(i)
      ph_after(i,nz) = ph_after(i,nz-1) * coef_top + pnt(i)
    enddo

!        call output_regular(nr,nz,r,z,ph,rmu,ph_after,rmu_after,0.0,'testin6')
!        print *,'output test 6'
!        pause


! TODO:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! update ghost points for \phi
!  do i=0,nr
!     k=0
!     ph_after(i,k-1)=ph_after(i,k+1) + (pnb(i)-alpha_s*ph_after(i,k))*2*dz/I_thickness - 2*dz*xld/vs*( (rmu_after(i+1,k)&
!& +rmu_after(i-1,k)-2*rmu_after(i,k))/dr2 + (rmu_after(i,k+1)+rmu_after(i,k-1)-2*rmu_after(i,k))/dz2 )/I_thickness
!
!     k=nz
!     ph_after(i,k+1)=ph_after(i,k-1) + (pnt(i)-alpha_s*ph_after(i,k))*2*dz/I_thickness - 2*dz*xld/vs*( (rmu_after(i+1,k)&
!& +rmu_after(i-1,k)-2*rmu_after(i,k))/dr2 + (rmu_after(i,k+1)+rmu_after(i,k-1)-2*rmu_after(i,k))/dz2 )/I_thickness
!  enddo


  end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Matrix operator of linear system for CH equation
  subroutine  MATRIX_VECTOR_PH_MU_CG(alpha,xx,beta,yy)
  use global
  implicit none
  real xx(2*nr*nz), zz(2*nr*nz)
  real yy(2*nr*nz),alpha,beta
  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
  real ff1(-1:nr,-1:nz),ff2(-1:nr,-1:nz)
  integer i,j

  call btoA_ph_mu(nr,nz,f1,f2,xx)

  call op_double_gradient_Neum(nr, nz, dr2, dz2, f1, ph_M_inner, ff1)

  do j=0,nz-1
     do i=0,nr-1
        ! ff1(i,j) = s_implicit*f1(i,j)/I_thickness - I_thickness*( &
        !  &(f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)-2*f1(i,j))/dz2 ) - f2(i,j)

        !ff1(i,j) = s_implicit*f1(i,j)/I_thickness*bdr_ph_used(i,j) - I_thickness*ff1(i,j) - f2(i,j) * bdr_ph_used(i,j)
        ff1(i,j) = s_implicit*f1(i,j)/I_thickness - I_thickness*ff1(i,j) - f2(i,j) * bdr_ph_used(i,j)
      enddo
  enddo

  ! modify operators to satisfy non-homo boundary conditions
  do i=0,nr-1
    j=0
    ff1(i,j) = ff1(i,j) + I_thickness/dz2*( f1(i,j-1) - coef_bot*f1(i,j) )

    j=nz-1
    ff1(i,j) = ff1(i,j) + I_thickness/dz2*( f1(i,j+1) - coef_top*f1(i,j) )
 enddo


!  do i=0,nr
!     j=0
!     ff1(i,j) = ff1(i,j) + 2*alpha_s*f1(i,j)/dz + 2*xld/vs/dz*(&
!       &(f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 + (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
!     j=nz
!     ff1(i,j) = ff1(i,j) + 2*alpha_s*f1(i,j)/dz + 2*xld/vs/dz*(&
!       &(f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 + (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
!  enddo

  call op_double_gradient_Neum(nr, nz, dr2, dz2, f2, ph_M_0, ff2)

  do j=0,nz-1
     do i=0,nr-1
        ! ff2(i,j) = f1(i,j) - dt*xld*( (f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 + (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
        ff2(i,j) = f1(i,j) * (bdr_ph_used(i,j) + bdr_ph_used_adjust) &
            & - dt*xld*ff2(i,j)
      enddo
  enddo



!  do j=0,nz
!     do i=0,nr
!        ff1(i,j) = f1(i,j) - (s_implicit*f2(i,j)/I_thickness - I_thickness*( &
!         &(f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 + (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 ))
!     enddo
!  enddo
!
!!  do i=0,nr
!!     j=0
!!     ff1(i,j) = ff1(i,j) + 2*alpha_s*f1(i,j)/dz + 2*xld/vs/dz*(&
!!       &(f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 + (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
!!     j=nz
!!     ff1(i,j) = ff1(i,j) + 2*alpha_s*f1(i,j)/dz + 2*xld/vs/dz*(&
!!       &(f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 + (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
!!  enddo
!
!
!  do j=0,nz
!     do i=0,nr
!        ff2(i,j) = f2(i,j) - dt*xld*( (f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)-2*f1(i,j))/dz2 )
!      enddo
!  enddo



  !TODO: xx cannot be changed
    call Atob_ph_mu(nr,nz,ff1, ff2, zz)
    yy = alpha*zz+beta*yy


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

  subroutine PSolve_PH_MU_CG(yy,xx)
  use global
  implicit none
  real xx(2*nr*nz),yy(2*nr*nz)
  real f1(-1:nr,-1:nz),f2(-1:nr,-1:nz)
  real ff1(-1:nr,-1:nz),ff2(-1:nr,-1:nz)

  real sol(-1:nr,-1:nz)
  integer i,j,k
  real tmpx, tmpz

!  yy = xx
!  return

  ! Jacobi preconditioner
  call btoA_ph_mu(nr,nz,f1,f2,xx)

  do j=0,nz-1
     do i=0,nr-1
          if ( i == 0 .or. j == nr-1 ) then
            tmpx = 1.0/dr2
          else
            tmpx = 2.0/dr2
          endif
          if ( i == 0 .or. j == nz-1 ) then
            tmpz = 1.0/dz2
          else
            tmpz = 2.0/dz2
          endif

          ff1(i,j) = f1(i,j) / ( (tmpx+tmpz)*I_thickness + s_implicit/I_thickness )
          !ff1(i,j) = f1(i,j) / ( (tmpx+tmpz)*I_thickness + s_implicit/I_thickness* bdr_ph_used(i,j) )
     enddo
  enddo

  do j=0,nz-1
     do i=0,nr-1
          if ( i == 0 .or. j == nr-1 ) then
            tmpx = 1.0/dr2
          else
            tmpx = 2.0/dr2
          endif
          if ( i == 0 .or. j == nz-1 ) then
            tmpz = 1.0/dz2
          else
            tmpz = 2.0/dz2
          endif
          ff2(i,j) = f2(i,j) / ( dt*xld*(tmpx+tmpz) )
      enddo
  enddo

  call Atob_ph_mu(nr,nz,ff1, ff2, yy)



! FFT precondioner
!  call btoA_ph_mu(nr,nz,f1,f2,xx)
!
!    do j=0,nz
!       do i=0,nr
!           f1(i,j) = -f1(i,j) + s_implicit*f2(i,j)/I_thickness - I_thickness*( (f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 &
!              &+ (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
!       enddo
!    enddo
!
!  call DCT(nr,nz,f1,sol,Ph_Lcc,wsavexc,wsaveyc)
!
!  do k=0,nz
!     do i=0,nr
!         f1(i,k) =  f2(i,k)  + dt*xld*( (sol(i+1,k)+sol(i-1,k)-2*sol(i,k))/dr2 &
!          &+ (sol(i,k+1)+sol(i,k-1)-2*sol(i,k))/dz2 )
!     enddo
!  enddo
!
!  call Atob_ph_mu(nr,nz,f1,sol, yy)

  end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Matrix operator of linear system for CH equation
  subroutine  MATRIX_VECTOR_PH_MU_fgmres(xx,yy)
  use global
  implicit none
  real xx(2*(nr+1)*(nz+1))
  real yy(2*(nr+1)*(nz+1))
  real f1(-1:nr+1,-1:nz+1),f2(-1:nr+1,-1:nz+1)
  real ff1(-1:nr+1,-1:nz+1),ff2(-1:nr+1,-1:nz+1)
  integer i,j

  call btoA_ph_mu(nr,nz,f1,f2,xx)

  do j=0,nz
     do i=0,nr
        ff1(i,j) = s_implicit*f1(i,j)/I_thickness - I_thickness*( &
         &(f1(i+1,j)+f1(i-1,j)-2*f1(i,j))/dr2 + (f1(i,j+1)+f1(i,j-1)-2*f1(i,j))/dz2 ) - f2(i,j)
     enddo
  enddo

!  do i=0,nr
!     j=0
!     ff1(i,j) = ff1(i,j) + 2*alpha_s*f1(i,j)/dz + 2*xld/vs/dz*(&
!       &(f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 + (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
!     j=nz
!     ff1(i,j) = ff1(i,j) + 2*alpha_s*f1(i,j)/dz + 2*xld/vs/dz*(&
!       &(f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 + (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
!  enddo


  do j=0,nz
     do i=0,nr
        ff2(i,j) = f1(i,j) - dt*xld*( (f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 + (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
      enddo
  enddo

    call Atob_ph_mu(nr,nz,ff1, ff2, yy)


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

  subroutine PSolve_PH_MU_fgmres(xx, yy)
  use global
  implicit none
  real xx(2*(nr+1)*(nz+1)),yy(2*(nr+1)*(nz+1))
  real f1(-1:nr+1,-1:nz+1),f2(-1:nr+1,-1:nz+1)
  real sol(-1:nr+1,-1:nz+1)
  integer i,j,k

  call btoA_ph_mu(nr,nz,f1,f2,xx)

    do j=0,nz
       do i=0,nr
           f1(i,j) = -f1(i,j) + s_implicit*f2(i,j)/I_thickness - I_thickness*( (f2(i+1,j)+f2(i-1,j)-2*f2(i,j))/dr2 &
              &+ (f2(i,j+1)+f2(i,j-1)-2*f2(i,j))/dz2 )
       enddo
    enddo

  call DCT(nr,nz,f1,sol,Ph_Lcc,wsavexc,wsaveyc)

  do k=0,nz
     do i=0,nr
         f1(i,k) =  f2(i,k)  + dt*xld*( (sol(i+1,k)+sol(i-1,k)-2*sol(i,k))/dr2 &
          &+ (sol(i,k+1)+sol(i,k-1)-2*sol(i,k))/dz2 )
     enddo
  enddo

  call Atob_ph_mu(nr,nz,f1,f2, yy)

!  yy = xx

!  call btoA(nr,nz,f1,xx)
!
!  call DCT(nr,nz,f1,f2,Ph_Lcc,wsavexc,wsaveyc)
!
!  call Atob(nr,nz,f2,yy)

  end subroutine
